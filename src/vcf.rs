// The original port was licensed under Apache-2.0. See LICENSE for the 0.2.0 changes.
// Rust adaptation of HaploGrouper genotype extraction and reference inference.
use crate::{
    Args,
    db::{Regions, Tree},
};
use anyhow::{Context, Result, ensure};
use flate2::read::MultiGzDecoder;
use std::{
    collections::{BTreeMap, HashMap, HashSet},
    fs::{self, File},
    io::{BufRead, BufReader},
    path::Path,
};

// A zero byte represents a base deleted by an ALT allele; N is a missing call.
pub enum Calls {
    Genotypes(Vec<u8>),
    Uniform(u8),
    Biallelic { alleles: [u8; 2], bits: Vec<u8> },
}
impl Calls {
    fn from_bases(bases: Vec<u8>) -> Self {
        let first = bases[0];
        let Some(&second) = bases.iter().find(|&&b| b != first) else {
            return Self::Uniform(first);
        };
        if bases.len() < 16 || bases.iter().any(|&b| b != first && b != second) {
            return Self::Genotypes(bases);
        }
        let mut bits = vec![0; bases.len().div_ceil(8)];
        for (i, &base) in bases.iter().enumerate() {
            if base == second {
                bits[i / 8] |= 1 << (i % 8);
            }
        }
        Self::Biallelic {
            alleles: [first, second],
            bits,
        }
    }
    pub fn at(&self, sample: usize) -> u8 {
        match self {
            Self::Genotypes(bases) => bases[sample],
            Self::Uniform(base) => *base,
            Self::Biallelic { alleles, bits } => {
                alleles[((bits[sample / 8] >> (sample % 8)) & 1) as usize]
            }
        }
    }
}
pub struct Site {
    pub pos: usize,
    pub calls: Calls,
}
pub struct Cohort {
    pub ids: Vec<String>,
    pub sites: Vec<Site>,
    pub latest: HashMap<usize, usize>,
}
fn reader(path: &Path) -> Result<Box<dyn BufRead>> {
    let mut input = BufReader::new(File::open(path)?);
    if input.fill_buf()?.starts_with(&[0x1f, 0x8b]) {
        Ok(Box::new(BufReader::new(MultiGzDecoder::new(input))))
    } else {
        Ok(Box::new(input))
    }
}

fn sample_columns(header: &str, ids: Option<&Path>) -> Result<(Vec<usize>, Vec<String>, usize)> {
    let columns: Vec<_> = header.split('\t').collect();
    let fixed = [
        "#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT",
    ];
    ensure!(
        columns.starts_with(&fixed) && columns.len() > fixed.len(),
        "expected a VCF #CHROM header with FORMAT and sample IDs"
    );
    let list = ids.map(fs::read_to_string).transpose()?;
    let wanted = list.as_ref().filter(|s| !s.is_empty()).map(|s| {
        s.lines()
            .filter_map(|l| l.split('\t').next())
            .collect::<HashSet<_>>()
    });
    let (indices, names): (Vec<_>, Vec<_>) = columns[9..]
        .iter()
        .enumerate()
        .filter(|(_, id)| wanted.as_ref().is_none_or(|set| set.contains(**id)))
        .map(|(i, id)| (i, (*id).to_owned()))
        .unzip();
    ensure!(!names.is_empty(), "no samples selected from this VCF");
    Ok((indices, names, columns.len() - 9))
}

impl Cohort {
    pub fn read(args: &Args, tree: &Tree, regions: &Regions) -> Result<Self> {
        let mut stream = reader(&args.vcf)?;
        let mut line = String::new();
        loop {
            line.clear();
            ensure!(stream.read_line(&mut line)? > 0, "missing #CHROM header");
            if !line.starts_with("##") {
                break;
            }
        }
        let (columns, ids, total) =
            sample_columns(line.trim_end_matches(['\r', '\n']), args.ids.as_deref())?;
        let mut cohort = Self {
            ids,
            sites: vec![],
            latest: HashMap::new(),
        };
        let mut pending: BTreeMap<usize, Vec<u8>> = BTreeMap::new();
        let mut current_chrom = String::new();
        loop {
            line.clear();
            if stream.read_line(&mut line)? == 0 {
                break;
            }
            let record = line.trim_end_matches(['\r', '\n']);
            if record.is_empty() {
                continue;
            }
            let mut parts = record.splitn(10, '\t');
            let chrom = parts.next().unwrap();
            if args.chrom.as_deref().is_some_and(|c| c != chrom) {
                continue;
            }
            if current_chrom.is_empty() {
                current_chrom = chrom.into();
            }
            ensure!(
                current_chrom == chrom,
                "VCF has multiple chromosomes; select one with -c"
            );
            let start: usize = parts
                .next()
                .context("missing POS")?
                .parse()
                .context("invalid VCF POS")?;
            let raw_ref = parts.nth(1).context("missing REF")?;
            let reference = if raw_ref == "." { "N" } else { raw_ref };
            ensure!(
                start > 0 && !reference.is_empty() && reference.is_ascii(),
                "invalid VCF POS or REF at {chrom}:{start}"
            );
            let end = start
                .checked_add(reference.len() - 1)
                .context("VCF position overflow")?;
            if !regions.contains(start) || !regions.contains(end) {
                continue;
            }
            let overlaps: Vec<_> = (start..=end)
                .filter(|p| tree.by_pos.contains_key(p))
                .collect();
            if !overlaps.is_empty() {
                let alt = parts.next().context("missing ALT")?;
                let format = parts.nth(3).context("missing FORMAT")?;
                let gt = format
                    .split(':')
                    .position(|x| x == "GT")
                    .with_context(|| format!("missing GT in FORMAT at {chrom}:{start}"))?;
                let alleles: Vec<_> = std::iter::once(reference).chain(alt.split(',')).collect();
                let fields: Vec<_> = parts
                    .next()
                    .context("missing VCF genotypes")?
                    .split('\t')
                    .collect();
                ensure!(
                    fields.len() == total,
                    "wrong sample count at {chrom}:{start}"
                );
                let mut genotypes = Vec::with_capacity(columns.len());
                for &col in &columns {
                    let first = fields[col]
                        .split(':')
                        .nth(gt)
                        .unwrap_or(".")
                        .split(['|', '/'])
                        .next()
                        .unwrap();
                    let index = if first == "." {
                        None
                    } else {
                        let idx: usize = first
                            .parse()
                            .with_context(|| format!("invalid GT at {chrom}:{start}"))?;
                        ensure!(
                            idx < alleles.len(),
                            "allele index out of range at {chrom}:{start}"
                        );
                        Some(idx)
                    };
                    genotypes.push(index);
                }
                for pos in overlaps {
                    let offset = pos - start;
                    let existing = pending.contains_key(&pos);
                    let calls = pending
                        .entry(pos)
                        .or_insert_with(|| vec![b'N'; columns.len()]);
                    for (call, idx) in calls.iter_mut().zip(&genotypes) {
                        if let Some(i) = *idx {
                            let allele = alleles[i];
                            if existing && (i == 0 || allele == "<NON_REF>") {
                                continue;
                            }
                            *call = if allele == "." || allele.starts_with('<') || allele == "*" {
                                b'N'
                            } else {
                                allele.as_bytes().get(offset).copied().unwrap_or(0)
                            };
                        }
                    }
                }
            }
            // Downstream bases of an indel may be updated by overlapping VCF records.
            while pending
                .first_key_value()
                .is_some_and(|(&pos, _)| pos <= start)
            {
                let (pos, calls) = pending.pop_first().unwrap();
                cohort.add(pos, Calls::from_bases(calls));
            }
        }
        for (pos, calls) in pending {
            cohort.add(pos, Calls::from_bases(calls));
        }
        ensure!(
            !cohort.sites.is_empty(),
            "VCF and selected reference loci have no overlapping positions"
        );
        Ok(cohort)
    }

    fn add(&mut self, pos: usize, calls: Calls) {
        self.latest.insert(pos, self.sites.len());
        self.sites.push(Site { pos, calls });
    }
    pub fn fill_reference(&mut self, tree: &Tree, fasta: &Path) -> Result<()> {
        let missing: Vec<_> = tree
            .positions
            .iter()
            .copied()
            .filter(|p| !self.latest.contains_key(p))
            .collect();
        let mut order: Vec<_> = (0..missing.len()).collect();
        order.sort_unstable_by_key(|&i| missing[i]);
        let mut bases = vec![0; missing.len()];
        let mut input = BufReader::new(File::open(fasta)?);
        let mut line = String::new();
        let mut header = false;
        let mut length = 0;
        let mut next = order.into_iter().peekable();
        loop {
            line.clear();
            if input.read_line(&mut line)? == 0 {
                break;
            }
            let seq = line.trim();
            if seq.is_empty() {
                continue;
            }
            if seq.starts_with('>') {
                if header {
                    break;
                }
                header = true;
                continue;
            }
            ensure!(header && seq.is_ascii(), "invalid first FASTA sequence");
            let end = length + seq.len();
            while let Some(&i) = next.peek() {
                if missing[i] > end {
                    break;
                }
                bases[i] = seq.as_bytes()[missing[i] - length - 1].to_ascii_uppercase();
                next.next();
            }
            length = end;
        }
        ensure!(length > 0, "missing first FASTA sequence");
        if let Some(i) = next.next() {
            anyhow::bail!(
                "locus {} exceeds FASTA sequence length {length}",
                missing[i]
            );
        }
        // Keep the locus-file addition order for floating-point score sums.
        for (pos, base) in missing.into_iter().zip(bases) {
            self.add(pos, Calls::Uniform(base));
        }
        Ok(())
    }
}
