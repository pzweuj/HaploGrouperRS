// The original port was licensed under Apache-2.0. See LICENSE for the 0.2.0 changes.
// Rust adaptation of the HaploGrouper CLI. See NOTICE.
mod db;
mod report;
mod vcf;

use anyhow::{Context, Result, ensure};
use clap::Parser;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use tempfile::NamedTempFile;

#[derive(Parser)]
#[command(
    version,
    infer_long_args = true,
    args_override_self = true,
    about = "Assign haplogroups from haploid genotypes in a VCF"
)]
struct Args {
    /// Input VCF (plain text, gzip or BGZF)
    #[arg(short = 'v', long = "vcfFile")]
    vcf: PathBuf,
    /// Output assignment table (TSV)
    #[arg(short = 'o', long = "outFile")]
    output: PathBuf,
    /// Branch mutations (upstream locus TSV)
    #[arg(short = 'l', long = "hGrpLocusFile")]
    loci: PathBuf,
    /// Phylogeny (node and parent TSV)
    #[arg(short = 't', long = "hGrpTreeFile")]
    tree: PathBuf,
    /// File of sample IDs, one per line (first tab-separated field)
    #[arg(short = 'i', long = "IDListFile")]
    ids: Option<PathBuf>,
    /// Inclusive intervals, e.g. 0-600,16024-16365
    #[arg(short = 'r', long)]
    regions: Option<String>,
    /// Exact VCF chromosome name
    #[arg(short = 'c', long)]
    chrom: Option<String>,
    /// Infer unreported sites from the first sequence in this FASTA
    #[arg(short = 'f', long = "referenceFasta")]
    fasta: Option<PathBuf>,
    /// Mutation weights (position, ancestral allele, derived allele, weight)
    #[arg(short = 'w', long = "weightFile")]
    weights: Option<PathBuf>,
    /// Write a per-node score matrix to this TSV
    #[arg(short = 'x', long = "verboseFile")]
    verbose: Option<PathBuf>,
    /// Include back mutation and third allele details
    #[arg(short = 'm', long = "mismatchDetails")]
    details: bool,
}

fn parent(path: &Path) -> &Path {
    path.parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."))
}

fn output_path(path: &Path) -> Result<PathBuf> {
    if path.exists() {
        ensure!(path.is_file(), "{} is not a file", path.display());
        Ok(path.canonicalize()?)
    } else {
        Ok(parent(path)
            .canonicalize()?
            .join(path.file_name().context("missing filename")?))
    }
}

fn run(mut args: Args) -> Result<()> {
    // Upstream argparse represents omitted optional values as empty strings.
    for option in [
        &mut args.ids,
        &mut args.fasta,
        &mut args.weights,
        &mut args.verbose,
    ] {
        if option.as_deref() == Some(Path::new("")) {
            *option = None;
        }
    }
    args.regions = args.regions.filter(|s| !s.is_empty());
    args.chrom = args.chrom.filter(|s| !s.is_empty());
    let output = output_path(&args.output).with_context(|| args.output.display().to_string())?;
    let verbose = args.verbose.as_deref().map(output_path).transpose()?;
    ensure!(
        Some(&output) != verbose.as_ref(),
        "assignment and score matrix require separate paths"
    );
    for input in [
        Some(&args.vcf),
        Some(&args.tree),
        Some(&args.loci),
        args.ids.as_ref(),
        args.fasta.as_ref(),
        args.weights.as_ref(),
    ]
    .into_iter()
    .flatten()
    {
        let path = input
            .canonicalize()
            .with_context(|| input.display().to_string())?;
        ensure!(
            path != output && Some(&path) != verbose.as_ref(),
            "output would overwrite input {}",
            input.display()
        );
    }
    let regions = db::Regions::parse(args.regions.as_deref())?;
    let tree = db::Tree::load(&args.tree, &args.loci, args.weights.as_deref(), &regions)?;
    let mut cohort = vcf::Cohort::read(&args, &tree, &regions)?;
    if let Some(fasta) = &args.fasta {
        cohort.fill_reference(&tree, fasta)?;
    }
    let mut assignments = NamedTempFile::new_in(parent(&output))?;
    let mut matrix = verbose
        .as_ref()
        .map(|p| NamedTempFile::new_in(parent(p)))
        .transpose()?;
    {
        let mut out = BufWriter::new(assignments.as_file_mut());
        let mut scores = matrix.as_mut().map(|f| BufWriter::new(f.as_file_mut()));
        report::write(&tree, &cohort, args.details, &mut out, scores.as_mut())?;
        out.flush()?;
        if let Some(writer) = scores.as_mut() {
            writer.flush()?;
        }
    }
    if let (Some(tmp), Some(path)) = (matrix, verbose) {
        tmp.persist(path)?;
    }
    assignments.persist(output)?;
    eprintln!("Assigned {} samples to haplogroups.", cohort.ids.len());
    Ok(())
}

fn main() -> Result<()> {
    run(Args::parse())
}
