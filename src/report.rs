// The original port was licensed under Apache-2.0. See LICENSE for the 0.2.0 changes.
// Rust adaptation of HaploGrouper scoring, tie resolution and TSV reports.
use crate::{db::Tree, vcf::Cohort};
use anyhow::Result;
use std::{
    collections::HashSet,
    io::{BufReader, Read, Seek, SeekFrom, Write},
};

#[derive(Clone, Copy, Default)]
struct Score {
    support: f64,
    against: f64,
}
impl Score {
    fn net(self) -> f64 {
        self.support - self.against
    }
    fn pair(self) -> String {
        format!("{}-{}", number(self.support), number(self.against))
    }
}

// Match the five significant digits used in the upstream TSVs (%.5g).
fn number(x: f64) -> String {
    if !x.is_finite() {
        return x.to_string();
    }
    if x.fract() == 0.0 && x.abs() < 100000.0 {
        return format!("{x:.0}");
    }
    let scientific = format!("{x:.4e}");
    let (mantissa, exp) = scientific.split_once('e').unwrap();
    let exp: i32 = exp.parse().unwrap();
    if (-4..5).contains(&exp) {
        let decimals = (4 - exp) as usize;
        let s = format!("{x:.decimals$}");
        if decimals == 0 {
            s
        } else {
            s.trim_end_matches('0').trim_end_matches('.').into()
        }
    } else {
        format!(
            "{}e{exp:+03}",
            mantissa.trim_end_matches('0').trim_end_matches('.')
        )
    }
}

fn allele(base: &u8) -> &[u8] {
    if *base == 0 {
        &[]
    } else {
        std::slice::from_ref(base)
    }
}

struct SampleScores {
    direct: Vec<Score>,
    path: Vec<Score>,
}
impl SampleScores {
    fn calculate(tree: &Tree, cohort: &Cohort, sample: usize, active: &[bool]) -> Self {
        let mut direct = vec![Score::default(); tree.nodes.len()];
        for site in &cohort.sites {
            let base = site.calls.at(sample);
            let gt = allele(&base);
            for &idx in &tree.by_pos[&site.pos] {
                let m = &tree.mutations[idx];
                let score = &mut direct[m.node];
                if gt == m.derived.as_bytes() {
                    score.support += m.weight;
                }
                if gt == m.ancestral.as_bytes() {
                    score.against += m.weight;
                }
            }
        }
        let mut path = direct.clone();
        for &i in &tree.order {
            if active[i]
                && let Some(p) = tree.nodes[i].parent
            {
                path[i].support = path[p].support + direct[i].support;
                path[i].against = path[p].against + direct[i].against;
            }
        }
        Self { direct, path }
    }
    fn leaders(&self, metric: fn(Score) -> f64) -> Vec<usize> {
        let best = self
            .path
            .iter()
            .map(|s| metric(*s))
            .fold(f64::NEG_INFINITY, f64::max);
        self.path
            .iter()
            .enumerate()
            .filter(|(i, s)| metric(**s) == best && self.direct[*i].support > 0.0)
            .map(|(i, _)| i)
            .collect()
    }
    fn describe(&self, tree: &Tree, nodes: &[usize]) -> String {
        nodes
            .iter()
            .map(|&i| format!("{}[{}]", tree.nodes[i].name, self.path[i].pair()))
            .collect::<Vec<_>>()
            .join(":")
    }
}

fn ancestors(tree: &Tree, node: usize) -> impl Iterator<Item = usize> + '_ {
    std::iter::successors(Some(node), |&i| tree.nodes[i].parent)
}
fn resolve(tree: &Tree, leaders: &[usize], path: &[Score]) -> (usize, String) {
    if leaders.is_empty() {
        return (tree.order[0], "NA".into());
    }
    let min = leaders
        .iter()
        .map(|&i| path[i].against)
        .fold(f64::INFINITY, f64::min);
    let tied: Vec<_> = leaders
        .iter()
        .copied()
        .filter(|&i| path[i].against == min)
        .collect();
    let mut pruning = String::new();
    if tied.len() < leaders.len() {
        pruning.push('A');
    }
    let mut chosen = tied[0];
    if tied.len() > 1 {
        pruning.push('M');
        let mut shared: HashSet<_> = ancestors(tree, chosen).collect();
        for &i in &tied[1..] {
            let chain: HashSet<_> = ancestors(tree, i).collect();
            shared.retain(|n| chain.contains(n));
        }
        chosen = ancestors(tree, chosen)
            .find(|i| shared.contains(i))
            .unwrap();
    }
    if pruning.is_empty() {
        pruning.push_str("NA");
    }
    (chosen, pruning)
}
fn list(items: &[String]) -> String {
    if items.is_empty() {
        "NA".into()
    } else {
        items.join(":")
    }
}
fn assignment(
    tree: &Tree,
    cohort: &Cohort,
    i: usize,
    scores: &SampleScores,
    details: bool,
    out: &mut impl Write,
) -> Result<()> {
    let net = scores.leaders(Score::net);
    let support = scores.leaders(|s| s.support);
    let (chosen, pruning) = resolve(tree, &net, &scores.path);
    let mut mismatch = Vec::new();
    let mut back = Vec::new();
    let mut third = Vec::new();
    let mut seen = HashSet::new();
    for n in ancestors(tree, chosen) {
        for &idx in &tree.nodes[n].mutations {
            let m = &tree.mutations[idx];
            let Some(&site) = cohort.latest.get(&m.pos) else {
                continue;
            };
            let base = cohort.sites[site].calls.at(i);
            let gt = allele(&base);
            let first = seen.insert(m.pos);
            if gt == m.derived.as_bytes() {
                continue;
            }
            let label = format!(
                "{}[{}]{}{}{}",
                String::from_utf8_lossy(gt),
                tree.nodes[n].name,
                m.ancestral,
                m.pos,
                m.derived
            );
            if gt != m.ancestral.as_bytes() {
                third.push(label);
            } else if first {
                mismatch.push(label);
            } else {
                back.push(label);
            }
        }
    }
    let s = scores.path[chosen];
    write!(
        out,
        "{}\t{}\t{}\t{}\t{}\t{}\t{}",
        cohort.ids[i],
        tree.nodes[chosen].name,
        number(s.net()),
        number(s.support),
        number(s.against),
        list(&mismatch),
        back.len()
    )?;
    if details {
        let other = (s.against - mismatch.len() as f64 - back.len() as f64).trunc() as i64;
        write!(
            out,
            "\t{}\t{}\t{}\t{}",
            list(&back),
            third.len(),
            list(&third),
            other
        )?;
    } else {
        write!(out, "\t{}", third.len())?;
    }
    writeln!(
        out,
        "\t{}\t{}\t{}",
        pruning,
        scores.describe(tree, &net),
        scores.describe(tree, &support)
    )?;
    Ok(())
}

pub fn write(
    tree: &Tree,
    cohort: &Cohort,
    details: bool,
    out: &mut impl Write,
    matrix: Option<&mut impl Write>,
) -> Result<()> {
    let mut active = vec![false; tree.nodes.len()];
    for site in &cohort.sites {
        for &idx in &tree.by_pos[&site.pos] {
            active[tree.mutations[idx].node] = true;
        }
    }
    for &i in tree.order.iter().rev() {
        if active[i]
            && let Some(p) = tree.nodes[i].parent
        {
            active[p] = true;
        }
    }
    let nodes: Vec<_> = (0..active.len()).filter(|&i| active[i]).collect();
    write!(
        out,
        "ID\tHaplogroup\tnetScore\tmatchScore\tmismatchScore\tmismatchLoci\tbackMutCnt"
    )?;
    if details {
        write!(
            out,
            "\tbackMutLoci\tthirdAlleleCnt\tthirdAlleleLoci\totherDiff"
        )?;
    } else {
        write!(out, "\tthirdAlleleCnt")?;
    }
    writeln!(out, "\tpruning\tallMaxNetScore\tallMaxScore")?;
    let mut score_sample = |i| -> Result<_> {
        let scores = SampleScores::calculate(tree, cohort, i, &active);
        assignment(tree, cohort, i, &scores, details, out)?;
        Ok(scores.path)
    };
    if let Some(matrix) = matrix {
        write_matrix(
            tree,
            &cohort.ids,
            &nodes,
            matrix,
            score_sample,
            4 * 1024 * 1024,
        )?;
    } else {
        for i in 0..cohort.ids.len() {
            score_sample(i)?;
        }
    }
    Ok(())
}

// Calculate one sample at a time, then emit rows in tree-file order. Large
// matrices are transposed into a temporary file in contiguous row fragments.
fn write_matrix(
    tree: &Tree,
    ids: &[String],
    nodes: &[usize],
    matrix: &mut impl Write,
    mut score_sample: impl FnMut(usize) -> Result<Vec<Score>>,
    budget: usize,
) -> Result<()> {
    const CELL: usize = 16; // support and against, two lossless f64 values
    let width = (budget / CELL / nodes.len()).clamp(1, ids.len());
    let mut block = vec![0u8; nodes.len() * width * CELL];
    let mut spool = if width < ids.len() {
        Some(tempfile::tempfile()?)
    } else {
        None
    };
    for start in (0..ids.len()).step_by(width) {
        let count = width.min(ids.len() - start);
        for col in 0..count {
            let scores = score_sample(start + col)?;
            for (&node, row) in nodes.iter().zip(block.chunks_exact_mut(width * CELL)) {
                let cell = &mut row[col * CELL..(col + 1) * CELL];
                cell[..8].copy_from_slice(&scores[node].support.to_le_bytes());
                cell[8..].copy_from_slice(&scores[node].against.to_le_bytes());
            }
        }
        if let Some(file) = &mut spool {
            for (row, fragment) in block.chunks_exact(width * CELL).enumerate() {
                let offset = (row as u64 * ids.len() as u64 + start as u64) * CELL as u64;
                file.seek(SeekFrom::Start(offset))?;
                file.write_all(&fragment[..count * CELL])?;
            }
        }
    }
    if let Some(file) = &mut spool {
        file.rewind()?;
    }
    let mut reader = spool.map(BufReader::new);
    write!(matrix, "hGrp\ttermMuts")?;
    for id in ids {
        write!(matrix, "\t{id}")?;
    }
    writeln!(matrix)?;
    for (&node, row) in nodes.iter().zip(block.chunks_exact(width * CELL)) {
        let n = &tree.nodes[node];
        write!(matrix, "{}\t{}", n.name, n.mutations.len())?;
        let mut memory = row;
        let source: &mut dyn Read = match &mut reader {
            Some(file) => file,
            None => &mut memory,
        };
        let mut previous = None;
        let mut pair = String::new();
        for _ in ids {
            let mut cell = [0; CELL];
            source.read_exact(&mut cell)?;
            if previous != Some(cell) {
                pair = Score {
                    support: f64::from_le_bytes(cell[..8].try_into().unwrap()),
                    against: f64::from_le_bytes(cell[8..].try_into().unwrap()),
                }
                .pair();
                previous = Some(cell);
            }
            write!(matrix, "\t{pair}")?;
        }
        writeln!(matrix)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::db::Node;
    use std::collections::HashMap;

    #[test]
    fn matrix_with_partial_last_disk_block_matches_in_memory() {
        let nodes = [0, 1, 2];
        let tree = Tree {
            nodes: nodes
                .iter()
                .map(|i| Node {
                    name: format!("N{i}"),
                    parent: None,
                    mutations: vec![],
                })
                .collect(),
            order: nodes.to_vec(),
            mutations: vec![],
            by_pos: HashMap::new(),
            positions: vec![],
        };
        let ids: Vec<_> = (0..5).map(|i| format!("S{i}")).collect();
        let render = |budget| {
            let mut output = Vec::new();
            write_matrix(
                &tree,
                &ids,
                &nodes,
                &mut output,
                |sample| {
                    Ok(nodes
                        .iter()
                        .map(|&n| Score {
                            support: sample as f64 / 2.0,
                            against: n as f64,
                        })
                        .collect())
                },
                budget,
            )
            .unwrap();
            output
        };
        let in_memory = render(16 * 3 * 5);
        assert_eq!(render(16 * 3 * 2), in_memory);
        assert_eq!(render(1), in_memory);
    }
}
