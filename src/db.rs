// SPDX-License-Identifier: Apache-2.0
// Rust adaptation of HaploGrouper tree, locus and weight handling.
use anyhow::{Context, Result, ensure};
use std::{collections::HashMap, fs, path::Path};

pub struct Regions(Vec<(usize, usize)>);
impl Regions {
    pub fn parse(s: Option<&str>) -> Result<Self> {
        let mut intervals = Vec::new();
        if let Some(s) = s {
            for interval in s.split(',') {
                let (a, b) = interval
                    .split_once('-')
                    .context("expected start-end for each region")?;
                let (a, b) = (a.parse()?, b.parse()?);
                ensure!(a <= b, "invalid region {interval}");
                intervals.push((a, b));
            }
        }
        Ok(Self(intervals))
    }
    pub fn contains(&self, pos: usize) -> bool {
        self.0.is_empty() || self.0.iter().any(|&(a, b)| a <= pos && pos <= b)
    }
}

pub struct Node {
    pub name: String,
    pub parent: Option<usize>,
    pub mutations: Vec<usize>,
}
pub struct Mutation {
    pub node: usize,
    pub pos: usize,
    pub derived: String,
    pub ancestral: String,
    pub weight: f64,
}
pub struct Tree {
    pub nodes: Vec<Node>,
    pub order: Vec<usize>, // parents before children; reporting still uses file order
    pub mutations: Vec<Mutation>,
    pub by_pos: HashMap<usize, Vec<usize>>,
    pub positions: Vec<usize>, // first occurrence order in the locus table
}

impl Tree {
    pub fn load(
        tree: &Path,
        loci: &Path,
        weights: Option<&Path>,
        regions: &Regions,
    ) -> Result<Self> {
        let text = fs::read_to_string(tree)?;
        let mut nodes = Vec::new();
        let mut names = HashMap::new();
        let mut parents = Vec::new();
        for (line, row) in text.lines().enumerate() {
            let mut cols = row.split('\t');
            let name = cols.next().unwrap();
            let parent = cols.next().with_context(|| {
                format!("{}:{}: missing parent column", tree.display(), line + 1)
            })?;
            ensure!(
                !name.is_empty() && names.insert(name, line).is_none(),
                "empty or duplicate node name at tree line {}",
                line + 1
            );
            nodes.push(Node {
                name: name.into(),
                parent: None,
                mutations: vec![],
            });
            parents.push(parent);
        }
        let mut children = vec![vec![]; nodes.len()];
        let mut order = Vec::new();
        for (i, parent) in parents.into_iter().enumerate() {
            if let Some(&p) = names.get(parent) {
                nodes[i].parent = Some(p);
                children[p].push(i);
            } else {
                // Upstream trees use a parent absent from the node table for the root.
                order.push(i);
            }
        }
        ensure!(order.len() == 1, "expected one root in the tree");
        let mut i = 0;
        while i < order.len() {
            order.extend_from_slice(&children[order[i]]);
            i += 1;
        }
        ensure!(order.len() == nodes.len(), "cycle in tree");
        let mut db = Self {
            nodes,
            order,
            mutations: vec![],
            by_pos: HashMap::new(),
            positions: vec![],
        };
        let mut skipped = 0;
        for (line, row) in fs::read_to_string(loci)?.lines().enumerate().skip(1) {
            let cols: Vec<_> = row.split('\t').collect();
            ensure!(cols.len() >= 5, "short locus record at line {}", line + 1);
            let (Ok(pos), Some(&node)) = (cols[1].parse::<usize>(), names.get(cols[3])) else {
                skipped += 1;
                continue;
            };
            ensure!(pos > 0, "locus at position zero on line {}", line + 1);
            if !regions.contains(pos) {
                continue;
            }
            let idx = db.mutations.len();
            let site = db.by_pos.entry(pos).or_default();
            if site.is_empty() {
                db.positions.push(pos);
            }
            site.push(idx);
            db.nodes[node].mutations.push(idx);
            db.mutations.push(Mutation {
                node,
                pos,
                derived: cols[2].to_ascii_uppercase(),
                ancestral: cols[4].to_ascii_uppercase(),
                weight: 1.0,
            });
        }
        if skipped > 0 {
            eprintln!("Skipped {skipped} loci with a nonnumeric position or unknown branch.");
        }
        if let Some(path) = weights {
            for (line, row) in fs::read_to_string(path)?.lines().enumerate() {
                let cols: Vec<_> = row.split('\t').collect();
                ensure!(cols.len() >= 4, "short weight record at line {}", line + 1);
                let pos: usize = cols[0].parse()?;
                let weight: f64 = cols[3].parse()?;
                ensure!(weight.is_finite(), "nonfinite weight at line {}", line + 1);
                if let Some(indices) = db.by_pos.get(&pos) {
                    for &idx in indices {
                        let m = &mut db.mutations[idx];
                        if m.ancestral == cols[1] && m.derived == cols[2] {
                            m.weight = weight;
                        }
                    }
                }
            }
        }
        Ok(db)
    }
}
