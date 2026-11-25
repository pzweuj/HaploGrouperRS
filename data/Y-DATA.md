# Y-chromosome reference data

These files are byte-for-byte copies from the HaploGrouper snapshot used for the Rust port. The International Society of Genetic Genealogy (ISOGG) compiled the Y-DNA Haplogroup Tree. The HaploGrouper authors prepared the tree and marker tables. The deCODE fork revised marker assignments within the 2019 tree. [ISOGG tree](https://isogg.org/tree/) · [Upstream HaploGrouper data](https://gitlab.com/bio_anth_decode/haploGrouper/-/tree/master/data)

| Tree file | Matching locus files | Origin |
| --- | --- | --- |
| `chrY_isogg2019_tree.txt` | `chrY_isogg2019_loci_b37.txt`, `chrY_isogg2019_loci_b38.txt` | ISOGG Y-DNA Haplogroup Tree 2019, version 14.196, formatted by HaploGrouper. |
| `chrY_isogg2019_tree.txt` | `chrY_isogg2019-decode1_loci_b37.txt`, `chrY_isogg2019-decode1_loci_b38.txt` | deCODE genetics revised the placement of some markers and removed others. The tree topology and names follow ISOGG 14.196. [Original change notes](chrY_isogg2019-decode1_loci_README.txt) |
| `chrY_isogg2016-deprecated_tree.txt` | `chrY_isogg2016-deprecated_loci_b37.txt` | Historical ISOGG Y-DNA Haplogroup Tree 2016 supplied for reproducing older analyses. |

ISOGG lists the Y-SNP Tree under Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported. Commercial use of the ISOGG-derived material is prohibited under this license. When sharing, credit ISOGG and the named upstream contributors, include the license link, identify changes and retain the notices. Adaptations carry the applicable ShareAlike terms. The upstream HaploGrouper license is also included with this distribution. [ISOGG terms](https://www.isogg.org/wiki/ISOGG_Wiki:Terms_of_Use) · [CC BY-NC-SA 3.0 legal code](https://creativecommons.org/licenses/by-nc-sa/3.0/legalcode) · [Upstream license](../LICENSES/Apache-2.0.txt)

For publications, cite the ISOGG tree year and version, give the tree URL and access date, and identify the deCODE fork if used. Also cite Jagadeesan et al., *HaploGrouper: a generalized approach to haplogroup classification*, Bioinformatics 37, 570–572. [ISOGG citation guidance](https://www.isogg.org/wiki/YTree%3ACan_I_use_ISOGG_YSNP_Tree_data%3F) · [HaploGrouper paper](https://doi.org/10.1093/bioinformatics/btaa729)
