# HaploGrouperRS

HaploGrouperRS assigns haplogroups to samples in a haploid VCF using a phylogenetic tree and its branch mutations. This Rust port of HaploGrouper reads plain VCF, gzip and BGZF. The repository includes ISOGG Y-chromosome reference tables for GRCh37 and GRCh38. Mitochondrial and custom trees use the same interface.

## Run an existing HaploGrouper command

Replace `python haploGrouper.py` with `haplogrouper-rs` (or `haplogrouper-rs.exe` on Windows). Keep the original arguments, filename case and working directory. Paths to the input and output files are relative to that directory.

~~~sh
haplogrouper-rs -v sample.vcf.gz -t data/chrY_isogg2019_tree.txt -l data/chrY_isogg2019-decode1_loci_b37.txt -c Y -m -o calls.tsv
~~~

Use the chromosome label from your VCF, such as `Y` or `chrY`. For GRCh38 input, change the locus file to `chrY_isogg2019-decode1_loci_b38.txt`. For mtDNA, point the tree and locus options at the upstream PhyloTree Build 17 tables. Upstream example VCFs, the mitochondrial references and callable-region BED files can be obtained from an existing [HaploGrouper](https://gitlab.com/bio_anth_decode/haploGrouper) checkout. [Upstream data](https://gitlab.com/bio_anth_decode/haploGrouper/-/tree/master/data) · [Included Y references](data/Y-DATA.md)

## Build or download

~~~sh
cargo build --release --locked
cargo test --locked
~~~

The executable is `target/release/haplogrouper-rs` on Unix or `target/release/haplogrouper-rs.exe` on Windows. Run it with `--help` to inspect the options.

The **Build HaploGrouperRS** GitHub Actions workflow runs on pushes, pull requests and manual dispatch. Open a successful run and download the artifact for your target. Extract the platform archive and run the executable from the extracted directory to use the bundled `data/` paths. Each archive includes the Y references, attribution, software license and dependency notices.

| Platform | Artifact suffix | Archive |
| --- | --- | --- |
| Windows x64 | `x86_64-pc-windows-msvc` | ZIP |
| Linux x64 (built on Ubuntu 22.04) | `x86_64-unknown-linux-gnu` | tar.gz |
| macOS Apple Silicon | `aarch64-apple-darwin` | tar.gz |
| macOS Intel | `x86_64-apple-darwin` | tar.gz |

## Command-line options

| Short | Upstream long option | Input or output |
| --- | --- | --- |
| `-v` | `--vcfFile` | Required. Coordinate-sorted VCF with a FORMAT column and at least one sample. |
| `-t` | `--hGrpTreeFile` | Required. Tab-separated node and parent names, without a header. The root has no parent in the table. |
| `-l` | `--hGrpLocusFile` | Required. Tab-separated markers with a header. The first five columns are marker name, 1-based position, derived allele, branch node and ancestral allele. |
| `-o` | `--outFile` | Required. Output assignment TSV. |
| `-i` | `--IDListFile` | Sample IDs, one per line, in the first tab-separated field. Output follows VCF sample order. Defaults to all samples. |
| `-r` | `--regions` | Comma-separated inclusive intervals, for example `0-600,16024-16365`. Coordinates also apply to reference inference. |
| `-c` | `--chrom` | Exact chromosome name. Supply this for a multi-chromosome VCF. |
| `-f` | `--referenceFasta` | Use the first FASTA sequence to infer reference alleles at sites absent from the VCF. Restrict inference to sequenced regions. A missing GT in an existing record stays missing. |
| `-w` | `--weightFile` | Four tab-separated columns without a header giving position, ancestral allele, derived allele and weight. Default weight is 1. |
| `-m` | `--mismatchDetails` | Include back mutations, third alleles and unexplained differences. |
| `-x` | `--verboseFile` | Write a second TSV with a row per scored node or ancestor and a column per selected sample. |
| `-h` | `--help` | Print usage. |

## Interpreting the output

Use a tree, marker table and VCF with matching coordinates. For a GT with multiple alleles the program uses the first allele, following the upstream haploid convention. It sums support and opposition along each path from the root. The net score is support minus opposition. A top-scoring candidate must have support on its terminal branch. Ties are resolved by minimum opposition and then by the most recent common ancestor. When no candidate qualifies, the reported haplogroup is the root.

The upstream HG00096 example produces **H16a1** with a net score of **45** for mtDNA, and **R1b1a1b1a1a2c1a1f1** with a net score of **952** for the deCODE GRCh37 Y reference. The Rust port also handles the trailing bases of the final VCF indel and returns a nonzero exit code for invalid inputs.

## Memory use

Uniform sites store one allele for the cohort. Biallelic sites with at least 16 selected samples store one bit per sample, and other sites retain the observed allele bytes. Reference inference reads FASTA lines in sequence. A large `-x` matrix uses a 4 MiB score buffer (or one sample column for very large trees) and an automatically removed file in the system temporary directory. Allow **16 bytes × emitted nodes × selected samples** of temporary disk space when the matrix spills. An assignment-only run needs no matrix scratch file.

Peak process working set measured on Windows x64 using release builds, the same upstream VCFs and all samples:

| Dataset | Samples | Output | Previous port | 0.2.0 |
| --- | ---: | --- | ---: | ---: |
| GRCh37 Y, ISOGG 2019 deCODE | 1,233 | Assignments | 50.36 MiB | 29.12 MiB |
| GRCh37 Y, ISOGG 2019 deCODE | 1,233 | Assignments and score matrix | 119.63 MiB | 32.40 MiB |
| Phase 3 mtDNA, Build 17 with rCRS | 2,534 | Assignments | 16.87 MiB | 10.88 MiB |
| Phase 3 mtDNA, Build 17 with rCRS | 2,534 | Assignments and score matrix | 227.02 MiB | 14.96 MiB |

The four comparisons produced identical SHA-256 digests for the corresponding output files.

## License and citation

Commercial use of the HaploGrouperRS 0.2.0 contributions is prohibited under the project license. This is a source-available, noncommercial distribution. Upstream HaploGrouper and earlier Apache-licensed contributions retain their Apache-2.0 grants. The ISOGG-derived data carries attribution, noncommercial and ShareAlike terms under CC BY-NC-SA 3.0. [Project license](LICENSE) · [Original Apache license](LICENSES/Apache-2.0.txt) · [Y-data attribution](data/Y-DATA.md) · [NOTICE](NOTICE)

For research using this implementation, cite Jagadeesan et al., *HaploGrouper: a generalized approach to haplogroup classification*, Bioinformatics 37, 570–572, and the source of your phylogeny. [HaploGrouper paper](https://doi.org/10.1093/bioinformatics/btaa729) · [PhyloTree citation](https://www.phylotree.org/)
