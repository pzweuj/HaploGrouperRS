use flate2::{Compression, write::GzEncoder};
use std::{
    fs,
    io::Write,
    path::Path,
    process::{Command, Output},
};
use tempfile::{TempDir, tempdir};

const HEADER: &str = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\n";

fn dataset(tree: &str, loci: &str, records: &str) -> TempDir {
    let dir = tempdir().unwrap();
    for (name, text) in [
        ("tree", tree.to_owned()),
        ("loci", format!("name\tpos\tder\tnode\tanc\n{loci}")),
        ("input", format!("{HEADER}{records}")),
    ] {
        fs::write(dir.path().join(name), text).unwrap();
    }
    dir
}
fn run(dir: &Path, args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_haplogrouper-rs"))
        .current_dir(dir)
        .args(["-t", "tree", "-l", "loci", "-v", "input", "-o", "calls"])
        .args(args)
        .output()
        .unwrap()
}
fn success(output: Output) {
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}

#[test]
fn tied_siblings_report_their_actual_common_ancestor() {
    let dir = dataset(
        "left\tclade\nroot\t\nright\tclade\nclade\troot\n",
        "l\t1\tT\tleft\tC\nr\t2\tA\tright\tG\n",
        "Y\t1\t.\tC\tT\t.\t.\t.\tGT\t1|0\t.\nY\t2\t.\tG\tA\t.\t.\t.\tGT\t1\t.\n",
    );
    success(run(dir.path(), &["-x", "matrix"]));
    let text = fs::read_to_string(dir.path().join("calls")).unwrap();
    assert_eq!(
        text.lines().nth(1).unwrap(),
        "S1\tclade\t0\t0\t0\tNA\t0\t0\tM\tleft[1-0]:right[1-0]\tleft[1-0]:right[1-0]"
    );
    assert_eq!(
        text.lines().nth(2).unwrap(),
        "S2\troot\t0\t0\t0\tNA\t0\t0\tNA\t\t"
    );
    assert_eq!(
        fs::read_to_string(dir.path().join("matrix")).unwrap(),
        "hGrp\ttermMuts\tS1\tS2\nleft\t1\t1-0\t0-0\nroot\t0\t0-0\t0-0\nright\t1\t1-0\t0-0\nclade\t0\t0-0\t0-0\n"
    );
}

#[test]
fn final_deletion_in_concatenated_gzip_and_absent_reference_site() {
    let dir = dataset(
        "root\t#\nA\troot\n",
        "del\t3\t\tA\tG\nsnp\t4\tT\tA\tC\n",
        "Y\t2\t.\tAG\tA\t.\t.\t.\tDP:GT\t5:1/0\t4:0\n",
    );
    let p = dir.path();
    fs::write(p.join("ids"), "S1\n").unwrap();
    fs::write(p.join("ref"), ">Y\nAAGT\n").unwrap();
    fs::write(p.join("weights"), "3\tG\t\t2.5\n").unwrap();
    let plain = fs::read(p.join("input")).unwrap();
    let mut compressed = Vec::new();
    for part in plain.chunks(30) {
        let mut gz = GzEncoder::new(&mut compressed, Compression::fast());
        gz.write_all(part).unwrap();
        gz.finish().unwrap();
    }
    fs::write(p.join("input"), compressed).unwrap();
    success(run(
        p,
        &[
            "-m", "-f", "ref", "-w", "weights", "-i", "ids", "-c", "Y", "-r", "2-4",
        ],
    ));
    let text = fs::read_to_string(p.join("calls")).unwrap();
    assert_eq!(text.lines().count(), 2);
    assert_eq!(
        text.lines().last().unwrap(),
        "S1\tA\t3.5\t3.5\t0\tNA\t0\tNA\t0\tNA\t0\tNA\tA[3.5-0]\tA[3.5-0]"
    );
}

#[test]
fn malformed_genotype_does_not_replace_an_existing_report() {
    let dir = dataset(
        "root\t\nA\troot\n",
        "m\t1\tG\tA\tC\n",
        "Y\t1\t.\tC\tG\t.\t.\t.\tGT\t8\t0\n",
    );
    fs::write(dir.path().join("calls"), "previous result").unwrap();
    let out = run(dir.path(), &[]);
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("allele index out of range"));
    assert_eq!(
        fs::read_to_string(dir.path().join("calls")).unwrap(),
        "previous result"
    );
}

#[test]
fn upstream_long_and_short_options_produce_identical_tables() {
    let dir = dataset(
        "root\t\nA\troot\n",
        "m\t1\tG\tA\tC\n",
        "Y\t1\t.\tC\tG\t.\t.\t.\tGT\t1\t0\n",
    );
    let p = dir.path();
    fs::write(p.join("ids"), "S1\n").unwrap();
    fs::write(p.join("ref"), ">Y\nC\n").unwrap();
    fs::write(p.join("weights"), "1\tC\tG\t2\n").unwrap();
    success(run(
        p,
        &[
            "-i", "ids", "-r", "1-1", "-c", "Y", "-f", "ref", "-w", "weights", "-x", "matrix", "-m",
        ],
    ));
    success(
        Command::new(env!("CARGO_BIN_EXE_haplogrouper-rs"))
            .current_dir(p)
            .args([
                "--vcfFile",
                "input",
                "--outFile",
                "long",
                "--hGrpLocusFile",
                "loci",
                "--hGrpTreeFile",
                "tree",
                "--IDListFile",
                "ids",
                "--regions",
                "1-1",
                "--chrom",
                "Y",
                "--referenceFasta",
                "ref",
                "--weightFile",
                "weights",
                "--verboseFile",
                "longmatrix",
                "--mismatchDetails",
            ])
            .output()
            .unwrap(),
    );
    for (a, b) in [("calls", "long"), ("matrix", "longmatrix")] {
        assert_eq!(fs::read(p.join(a)).unwrap(), fs::read(p.join(b)).unwrap());
    }
}
