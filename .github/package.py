"""Package a release build from the repository root (also used by Actions)."""
import json
import subprocess
import sys
import tarfile
from pathlib import Path


def package(target):
    windows = "windows" in target
    binary = "haplogrouper-rs.exe" if windows else "haplogrouper-rs"
    executable = Path("target", target, "release", binary)
    subprocess.run([executable, "--version"], check=True)
    files = {binary: executable}
    for name in ["LICENSE", "NOTICE", "LICENSES/Apache-2.0.txt", "README.md"]:
        files[name] = Path(name)

    # Carry the license notices shipped by the locked crates for this target.
    metadata = json.loads(subprocess.check_output([
        "cargo", "+stable", "metadata", "--locked", "--format-version", "1",
        "--filter-platform", target,
    ]))
    notices = []
    for crate in metadata["packages"]:
        if crate["source"] is None:
            continue
        name = f"{crate['name']}-{crate['version']}"
        notices.append(f"{name}  {crate['license']}  {crate['repository'] or ''}")
        directory = Path(crate["manifest_path"]).parent
        for path in sorted(directory.iterdir()):
            if path.is_file() and path.name.startswith(("LICENSE", "COPYING", "NOTICE")):
                files[f"LICENSES/crates/{name}/{path.name}"] = path
    destination = Path("target/archives")
    destination.mkdir(exist_ok=True)
    index = executable.parent / "crate-notices.txt"
    index.write_text("\n".join(notices) + "\n", encoding="utf-8")
    files["LICENSES/crates/index.txt"] = index
    with tarfile.open(destination / f"haplogrouper-rs-{target}.tar.gz", "w:gz") as archive:
        for name, path in files.items():
            archive.add(path, arcname=name)


if __name__ == "__main__":
    package(sys.argv[1])
