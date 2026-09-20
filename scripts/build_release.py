"""Create a deterministic, checksummed archive of the public repository files."""

from pathlib import Path
import argparse
import hashlib
import re
import zipfile

ROOT = Path(__file__).resolve().parents[1]
PUBLIC_ROOT_ENTRIES = {
    "src",
    "scripts",
    "figures",
    "tables",
    "tests",
    "docs",
    "metadata",
    ".github",
    "README.md",
    "LICENSE",
    "DATA_LICENSE.md",
    "CITATION.cff",
    "pyproject.toml",
    "requirements-tested.txt",
    ".gitignore",
    ".gitattributes",
    ".zenodo.json",
}
SKIP = {
    ".git",
    ".venv",
    "__pycache__",
    ".pytest_cache",
    ".ruff_cache",
    "results",
    "work",
    "raw",
    "dist",
    "build",
}


def public_files():
    """Package the public tree only, excluding outputs, caches and local secrets."""
    return sorted(
        p
        for p in ROOT.rglob("*")
        if p.is_file()
        and p.relative_to(ROOT).parts[0] in PUBLIC_ROOT_ENTRIES
        and not any(
            x in SKIP or x.endswith(".egg-info") for x in p.relative_to(ROOT).parts
        )
        and p.name != "SHA256SUMS.txt"
        and not p.name.startswith(".env")
        and p.suffix.lower() not in {".pem", ".key", ".nc", ".parquet", ".bin", ".zip"}
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=ROOT / "dist")
    args = parser.parse_args()
    files = public_files()
    sums = "".join(
        f"{hashlib.sha256(p.read_bytes()).hexdigest()}  {p.relative_to(ROOT).as_posix()}\n"
        for p in files
    )
    manifest = ROOT / "metadata/SHA256SUMS.txt"
    manifest.write_text(sums, encoding="utf-8")
    files.append(manifest)
    # The installation metadata is the single source of the software version.
    version = re.search(
        r'^version = "([^"]+)"',
        (ROOT / "pyproject.toml").read_text(encoding="utf-8"),
        re.MULTILINE,
    ).group(1)
    name = f"TCEP-latitude-analysis-v{version}"
    args.output_dir.mkdir(parents=True, exist_ok=True)
    target = args.output_dir / (name + ".zip")
    with zipfile.ZipFile(target, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for p in sorted(files):
            entry = zipfile.ZipInfo(
                name + "/" + p.relative_to(ROOT).as_posix(),
                date_time=(2026, 9, 19, 0, 0, 0),
            )
            entry.compress_type = zipfile.ZIP_DEFLATED
            archive.writestr(entry, p.read_bytes())
    digest = hashlib.sha256(target.read_bytes()).hexdigest()
    target.with_suffix(".zip.sha256").write_text(
        f"{digest}  {target.name}\n", encoding="utf-8"
    )
    print(
        f"{target.name}: {len(files)} files, {target.stat().st_size:,} bytes, SHA256 {digest}"
    )


if __name__ == "__main__":
    main()
