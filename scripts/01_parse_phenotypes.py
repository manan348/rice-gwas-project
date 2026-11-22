"""Stage parse_pheno: extract blast phenotypes from Lin et al. 2018 Table S1 (PDF)
and match them to genotyped accessions on NSFTV ID.

Outputs (work/pheno/):
  phenotypes_parsed.csv          every accession row parsed from the PDF
  genotype_phenotype_matched.csv parsed rows that are present in the genotype .fam
  summary.json                   counts used by later stages and the README
"""
import re

import pandas as pd
import pdfplumber

from gwas_utils import ensure_dir, load_params, read_fam, write_json

TRAIT_COLS = ["D41-2_LT", "D41-2_DLA", "12YL_LT", "12YL_DLA"]

# One table row: GSOR  NSFTV  Country  Subgroup  Name  LT(D41-2) DLA(D41-2) LT(12YL) DLA(12YL)
# Country and Name are lazy free text; the subgroup code anchors the split between them and
# the four trailing numeric/'-' fields anchor the end of the name.
ROW_RE = re.compile(
    r"^(\d{6})\s+(\d+)\s+(.+?)\s+(TEJ|IND|AUS|TRJ|ADMIX|AROMATIC)\s+(.+?)\s+"
    r"([\d.]+|-)\s+([\d.]+|-)\s+([\d.]+|-)\s+([\d.]+|-)\s*$"
)

# Counts from the original (published) run; a mismatch means the PDF or parser changed.
EXPECTED = {"parsed": 314, "D41-2": 304, "12YL": 228}


def main() -> None:
    params = load_params()
    out = ensure_dir("work/pheno")

    with pdfplumber.open(params["inputs"]["phenotype_pdf"]) as pdf:
        n_pages = len(pdf.pages)
        lines = []
        for page in pdf.pages:
            lines.extend((page.extract_text() or "").split("\n"))

    rows = [m.groups() for m in (ROW_RE.match(l.strip()) for l in lines) if m]
    df = pd.DataFrame(rows, columns=["GSOR", "NSFTV_ID", "Country", "Subgroup", "Name", *TRAIT_COLS])
    for col in TRAIT_COLS:
        df[col] = pd.to_numeric(df[col].replace("-", pd.NA), errors="coerce")
    df["NSFTV_ID"] = df["NSFTV_ID"].astype(int)

    dup_ids = df["NSFTV_ID"][df["NSFTV_ID"].duplicated()].tolist()
    if dup_ids:
        raise SystemExit(f"ERROR: NSFTV IDs parsed more than once: {dup_ids}")

    counts = {
        "pdf_pages": n_pages,
        "parsed": int(len(df)),
        "D41-2": int(df["D41-2_LT"].notna().sum()),
        "12YL": int(df["12YL_LT"].notna().sum()),
    }
    print(counts)
    for key, val in EXPECTED.items():
        if counts[key] != val:
            raise SystemExit(f"ERROR: expected {val} '{key}' rows, parsed {counts[key]}")

    fam = read_fam(params["inputs"]["genotype_prefix"] + ".fam")
    merged = fam.merge(df, on="NSFTV_ID", how="inner")
    counts["genotyped"] = int(len(fam))
    counts["matched"] = int(len(merged))
    print(f"Matched {len(merged)} of {len(fam)} genotyped accessions to phenotype data")

    df.to_csv(out / "phenotypes_parsed.csv", index=False)
    merged.to_csv(out / "genotype_phenotype_matched.csv", index=False)
    write_json(counts, out / "summary.json")


if __name__ == "__main__":
    main()
