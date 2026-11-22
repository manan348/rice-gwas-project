"""Stage qc: PLINK genotype QC (--geno / --maf / --mind) on the full 413-accession panel.

Outputs (work/qc/): sativa_qc.{bed,bim,fam,log,irem,nosex}, summary.json
"""
import pandas as pd

from gwas_utils import ensure_dir, find_binary, load_params, read_fam, run, write_json


def main() -> None:
    params = load_params()
    qc = params["qc"]
    prefix = params["inputs"]["genotype_prefix"]
    out = ensure_dir("work/qc")

    run([
        find_binary("plink"), "--file", prefix,
        "--geno", str(qc["geno"]), "--maf", str(qc["maf"]), "--mind", str(qc["mind"]),
        "--make-bed", "--out", str(out / "sativa_qc"),
    ])

    n_snps_in = sum(1 for _ in open(prefix + ".map", encoding="utf-8"))
    n_acc_in = len(read_fam(prefix + ".fam"))
    bim = pd.read_csv(out / "sativa_qc.bim", sep=r"\s+", header=None)
    fam = read_fam(out / "sativa_qc.fam")
    summary = {
        "snps_in": n_snps_in,
        "accessions_in": n_acc_in,
        "snps_pass": int(len(bim)),
        "accessions_pass": int(len(fam)),
    }
    print(summary)
    write_json(summary, out / "summary.json")


if __name__ == "__main__":
    main()
