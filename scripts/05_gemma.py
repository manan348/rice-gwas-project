"""Stage gemma: centred kinship + univariate LMM (-lmm 4) for each trait.

The phenotype file's first column is a dummy trait observed for every accession, so the
kinship matrix is estimated from all de-duplicated accessions (GEMMA only uses individuals
with a non-missing value in the selected phenotype column).

Outputs (work/gemma/): gemma_phenotypes.txt, gemma_covariates.txt, kinship.cXX.txt,
                      <slug>.assoc.txt / <slug>.log.txt for each trait
"""
import pandas as pd

from gwas_utils import ensure_dir, find_binary, load_params, read_fam, run


def main() -> None:
    params = load_params()
    g = params["gemma"]
    traits = params["traits"]
    gemma = find_binary("gemma")
    out = ensure_dir("work/gemma")
    bfile = "work/dedup/sativa_qc_dedup"

    fam = read_fam(bfile + ".fam")
    pheno = pd.read_csv("work/pheno/genotype_phenotype_matched.csv")
    trait_ids = [t["id"] for t in traits]
    ph = fam[["NSFTV_ID"]].merge(pheno[["NSFTV_ID", *trait_ids]], on="NSFTV_ID", how="left")
    assert len(ph) == len(fam) and (ph["NSFTV_ID"].values == fam["NSFTV_ID"].values).all()
    ph.insert(1, "dummy_all_valid", 1)
    ph[["dummy_all_valid", *trait_ids]].to_csv(
        out / "gemma_phenotypes.txt", sep="\t", header=False, index=False, na_rep="NA")

    n_pcs = g["n_pcs"]
    vec = pd.read_csv("work/pca/sativa_pca.eigenvec", sep=r"\s+", header=None)
    vec = vec.iloc[:, : 2 + n_pcs]
    vec.columns = ["FamilyID", "NSFTV_ID"] + [f"PC{i}" for i in range(1, n_pcs + 1)]
    cov = fam[["NSFTV_ID"]].merge(vec.drop(columns="FamilyID"), on="NSFTV_ID", how="left")
    assert cov.notna().all().all(), "every accession must have PCs"
    cov.insert(1, "intercept", 1)
    cov.drop(columns="NSFTV_ID").to_csv(out / "gemma_covariates.txt", sep="\t", header=False,
                                        index=False)

    common = ["-bfile", bfile, "-p", str(out / "gemma_phenotypes.txt"), "-outdir", str(out)]
    run([gemma, *common, "-n", "1", "-gk", "1", "-o", "kinship"])

    for i, t in enumerate(traits, start=2):  # column 1 is the dummy trait
        run([gemma, *common, "-k", str(out / "kinship.cXX.txt"),
             "-c", str(out / "gemma_covariates.txt"), "-n", str(i), "-lmm", "4",
             "-maf", str(g["maf"]), "-miss", str(g["miss"]), "-o", t["slug"]])


if __name__ == "__main__":
    main()
