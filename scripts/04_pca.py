"""Stage pca: LD-prune the de-duplicated panel and compute genotype PCs.

Outputs: work/pca/ (prune lists, sativa_pca.eigenval/.eigenvec, summary.json)
         results/figures/pca_subpopulation.png
"""
import pandas as pd

from gwas_utils import ensure_dir, find_binary, load_params, run, set_plot_style, write_json

SUBGROUP_COLORS = {
    "IND": "#D55E00", "AUS": "#E69F00", "TEJ": "#0072B2", "TRJ": "#56B4E9",
    "AROMATIC": "#009E73", "ADMIX": "#CC79A7",
}


def main() -> None:
    params = load_params()
    p = params["pca"]
    plink = find_binary("plink")
    out = ensure_dir("work/pca")
    dedup = "work/dedup/sativa_qc_dedup"

    run([plink, "--bfile", dedup, "--indep-pairwise", str(p["ld_window"]), str(p["ld_step"]),
         str(p["ld_r2"]), "--out", str(out / "sativa_dedup_pruned")])
    run([plink, "--bfile", dedup, "--extract", str(out / "sativa_dedup_pruned.prune.in"),
         "--pca", str(p["n_pcs_compute"]), "--out", str(out / "sativa_pca")])

    n_pc = p["n_pcs_compute"]
    eigenval = pd.read_csv(out / "sativa_pca.eigenval", header=None, names=["eigenval"])
    # Share of the variance captured by the computed PCs (not of total genetic variance).
    eigenval["pct_of_computed"] = 100 * eigenval["eigenval"] / eigenval["eigenval"].sum()
    eigenvec = pd.read_csv(out / "sativa_pca.eigenvec", sep=r"\s+", header=None,
                           names=["FamilyID", "NSFTV_ID"] + [f"PC{i}" for i in range(1, n_pc + 1)])

    pheno = pd.read_csv("work/pheno/genotype_phenotype_matched.csv")
    lab = eigenvec.merge(pheno[["NSFTV_ID", "Subgroup"]], on="NSFTV_ID", how="left")
    lab["Subgroup"] = lab["Subgroup"].fillna("Unphenotyped (no label)")

    set_plot_style()
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(8, 6))
    unl = lab[lab["Subgroup"].str.startswith("Unphenotyped")]
    ax.scatter(unl["PC1"], unl["PC2"], c="#BBBBBB", s=22, alpha=0.7, edgecolors="none",
               label=f"Unlabelled (n={len(unl)})")
    for sg, color in SUBGROUP_COLORS.items():
        g = lab[lab["Subgroup"] == sg]
        if len(g):
            ax.scatter(g["PC1"], g["PC2"], c=color, s=26, alpha=0.85, edgecolors="none",
                       label=f"{sg} (n={len(g)})")
    pct = eigenval["pct_of_computed"]
    ax.set_xlabel(f"PC1 ({pct[0]:.1f}% of variance in top {n_pc} PCs)")
    ax.set_ylabel(f"PC2 ({pct[1]:.1f}% of variance in top {n_pc} PCs)")
    ax.set_title(f"PCA of {len(lab)} accessions (LD-pruned SNPs)", fontweight="bold")
    ax.legend(title="Subgroup", frameon=False, fontsize=9, loc="best")
    plt.tight_layout()
    ensure_dir("results/figures")
    fig.savefig("results/figures/pca_subpopulation.png", dpi=300, bbox_inches="tight")

    n_pruned = sum(1 for _ in open(out / "sativa_dedup_pruned.prune.in", encoding="utf-8"))
    write_json({
        "ld_pruned_snps": n_pruned,
        "accessions": int(len(lab)),
        "eigenval": [round(float(v), 4) for v in eigenval["eigenval"]],
        "pct_of_computed": [round(float(v), 2) for v in pct],
    }, out / "summary.json")


if __name__ == "__main__":
    main()
