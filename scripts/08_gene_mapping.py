"""Stage gene_mapping: map suggestive-or-better SNPs to MSU gene models and check the
validation genes (Pita, Ptr/Pita2).

Coordinates: GEMMA positions ('ps') are MSU6 (IRGSP build 5), the build of the 44K .map file.
Gene models are MSU7 (IRGSP-1.0). Every SNP is first lifted to MSU7 with RiceDiversity's
official per-SNP table (column 'snp_pos_msu7'); all distances are MSU7 vs MSU7. The shift
between builds is ~2.6 kb on chr12, ~16 kb on chr10 and ~466 kb on chr11 at the hit SNPs, so
mixing builds would give wrong genes.

Rules (documented so the tables are reproducible, not hand-picked):
  * snp_gene_pairs: every gene overlapping +/- window_bp of each hit SNP.
  * top_candidate_genes: for each hit SNP, its nearest gene (distance 0 = SNP inside the gene);
    genes are aggregated over SNPs and keep their strongest p-value. Known validation genes
    within the window are always listed. Genes MSU flags as TE-related (or whose note mentions
    a transposon / retrotransposon) are flagged low-confidence (repeat regions are prone to genotyping and
    mapping artefacts).
  * validation_pita_ptr_hits: every hit SNP within window_bp of a validation gene, with its
    LD (r2) to that trait's lead SNP on the same chromosome, computed in the accessions
    analysed for that trait.

Outputs: work/mapping/ (snp_gene_pairs.csv, LD files)
         results/tables/top_candidate_genes.csv, results/tables/validation_pita_ptr_hits.csv
         results/figures/candidate_genes_chr12.png
"""
import re

import numpy as np
import pandas as pd

from gwas_utils import (BONF_COLOR, PALETTE, SUGG_COLOR, ensure_dir, find_binary, load_params,
                        read_assoc, read_fam, run, set_plot_style, thresholds)

TE_RE = re.compile(r"transpos|retrotranspos", re.IGNORECASE)


def load_genes(path: str) -> pd.DataFrame:
    """MSU7 locus table (one row per gene model) -> one row per locus spanning all models."""
    m = pd.read_csv(path, sep="\t", dtype={"chr": str})
    m["chr_str"] = m["chr"].str.replace("Chr", "", regex=False)
    m = m[m["chr_str"].str.isnumeric()].copy()
    m["chr"] = m["chr_str"].astype(int)
    rep = m.sort_values("is_representative", ascending=False)  # 'Y' first
    genes = rep.groupby("locus").agg(
        chr=("chr", "first"), start=("start", "min"), end=("stop", "max"),
        annotation=("annotation", "first"), is_TE=("is_TE", "first"),
    ).reset_index().rename(columns={"locus": "gene_id"})
    genes["is_TE"] = genes["is_TE"].eq("Y") | genes["annotation"].str.contains(TE_RE, na=False)
    return genes[["chr", "start", "end", "gene_id", "annotation", "is_TE"]]


def load_liftover(path: str) -> pd.Series:
    """RiceDiversity 44K table: SNP id -> MSU7 (IRGSP-1.0) position; NaN if unmappable."""
    t = pd.read_csv(path, sep="\t", header=None, skiprows=1,
                    names=["chr", "rs", "cM", "pos_v4", "pos_msu6", "pos_msu7", "x1", "x2"])
    return pd.to_numeric(t.set_index("rs")["pos_msu7"], errors="coerce")


def distance(start, end, pos):
    return np.where((start <= pos) & (pos <= end), 0,
                    np.minimum(np.abs(start - pos), np.abs(end - pos)))


def ld_to_lead(plink, bfile, keep_ids, fam, lead_rs, chrom, out_prefix) -> pd.DataFrame:
    keep = fam[fam["NSFTV_ID"].isin(keep_ids)][["FamilyID", "NSFTV_ID"]]
    keep.to_csv(f"{out_prefix}.keep", sep="\t", header=False, index=False)
    run([plink, "--bfile", bfile, "--keep", f"{out_prefix}.keep", "--chr", str(chrom),
         "--r2", "--ld-snp", lead_rs, "--ld-window", "999999", "--ld-window-kb", "100000",
         "--ld-window-r2", "0", "--out", out_prefix])
    ld = pd.read_csv(f"{out_prefix}.ld", sep=r"\s+")
    return ld[["SNP_B", "R2"]].rename(columns={"SNP_B": "rs", "R2": "r2_with_lead"})


def main() -> None:
    params = load_params()
    mp = params["mapping"]
    window = int(mp["window_bp"])
    known = {v: k for k, v in mp["known_genes"].items()}  # gene_id -> name
    out = ensure_dir("work/mapping")
    plink = find_binary("plink")
    bfile = "work/dedup/sativa_qc_dedup"

    genes = load_genes(params["inputs"]["gene_annotation"])
    lift = load_liftover(params["inputs"]["snp_liftover"])
    missing = set(known) - set(genes["gene_id"])
    if missing:
        raise SystemExit(f"ERROR: validation genes not in annotation: {missing}")
    hits = pd.read_csv("results/tables/suggestive_hits_all_traits.csv")
    hits["snp_pos_msu7"] = hits["rs"].map(lift)
    unmapped = hits[hits["snp_pos_msu7"].isna()]
    if len(unmapped):
        print(f"WARNING: {len(unmapped)} hit SNPs have no MSU7 position and are not mapped:",
              unmapped["rs"].tolist())
    hits = hits.dropna(subset=["snp_pos_msu7"]).astype({"snp_pos_msu7": int})

    # ---- all SNP-gene pairs within the window
    pairs = []
    for _, h in hits.iterrows():
        pos = h["snp_pos_msu7"]
        g = genes[(genes["chr"] == h["chr"]) & (genes["end"] >= pos - window)
                  & (genes["start"] <= pos + window)].copy()
        g["distance_bp"] = distance(g["start"].to_numpy(), g["end"].to_numpy(), pos)
        for c in ["trait", "rs", "ps", "snp_pos_msu7", "p_wald", "beta"]:
            g[c] = h[c]
        pairs.append(g)
    pairs = pd.concat(pairs).rename(columns={"ps": "snp_pos_msu6", "start": "gene_start",
                                             "end": "gene_end"})
    pairs = pairs[["trait", "chr", "rs", "snp_pos_msu6", "snp_pos_msu7", "p_wald", "beta",
                   "gene_id", "gene_start", "gene_end", "distance_bp", "is_TE", "annotation"]]
    pairs = pairs.sort_values(["chr", "snp_pos_msu7", "trait", "distance_bp"]).reset_index(drop=True)
    pairs.to_csv(out / "snp_gene_pairs.csv", index=False)
    print(f"{len(pairs)} SNP-gene pairs within {window:,} bp of {len(hits)} hit SNPs")

    # ---- candidate genes: nearest gene per SNP (+ validation genes in window)
    nearest = pairs.loc[pairs.groupby(["trait", "rs"])["distance_bp"].idxmin()]
    cand = pd.concat([nearest, pairs[pairs["gene_id"].isin(known)]])
    cand = cand.sort_values("p_wald")
    agg = cand.groupby("gene_id").agg(
        chr=("chr", "first"), gene_start=("gene_start", "first"), gene_end=("gene_end", "first"),
        best_p=("p_wald", "min"), best_snp=("rs", "first"),
        min_distance_bp=("distance_bp", "min"), n_hit_snps=("rs", "nunique"),
        traits=("trait", lambda s: ";".join(sorted(set(s)))),
        annotation=("annotation", "first"), is_TE=("is_TE", "first"),
    ).reset_index()
    agg["known_gene"] = agg["gene_id"].map(known).fillna("")
    agg["confidence"] = np.select(
        [agg["known_gene"] != "", agg["is_TE"]],
        ["known R gene", "low (transposon-annotated)"], default="candidate")
    agg = agg.sort_values("best_p").reset_index(drop=True)
    agg["best_p"] = agg["best_p"].map(lambda p: f"{p:.4g}")
    agg = agg.drop(columns="is_TE")
    agg.to_csv("results/tables/top_candidate_genes.csv", index=False)
    print(agg.to_string(index=False))

    # ---- validation genes with LD to each trait's lead SNP on that chromosome
    fam = read_fam(bfile + ".fam")
    pheno = pd.read_csv("work/pheno/genotype_phenotype_matched.csv")
    val = pairs[pairs["gene_id"].isin(known)].copy()
    val["known_gene"] = val["gene_id"].map(known)
    ld_parts = []
    for (trait, chrom), grp in val.groupby(["trait", "chr"]):
        th = hits[(hits["trait"] == trait) & (hits["chr"] == chrom)]
        lead = th.loc[th["p_wald"].idxmin(), "rs"]
        ids = pheno.loc[pheno[trait].notna(), "NSFTV_ID"]
        ld = ld_to_lead(plink, bfile, ids, fam, lead, chrom, str(out / f"ld_{trait}_chr{chrom}"))
        ld["trait"], ld["lead_snp"] = trait, lead
        ld_parts.append(ld)
    if ld_parts:
        val = val.merge(pd.concat(ld_parts), on=["trait", "rs"], how="left")
        val["r2_with_lead"] = val["r2_with_lead"].round(3)
    val = val.sort_values(["known_gene", "distance_bp", "trait"])
    val[["known_gene", "gene_id", "gene_start", "gene_end", "trait", "chr", "rs", "snp_pos_msu6",
         "snp_pos_msu7", "p_wald", "distance_bp", "lead_snp", "r2_with_lead", "annotation"]].to_csv(
        "results/tables/validation_pita_ptr_hits.csv", index=False)

    # ---- chr12 regional figure
    set_plot_style()
    import matplotlib.pyplot as plt

    lo, hi = mp["chr12_region"]
    alpha = params["thresholds"]["bonferroni_alpha"]
    kg = genes[genes["gene_id"].isin(known)]
    traits = [t for t in params["traits"]
              if ((hits["trait"] == t["id"]) & (hits["chr"] == 12)).any()]
    fig, axes = plt.subplots(len(traits), 1, figsize=(11, 3.2 * len(traits)), sharex=True,
                             squeeze=False)
    for ax, t in zip(axes[:, 0], traits):
        a = read_assoc(f"work/gemma/{t['slug']}.assoc.txt")
        thr = thresholds(len(a), alpha)
        a["pos_msu7"] = a["rs"].map(lift)
        r = a[(a["chr"] == 12) & a["pos_msu7"].between(lo, hi)]
        ax.scatter(r["pos_msu7"] / 1e6, -np.log10(r["p_wald"]), s=14, color=PALETTE[0],
                   edgecolors="none", alpha=0.8)
        ax.axhline(-np.log10(thr["bonferroni"]), color=BONF_COLOR, ls="--", lw=1)
        ax.axhline(-np.log10(thr["suggestive"]), color=SUGG_COLOR, ls=":", lw=1)
        for _, g in kg.iterrows():
            ax.axvspan(g["start"] / 1e6, g["end"] / 1e6, color=BONF_COLOR, alpha=0.25, lw=0)
            ax.axvline((g["start"] + g["end"]) / 2e6, color=BONF_COLOR, lw=0.6, alpha=0.6)
        ax.set_ylabel(r"$-\log_{10}(p)$")
        ax.set_title(t["label"], loc="left", fontsize=11, fontweight="bold")
        ax.grid(axis="y", alpha=0.15)
    top = axes[0, 0]
    top.set_title(traits[0]["label"], loc="left", fontsize=11, fontweight="bold", pad=34)
    ymax = top.get_ylim()[1]
    # Pita and Ptr are only ~200 kb apart: put the left gene's label to its left, the rest right.
    for i, (_, g) in enumerate(kg.sort_values("start").iterrows()):
        name = known[g["gene_id"]]
        ha, dx = ("right", -3) if i == 0 else ("left", 3)
        top.annotate(f"{name}\n{g['gene_id']}", xy=((g["start"] + g["end"]) / 2e6, ymax),
                     xytext=(dx, 4), textcoords="offset points", ha=ha, va="bottom",
                     fontsize=8, color=BONF_COLOR, fontweight="bold")
    axes[-1, 0].set_xlabel("Chromosome 12 position (Mb, MSU7 / IRGSP-1.0)")
    fig.suptitle("Chr12 association peak and the cloned blast-resistance genes", y=1.02,
                 fontweight="bold")
    fig.tight_layout()
    fig.savefig("results/figures/candidate_genes_chr12.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
