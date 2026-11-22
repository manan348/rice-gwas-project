"""Stage dedup: detect and remove near-duplicate accessions.

PLINK's PI_HAT (method-of-moments IBD) assumes outbred, HWE samples and saturates at ~1.0
for many unrelated pairs of inbred rice lines, so near-duplicates are called on raw
identity-by-state (DST) instead. Pairs with DST >= threshold where BOTH accessions are
phenotyped are joined into clusters (connected components); in each cluster the accession
with the lowest genotype missingness (ties -> lowest NSFTV ID) is kept.

Note: pairs involving an unphenotyped accession are not de-duplicated. Those accessions do
not enter the association tests (GEMMA drops missing phenotypes) but do enter the PCA.

Outputs: work/dedup/ (dedup fileset, genome/missing tables, summary.json)
         results/tables/duplicate_clusters.csv
"""
import networkx as nx
import pandas as pd

from gwas_utils import ensure_dir, find_binary, load_params, read_fam, run, write_json


def main() -> None:
    params = load_params()
    thr = params["dedup"]["dst_threshold"]
    plink = find_binary("plink")
    out = ensure_dir("work/dedup")
    qc = "work/qc/sativa_qc"

    run([plink, "--bfile", qc, "--genome", "--out", str(out / "sativa_ibd_full")])
    run([plink, "--bfile", qc, "--missing", "--out", str(out / "sativa_missing")])

    genome = pd.read_csv(out / "sativa_ibd_full.genome", sep=r"\s+")
    imiss = pd.read_csv(out / "sativa_missing.imiss", sep=r"\s+")
    qc_fam = read_fam(qc + ".fam")
    pheno = pd.read_csv("work/pheno/genotype_phenotype_matched.csv")
    pheno_ids = set(pheno["NSFTV_ID"])

    dups = genome[genome["DST"] >= thr].copy()
    dups["both_phenotyped"] = dups["IID1"].isin(pheno_ids) & dups["IID2"].isin(pheno_ids)
    pheno_dups = dups[dups["both_phenotyped"]]

    G = nx.Graph()
    G.add_edges_from(zip(pheno_dups["IID1"], pheno_dups["IID2"]))
    clusters = sorted((sorted(c) for c in nx.connected_components(G)), key=lambda c: c[0])

    rows, to_remove = [], []
    fmiss = imiss.set_index("IID")["F_MISS"]
    info = pheno.set_index("NSFTV_ID")
    for k, members in enumerate(clusters, start=1):
        keep = min(members, key=lambda iid: (fmiss[iid], iid))
        for iid in members:
            rows.append({
                "cluster": k, "NSFTV_ID": iid, "Name": info.at[iid, "Name"],
                "Country": info.at[iid, "Country"], "Subgroup": info.at[iid, "Subgroup"],
                "F_MISS": round(float(fmiss[iid]), 6), "kept": iid == keep,
            })
            if iid != keep:
                to_remove.append(iid)

    ensure_dir("results/tables")
    pd.DataFrame(rows).to_csv("results/tables/duplicate_clusters.csv", index=False)

    remove_df = qc_fam[qc_fam["NSFTV_ID"].isin(to_remove)][["FamilyID", "NSFTV_ID"]]
    remove_df.to_csv(out / "duplicates_to_remove.txt", sep="\t", header=False, index=False)
    run([plink, "--bfile", qc, "--remove", str(out / "duplicates_to_remove.txt"),
         "--make-bed", "--out", str(out / "sativa_qc_dedup")])

    dedup_fam = read_fam(out / "sativa_qc_dedup.fam")
    qc_ids = set(qc_fam["NSFTV_ID"])
    summary = {
        "pairs_total": int(len(genome)),
        "pairs_pi_hat_ge_0.99": int((genome["PI_HAT"] >= 0.99).sum()),
        "pairs_pi_hat_ge_0.5": int((genome["PI_HAT"] >= 0.5).sum()),
        "dst_threshold": thr,
        "pairs_dst_ge_threshold": int(len(dups)),
        "pairs_dst_both_phenotyped": int(len(pheno_dups)),
        "clusters": len(clusters),
        "accessions_in_clusters": int(sum(len(c) for c in clusters)),
        "removed": int(len(to_remove)),
        "removed_ids": sorted(int(i) for i in to_remove),
        "phenotyped_lost_to_mind": sorted(int(i) for i in pheno_ids - qc_ids),
        "phenotyped_after_qc": len(pheno_ids & qc_ids),
        "accessions_after_dedup": int(len(dedup_fam)),
        "phenotyped_after_dedup": int(dedup_fam["NSFTV_ID"].isin(pheno_ids).sum()),
    }
    print(summary)
    write_json(summary, out / "summary.json")


if __name__ == "__main__":
    main()
