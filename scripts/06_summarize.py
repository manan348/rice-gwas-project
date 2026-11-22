"""Stage summarize: per-trait run statistics and the suggestive-or-better hit table.

Thresholds are computed per trait from the number of SNPs GEMMA actually tested:
  Bonferroni = alpha / n_tested,  suggestive = 1 / n_tested.

Outputs: results/tables/gwas_run_summary.csv
         results/tables/suggestive_hits_all_traits.csv
         results/metrics.json  (pipeline-wide counts, tracked as DVC metrics)
"""
import re

import pandas as pd

from gwas_utils import (classify, ensure_dir, lambda_gc, load_params, read_assoc, read_json,
                        thresholds, write_json)

LOG_PATTERNS = {
    "n": r"number of analyzed individuals\s*=\s*(\d+)",
    "n_snps_tested": r"number of analyzed SNPs(?:/var)?\s*=\s*(\d+)",
    "PVE": r"pve estimate in the null model\s*=\s*([-\d.eE]+)",
    "PVE_SE": r"se\(pve\) in the null model\s*=\s*([-\d.eE]+)",
}


def parse_gemma_log(path: str) -> dict:
    text = open(path, encoding="utf-8").read()
    vals = {}
    for key, pat in LOG_PATTERNS.items():
        m = re.search(pat, text)
        if not m:
            raise SystemExit(f"ERROR: '{key}' not found in {path}")
        vals[key] = float(m.group(1))
    return vals


def main() -> None:
    params = load_params()
    alpha = params["thresholds"]["bonferroni_alpha"]
    ensure_dir("results/tables")

    summary_rows, hits = [], []
    for t in params["traits"]:
        assoc = read_assoc(f"work/gemma/{t['slug']}.assoc.txt")
        log = parse_gemma_log(f"work/gemma/{t['slug']}.log.txt")
        n_tested = len(assoc)
        if n_tested != int(log["n_snps_tested"]):
            raise SystemExit(f"ERROR: {t['id']}: assoc rows {n_tested} != log {log['n_snps_tested']}")
        thr = thresholds(n_tested, alpha)
        top_p = float(assoc["p_wald"].min())
        isolate, trait = t["label"].split(" ", 1)   # "12YL-DL-3-2 Lesion Type"
        summary_rows.append({
            "trait_id": t["id"], "isolate": isolate, "trait": trait,
            "n": int(log["n"]), "n_snps_tested": n_tested,
            "PVE": round(log["PVE"], 3), "PVE_SE": round(log["PVE_SE"], 3),
            "lambda_GC": round(lambda_gc(assoc["p_wald"]), 3),
            "top_p": f"{top_p:.3g}",
            "bonferroni_p": f"{thr['bonferroni']:.3g}", "suggestive_p": f"{thr['suggestive']:.3g}",
            "n_bonferroni_snps": int((assoc["p_wald"] <= thr["bonferroni"]).sum()),
            "n_suggestive_or_better_snps": int((assoc["p_wald"] <= thr["suggestive"]).sum()),
            "significance": classify(top_p, thr),
        })
        h = assoc[assoc["p_wald"] <= thr["suggestive"]].copy()
        h["trait"] = t["id"]
        h["bonferroni"] = h["p_wald"] <= thr["bonferroni"]
        hits.append(h)

    summary = pd.DataFrame(summary_rows)
    summary.to_csv("results/tables/gwas_run_summary.csv", index=False)
    print(summary.to_string(index=False))

    all_hits = pd.concat(hits)[["trait", "chr", "rs", "ps", "allele1", "allele0", "af",
                                "beta", "se", "p_wald", "bonferroni"]]
    all_hits = all_hits.sort_values(["chr", "ps", "trait"]).reset_index(drop=True)
    all_hits.to_csv("results/tables/suggestive_hits_all_traits.csv", index=False)
    print(f"{len(all_hits)} suggestive-or-better SNPs across {len(summary)} traits")

    metrics = {
        "phenotypes": read_json("work/pheno/summary.json"),
        "qc": read_json("work/qc/summary.json"),
        "dedup": {k: v for k, v in read_json("work/dedup/summary.json").items()
                  if not k.endswith("_ids") and k != "phenotyped_lost_to_mind"},
        "pca": {k: v for k, v in read_json("work/pca/summary.json").items() if k != "eigenval"},
        "gwas": {r["trait_id"]: {k: v for k, v in r.items() if k != "trait_id"}
                 for r in summary_rows},
        "suggestive_or_better_snps_total": int(len(all_hits)),
    }
    write_json(metrics, "results/metrics.json")


if __name__ == "__main__":
    main()
