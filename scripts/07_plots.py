"""Stage plots: Manhattan + QQ per trait, a 4-trait Manhattan panel, and a hits-per-chromosome bar.

All figures are written to results/figures/ at 300 dpi.
"""
import numpy as np
import pandas as pd

from gwas_utils import (BONF_COLOR, PALETTE, SUGG_COLOR, cumulative_positions,
                        expected_quantiles, lambda_gc, load_params, read_assoc,
                        set_plot_style, thresholds)

set_plot_style()
import matplotlib.pyplot as plt  # noqa: E402  (after backend selection)

FIG = "results/figures"


def draw_manhattan(ax, assoc, thr, title, compact=False):
    d, centers = cumulative_positions(assoc)
    for i, chrom in enumerate(sorted(d["chr"].unique())):
        sub = d[d["chr"] == chrom]
        ax.scatter(sub["pos_cum"], sub["neg_log_p"], c=PALETTE[i % 2], s=6 if compact else 10,
                   alpha=0.85, edgecolors="none", rasterized=True)
    ax.axhline(-np.log10(thr["bonferroni"]), color=BONF_COLOR, ls="--", lw=1.2,
               label=f"Bonferroni (p = {thr['bonferroni']:.1e})")
    ax.axhline(-np.log10(thr["suggestive"]), color=SUGG_COLOR, ls=":", lw=1.2,
               label=f"Suggestive (p = {thr['suggestive']:.1e})")
    ax.set_xticks([c for _, c in centers])
    ax.set_xticklabels([str(k) for k, _ in centers], fontsize=8 if compact else 10)
    ax.set_xlim(d["pos_cum"].min(), d["pos_cum"].max())
    top = max(d["neg_log_p"].max(), -np.log10(thr["bonferroni"]))
    ax.set_ylim(0, top + 0.8)
    ax.set_ylabel(r"$-\log_{10}(p)$")
    ax.grid(axis="y", alpha=0.15)
    if compact:
        ax.set_title(title, fontsize=11, fontweight="bold", loc="left")
    else:
        ax.set_title(f"GWAS of {title}", fontsize=14, fontweight="bold", pad=12)
        ax.set_xlabel("Chromosome")
    ax.legend(frameon=False, loc="upper left", fontsize=8 if compact else 10)


def draw_qq(assoc, title, path):
    p = np.sort(assoc["p_wald"].to_numpy())
    observed = -np.log10(p)
    expected = expected_quantiles(len(p))
    lam = lambda_gc(p)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(expected, observed, s=10, alpha=0.5, color=PALETTE[0], edgecolors="none")
    m = max(expected.max(), observed.max()) + 0.3
    ax.plot([0, m], [0, m], color=BONF_COLOR, ls="--", lw=1.2)
    ax.set_xlim(0, m)
    ax.set_ylim(0, m)
    ax.set_xlabel(r"Expected $-\log_{10}(p)$")
    ax.set_ylabel(r"Observed $-\log_{10}(p)$")
    ax.set_title(f"QQ Plot: {title}", fontsize=13, fontweight="bold", pad=12)
    ax.text(0.05, 0.92, r"$\lambda_{GC}$" + f" = {lam:.3f}", transform=ax.transAxes,
            bbox=dict(facecolor="white", edgecolor="#cccccc", boxstyle="round,pad=0.4"))
    ax.grid(alpha=0.15)
    fig.tight_layout()
    fig.savefig(path, dpi=300, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    params = load_params()
    alpha = params["thresholds"]["bonferroni_alpha"]
    traits = params["traits"]
    data = {t["id"]: read_assoc(f"work/gemma/{t['slug']}.assoc.txt") for t in traits}

    for t in traits:
        a = data[t["id"]]
        thr = thresholds(len(a), alpha)
        fig, ax = plt.subplots(figsize=(15, 5.5))
        draw_manhattan(ax, a, thr, t["label"])
        fig.tight_layout()
        fig.savefig(f"{FIG}/manhattan_{t['slug']}.png", dpi=300, bbox_inches="tight")
        plt.close(fig)
        draw_qq(a, t["label"], f"{FIG}/qq_{t['slug']}.png")

    fig, axes = plt.subplots(len(traits), 1, figsize=(15, 4 * len(traits)))
    for ax, t in zip(axes, traits):
        a = data[t["id"]]
        draw_manhattan(ax, a, thresholds(len(a), alpha), t["label"], compact=True)
    axes[-1].set_xlabel("Chromosome")
    fig.tight_layout()
    fig.savefig(f"{FIG}/manhattan_all_traits_panel.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    hits = pd.read_csv("results/tables/suggestive_hits_all_traits.csv")
    counts = (hits.groupby(["chr", "trait"]).size().unstack(fill_value=0)
              .reindex(index=range(1, 13), columns=[t["id"] for t in traits], fill_value=0))
    counts = counts.loc[counts.sum(axis=1) > 0]
    fig, ax = plt.subplots(figsize=(10, 5))
    counts.plot(kind="bar", stacked=True, ax=ax,
                color=["#3B5B92", "#8FA6C9", "#C0392B", "#E8A87C"][: counts.shape[1]])
    ax.set_xlabel("Chromosome")
    ax.set_ylabel("Suggestive-or-better SNPs")
    ax.set_title("Suggestive-or-better SNPs by chromosome and trait", fontweight="bold")
    labels = [f"{t['id']} ({int(counts[t['id']].sum())})" for t in traits]
    ax.legend(labels, title="Trait (n SNPs)", frameon=False)
    plt.xticks(rotation=0)
    fig.tight_layout()
    fig.savefig(f"{FIG}/hits_summary_by_chromosome.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
