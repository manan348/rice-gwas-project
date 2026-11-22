"""Shared helpers for the rice blast GWAS pipeline stages.

Every stage script is run from the repository root by `dvc repro` (see dvc.yaml).
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
from scipy import stats

REPO_ROOT = Path(__file__).resolve().parents[1]
FAM_COLS = ["FamilyID", "NSFTV_ID", "FatherID", "MotherID", "Sex", "Phenotype"]


# --------------------------------------------------------------------------- config / io
def load_params() -> dict:
    with open(REPO_ROOT / "params.yaml", encoding="utf-8") as fh:
        return yaml.safe_load(fh)


def ensure_dir(path: str | Path) -> Path:
    p = Path(path)
    p.mkdir(parents=True, exist_ok=True)
    return p


def write_json(obj: dict, path: str | Path) -> None:
    ensure_dir(Path(path).parent)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=2, sort_keys=False)
        fh.write("\n")


def read_json(path: str | Path) -> dict:
    with open(path, encoding="utf-8") as fh:
        return json.load(fh)


def read_fam(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep=r"\s+", header=None, names=FAM_COLS)


# --------------------------------------------------------------------------- binaries
def find_binary(name: str) -> str:
    """Locate PLINK/GEMMA: $PLINK / $GEMMA env var, then ./bin, then PATH."""
    env = os.environ.get(name.upper())
    if env:
        return env
    local = REPO_ROOT / "bin" / name
    if local.exists():
        return str(local)
    found = shutil.which(name)
    if found:
        return found
    sys.exit(f"ERROR: '{name}' not found. Run scripts/tools/get_binaries.sh or put it on PATH.")


def run(cmd: list[str]) -> None:
    """Run an external command, echoing it, and fail loudly on a non-zero exit."""
    print("+", " ".join(f'"{c}"' if " " in c else c for c in cmd), flush=True)
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    print(result.stdout, flush=True)
    if result.returncode != 0:
        sys.exit(f"ERROR: command failed with exit code {result.returncode}")


# --------------------------------------------------------------------------- statistics
def lambda_gc(pvalues: pd.Series | np.ndarray) -> float:
    """Genomic inflation factor: median observed chi2(1) / median expected chi2(1)."""
    chisq = stats.chi2.isf(np.asarray(pvalues, dtype=float), df=1)
    return float(np.median(chisq) / stats.chi2.ppf(0.5, df=1))


def expected_quantiles(n: int) -> np.ndarray:
    """Expected -log10 p under the null using (i - 0.5) / n plotting positions."""
    return -np.log10((np.arange(1, n + 1) - 0.5) / n)


def thresholds(n_tests: int, alpha: float) -> dict:
    return {"bonferroni": alpha / n_tests, "suggestive": 1.0 / n_tests}


def classify(top_p: float, thr: dict) -> str:
    if top_p <= thr["bonferroni"]:
        return "Bonferroni"
    if top_p <= thr["suggestive"]:
        return "suggestive only"
    return "no association"


# --------------------------------------------------------------------------- plotting
PALETTE = ["#3B5B92", "#8FA6C9"]
BONF_COLOR = "#C0392B"
SUGG_COLOR = "#2980B9"


def set_plot_style() -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#333333",
        "figure.facecolor": "white",
        "savefig.dpi": 300,
    })


def cumulative_positions(assoc: pd.DataFrame) -> tuple[pd.DataFrame, list[tuple[int, float]]]:
    """Add a genome-wide x coordinate ('pos_cum') for Manhattan plots."""
    d = assoc.sort_values(["chr", "ps"]).reset_index(drop=True).copy()
    d["pos_cum"] = 0
    offset = 0
    centers = []
    for chrom in sorted(d["chr"].unique()):
        mask = d["chr"] == chrom
        d.loc[mask, "pos_cum"] = d.loc[mask, "ps"] + offset
        centers.append((int(chrom), float(d.loc[mask, "pos_cum"].mean())))
        offset = d.loc[mask, "pos_cum"].max()
    d["neg_log_p"] = -np.log10(d["p_wald"])
    return d, centers


def read_assoc(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t")
