# Rice Blast Resistance GWAS

A genome-wide association study of blast disease (*Magnaporthe oryzae*) resistance in
the RDP1 rice diversity panel — a classical-inference counterpart to **NeuroCrop**, a
machine-learning maize yield predictor, in the same portfolio.

The whole analysis is a reproducible [DVC](https://dvc.org) pipeline: every table and
figure in `results/` is produced by the scripts in `scripts/` from the raw inputs in
`data/`, with all settings in `params.yaml`.

---

**Author:** Abdul Manan
*Plant Breeder | Machine & Deep Learning Researcher*
[abdulmanan2287@gmail.com](mailto:abdulmanan2287@gmail.com) | [LinkedIn](https://www.linkedin.com/in/abdul-manan-0aa546332/) | [GitHub](https://github.com/manan348)

**Last Updated:** November 2026

---

## Project Objective

> Can a standard GWAS pipeline, run entirely on public data, **recover genomic
> regions already known to confer blast resistance**?

The phenotype source paper reports resistance QTLs co-localizing with two cloned
resistance genes, *Pita* and *Ptr* (*Pita2*). Recovering an association peak over those
genes is the pipeline's correctness check. Note that this is a **re-analysis of the same
data** as the source paper, not an independent replication.

---

## Dataset Information

| | Genotype data | Phenotype data |
|---|---|---|
| **Source** | [RiceDiversity 44K SNP panel](http://www.ricediversity.org/data/sets/44kgwas/) | Lin et al. 2018, *Botanical Studies* 59:32 (Table S1) |
| **DOI** | — | [10.1186/s40529-018-0248-4](https://as-botanicalstudies.springeropen.com/articles/10.1186/s40529-018-0248-4) |
| **License** | Public research data (RiceDiversity.org) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |
| **Content** | 413 RDP1 accessions, 36,901 SNPs, positions on **MSU6** (IRGSP build 5) | Lesion type (LT) and diseased leaf area (DLA) for two isolates, D41-2 and 12YL-DL-3-2 |

Annotation inputs (also DVC-tracked under `data/annotation/`):
- MSU7 locus table `all.locus_brief_info.7.0` (gene models on IRGSP-1.0)
- RiceDiversity `44K.MSU6.SNP_Information.MSU7` — official per-SNP MSU6 → MSU7 coordinates

---

## Pipeline

```
parse_pheno → qc → dedup → pca → gemma → summarize → plots → gene_mapping
```

### 1 Phenotype extraction (`01_parse_phenotypes.py`)
Table S1 is parsed from the PDF with `pdfplumber` and a regex anchored on the
subpopulation codes. The script asserts the expected counts: **314** accessions parsed
(304 scored for D41-2, 228 for 12YL-DL-3-2), all 314 matching a genotyped NSFTV ID.

### 2 Genotype QC (`02_qc.py`, PLINK 1.9)
`--geno 0.1 --maf 0.05 --mind 0.1` → **30,122 SNPs, 383 accessions**.

### 3 Near-duplicate removal (`03_dedup.py`)
In fully inbred lines PLINK's IBD estimate (`PI_HAT`) saturates: 2,960 of 73,153 pairs
had `PI_HAT ≥ 0.99`, far too many to be real duplicates. Raw identity-by-state (`DST`)
is used instead. 56 pairs have `DST ≥ 0.99`; the 31 pairs where both accessions are
phenotyped form **18 clusters (41 accessions)**, e.g. two independent *Nipponbare*
submissions. The accession with the least missing data is kept per cluster (ties
broken by ID), removing **23** → **360 accessions** (267 phenotyped).
All clusters are listed in `results/tables/duplicate_clusters.csv`.

*Caveat:* near-duplicate pairs involving an unphenotyped accession are not removed;
they only affect the PCA, not the association tests.

### 4 Population structure (`04_pca.py`)
PCA on LD-pruned markers (`--indep-pairwise 50 5 0.2`, 1,314 SNPs). PC1 and PC2 explain
40.7 % and 19.1 % of the variance captured by the top 10 PCs; indica and aus cluster
together, japonica splits into temperate/tropical groups, and admixed accessions fall
between them.

### 5 Association testing (`05_gemma.py`, GEMMA 0.98.5)
Univariate LMM (`-lmm 4`, Wald test) per trait, with a centred kinship matrix
(`-gk 1`) and intercept + PC1–3 as covariates. MAF ≥ 0.05 and missingness ≤ 0.05 are
re-applied **within each trait's analysed subset** (`-maf 0.05 -miss 0.05`).

Thresholds are per trait, from the number of SNPs actually tested:
Bonferroni = 0.05 / n, suggestive = 1 / n (≈ 1.9 × 10⁻⁶ and ≈ 3.8 × 10⁻⁵).

| Trait | N | SNPs tested | SNP heritability (PVE ± SE) | λ_GC | Top p | Result |
|---|---|---|---|---|---|---|
| D41-2 Lesion Type | 257 | 26,059 | 0.748 ± 0.086 | 0.994 | 1.61 × 10⁻⁵ | suggestive only (6 SNPs) |
| D41-2 Diseased Leaf Area | 257 | 26,059 | 0.342 ± 0.173 | 0.892 | 1.34 × 10⁻⁶ | **Bonferroni** (1 SNP) |
| 12YL-DL-3-2 Lesion Type | 193 | 26,203 | 0.882 ± 0.062 | 0.911 | 3.56 × 10⁻⁷ | **Bonferroni** (9 SNPs) |
| 12YL-DL-3-2 Diseased Leaf Area | 193 | 26,203 | 0.295 ± 0.195 | 0.890 | 1.01 × 10⁻⁴ | **no association** |

Source: `results/tables/gwas_run_summary.csv`.

#### Observations
- **Two traits clear Bonferroni** (12YL-DL-3-2 LT and D41-2 DLA; the latter only just:
  1.34 × 10⁻⁶ vs 1.92 × 10⁻⁶). D41-2 LT reaches the suggestive level only, and
  12YL-DL-3-2 DLA has no SNP even at the suggestive threshold.
- **λ_GC is below 1** for three traits (0.89–0.91). This is *deflation*, not good
  control: the test statistics are conservative. Likely contributors are the small
  samples, the coarse ordinal/semi-quantitative scores, and an LMM that already
  absorbs much of the polygenic signal (high PVE). The practical consequence is lower
  power; no genomic-control correction was applied.
- LT has higher and more precisely estimated SNP heritability than DLA for both
  isolates, consistent with LT being the more direct readout of major-gene resistance.
  The PVE standard errors for DLA are large, so this comparison is indicative only.

### 6 Candidate gene mapping (`08_gene_mapping.py`)
The 31 suggestive-or-better SNPs are **lifted from MSU6 to MSU7** with RiceDiversity's
official per-SNP table before any distance is computed, so SNPs and gene models share
one coordinate system. (The shift between builds is ~2.6 kb on chr12, ~16 kb on chr10
and ~466 kb on chr11 at the hit SNPs — far from negligible.) Genes within ±200 kb are
recorded (`work/mapping/snp_gene_pairs.csv`); `top_candidate_genes.csv` lists the
nearest gene to each hit SNP plus the validation genes, with the strongest p-value per
gene. Genes annotated as transposon/retrotransposon-related are marked low-confidence.

** Validation: the chr12 peak covers *Pita* and *Ptr***
- ***Pita*** (`LOC_Os12g18360`, NB-ARC domain protein): two 12YL-DL-3-2 LT SNPs lie
  **inside** the gene (`id12004093`, p = 1.08 × 10⁻⁶, Bonferroni; `fd2`, p = 2.71 × 10⁻⁶).
- ***Ptr / Pita2*** (`LOC_Os12g18729`): 12YL-DL-3-2 LT SNPs 29–37 kb away
  (best p = 5.89 × 10⁻⁷, Bonferroni) and D41-2 LT SNPs 20–58 kb away.
- **Limits of this check.** All 12YL-DL-3-2 LT hits from *Pita* (10.6 Mb) to 12.4 Mb are
  in near-perfect LD with the lead SNP (r² = 0.94–1.00; `validation_pita_ptr_hits.csv`).
  The signal is therefore one long haplotype block in the low-recombination region
  around the chr12 centromere: it confirms that the known resistance locus is recovered,
  but the 44K panel **cannot tell *Pita*, *Ptr* or any other gene on the block apart**.
  The lead SNP itself (`id12004284`, 11.19 Mb) sits 0.9 kb from `LOC_Os12g19304`, which
  is best read as a passenger on the same haplotype, not a new candidate.

** Signals outside chr12 (single-trait, D41-2 DLA)**
- Chr11, `id11006398` (p = 1.34 × 10⁻⁶, the only Bonferroni hit outside chr12) lies
  inside `LOC_Os11g31330` (*no apical meristem* / NAC-family protein) on MSU7.
   An earlier version of this README named `LOC_Os11g30600` here; that came from
  comparing MSU6 SNP positions with MSU7 genes (a ~466 kb shift on this chromosome) and
  is retracted.
- Chr10 (3 suggestive SNPs, 6.2–7.1 Mb on MSU6): nearest genes include a protein kinase
  (`LOC_Os10g12620`, SNP inside the gene) and a retrotransposon-annotated gene
  (low confidence).

These are single-trait, single-SNP signals at or near the threshold and should be read
as leads for replication, not as discoveries.

---

## Outputs

`results/tables/`
| File | Content |
|---|---|
| `gwas_run_summary.csv` | Per-trait N, SNPs tested, PVE, λ_GC, thresholds, top p, significance class |
| `suggestive_hits_all_traits.csv` | All 31 SNPs with p ≤ 1/n (full-precision GEMMA output, MSU6 positions) |
| `top_candidate_genes.csv` | Nearest gene per hit SNP + validation genes (MSU7) |
| `validation_pita_ptr_hits.csv` | Hit SNPs within 200 kb of *Pita*/*Ptr*, with MSU6 and MSU7 positions and r² to the lead SNP |
| `duplicate_clusters.csv` | Near-duplicate clusters and which accession was kept |

`results/figures/` — 300-dpi PNGs (see `results/figures/README.md`).
`results/metrics.json` — pipeline-wide counts (`dvc metrics show`).
`notebooks/rice_gwas_walkthrough.ipynb` — read-only tour of the results.

---

## Limitations

- Modest sample sizes (193–257 per trait) and only 44K markers: limited power, and
  the chr12 signal cannot be fine-mapped below the ~2 Mb haplotype block.
- λ_GC < 1 for three traits: the tests are conservative and some weaker true signals
  are likely missed.
- LT and DLA are ordinal/semi-quantitative scores analysed with a Gaussian LMM.
- Validation re-uses the source paper's data, so it checks the pipeline, not the
  biology.
- Near-duplicates involving unphenotyped accessions are left in the PCA panel.

---

## Reproducing the analysis

Requires Linux, WSL or Colab (PLINK and GEMMA are Linux binaries).

```bash
# 1. Environment (either)
conda env create -f environment.yml && conda activate rice-gwas
# or: python -m venv .venv && . .venv/bin/activate && pip install -r requirements.txt

# 2. Pinned binaries (PLINK 1.90b7.2, GEMMA 0.98.5) into ./bin  (skip if using conda)
bash scripts/tools/get_binaries.sh

# 3. Data and pipeline
dvc pull          # genotype, phenotype PDF and annotation files from the DagsHub remote
dvc repro         # runs only the stages whose inputs or params changed
dvc metrics show
```

Change any setting in `params.yaml` (QC filters, number of PCs, thresholds, mapping
window, …) and re-run `dvc repro`; only the affected stages re-execute.

## Repository structure

```
├── data/                    # DVC-tracked raw inputs (data.dvc)
│   ├── genotype/            #   sativas413.ped/.map/.fam
│   ├── phenotype/           #   Lin et al. 2018 Table S1 (PDF)
│   └── annotation/          #   MSU7 locus table, MSU6→MSU7 SNP table
├── scripts/                 # 01–08 pipeline stages + gwas_utils.py, tools/get_binaries.sh
├── notebooks/               # results walkthrough
├── results/
│   ├── figures/             # 300-dpi Manhattan, QQ, PCA, chr12 figures
│   ├── tables/              # summary, hits, candidate genes, validation, duplicates
│   └── metrics.json
├── work/                    # DVC-cached intermediates (not in git)
├── dvc.yaml / dvc.lock      # pipeline definition and lock
├── params.yaml              # all analysis parameters
├── environment.yml / requirements.txt
└── README.md
```

## Data sources

- RiceDiversity 44K genotype data and SNP coordinate tables: http://www.ricediversity.org/data/sets/44kgwas/
- MSU Rice Genome Annotation Project, release 7: http://rice.uga.edu/
- Lin et al. 2018, *Botanical Studies* 59:32: https://as-botanicalstudies.springeropen.com/articles/10.1186/s40529-018-0248-4
