# 🌾 Rice Blast Resistance GWAS

A genome-wide association study identifying genomic regions linked to blast disease
(*Magnaporthe oryzae*) resistance in the RDP1 rice diversity panel — a classical-
inference counterpart to **NeuroCrop**, a machine-learning maize yield predictor, in
the same portfolio.

---

👤 **Author:** Abdul Manan
🧪 *Plant Breeder | Machine & Deep Learning Researcher*
📧 [abdulmanan2287@gmail.com](mailto:abdulmanan2287@gmail.com) | 🔗 [LinkedIn](https://www.linkedin.com/in/abdul-manan-0aa546332/) | 💻 [GitHub](https://github.com/manan348)

🗓️ **Last Updated:** September 2026

---

## 🎯 Project Objective

> Can a standard GWAS pipeline, run entirely on public data, **recover genomic
> regions already known to confer blast resistance** — and along the way, surface
> novel candidates worth a closer look?

The pipeline was deliberately designed with a built-in validation checkpoint: the
phenotype source paper reports resistance QTLs co-localizing with two previously
cloned resistance genes, *Pita* and *Ptr*. Recovering hits near those genes from an
independently-run pipeline is the project's primary correctness check.

---

## 📂 Dataset Information

| | Genotype data | Phenotype data |
|---|---|---|
| **Source** | [RiceDiversity 44K SNP panel](http://www.ricediversity.org/data/sets/44kgwas/) | Lin et al. 2018, *Botanical Studies* 59:32 |
| **DOI** | — | [10.1186/s40529-018-0248-4](https://as-botanicalstudies.springeropen.com/articles/10.1186/s40529-018-0248-4) |
| **License** | Public research data (RiceDiversity.org) | [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/) |
| **Content** | 413 RDP1 accessions, MSU6 reference, ~36,900 SNP markers | Lesion type (LT) and diseased leaf area (DLA) scores for two *M. oryzae* isolates, D41-2 and 12YL-DL-3-2 |

---

## 🧬 Pipeline

### 1️⃣ Data acquisition and merge
Phenotype scores were extracted from the source paper's supplementary PDF (Table S1)
via `pdfplumber` and a regex parser anchored on rice subpopulation codes
(TEJ/IND/AUS/TRJ/ADMIX/AROMATIC), then matched to genotypes on NSFTV accession ID —
314 of 413 genotyped accessions had usable phenotype data.

### 2️⃣ Genotype QC (PLINK 1.9)
- Standard filters: `--geno 0.1 --maf 0.05 --mind 0.1` → 30,122 SNPs, 383 accessions
- **Duplicate/near-duplicate detection:** PLINK's default IBD estimator (`--genome`,
  `PI_HAT`) proved unreliable on this largely homozygous, self-pollinating species —
  it collapsed to `PI_HAT ≈ 1.0` for hundreds of genotypically dissimilar pairs.
  Switched to raw identity-by-state (`DST`), which produced a clean, interpretable
  distribution and correctly isolated 18 genuine near-duplicate clusters (41
  accessions total, including known cases like two independent *Nipponbare*
  submissions). Kept the least-missing-data accession per cluster, removing 23.
- **Final QC'd set: 360 accessions × 30,122 SNPs** (267 with usable phenotypes)
- PCA on LD-pruned markers confirmed the panel's population structure matches
  established rice genetics: indica and aus cluster together, japonica splits into
  temperate/tropical subgroups, admixed accessions scatter between clusters.

### 3️⃣ Association testing (GEMMA mixed model)
Univariate linear mixed models (`-lmm 4`) were run separately for all four
phenotype/isolate combinations, using a centered kinship matrix plus the first three
genotype PCs as covariates.

| Trait | N | PVE (heritability) | λ_GC | Top hit p-value |
|---|---|---|---|---|
| D41-2 Lesion Type | 257 | 0.748 (SE 0.086) | 0.993 | 1.61 × 10⁻⁵ |
| D41-2 Diseased Leaf Area | 257 | 0.342 (SE 0.173) | 0.883 | 1.34 × 10⁻⁶ |
| 12YL-DL-3-2 Lesion Type | 193 | 0.882 (SE 0.062) | 0.910 | 3.56 × 10⁻⁷ |
| 12YL-DL-3-2 Diseased Leaf Area | 193 | 0.295 (SE 0.195) | 0.892 | 1.01 × 10⁻⁴ |

#### 📈 Observations
- QQ plots for all four traits tracked the null diagonal closely, with lift confined
  to the extreme tail — indicating population structure was well-controlled by the
  kinship matrix and PC covariates, not left as residual inflation.
- Lesion type (a direct qualitative resistance readout) showed consistently higher
  and more precisely estimated heritability than diseased leaf area (a continuous,
  more environmentally-influenced severity measure) across **both** isolates —
  consistent with LT being the more direct readout of major-gene resistance.
- Two of four traits (12YL-DL-3-2 LT, D41-2 DLA) cleared strict Bonferroni
  significance; the rest are suggestive-threshold results only.

### 4️⃣ Candidate gene mapping
All 31 SNPs clearing a suggestive threshold (p ≤ 1/n_tests) were mapped against the
MSU7 gene annotation (genes within ±200 kb, to accommodate minor MSU6→MSU7
coordinate differences).

**✅ Validation confirmed on both target genes:**
- ***Pita*** (`LOC_Os12g18360`) — a 12YL-DL-3-2 lesion type SNP sits **1,392 bp**
  from this gene. Independently, the MSU annotation describes it as an "NB-ARC
  domain containing protein" — the defining NBS-LRR structural signature — matching
  Pita's known molecular function with no input from this analysis.
- ***Ptr/Pita2*** (`LOC_Os12g18729`) — hits from both isolates' lesion type traits
  land 23–35 kb away, consistent with the extended linkage disequilibrium typical of
  this low-recombination pericentromeric region.

**🔍 Other candidates:**
- `LOC_Os11g30600` — annotated only as a "hypothetical protein" (no known function),
  sitting 3.6 kb from the single strongest p-value in the full hit table
  (1.34 × 10⁻⁶, D41-2 diseased leaf area). Not a known resistance gene, but its
  statistical strength makes it worth flagging for future characterization.
- Three additional hits (chromosomes 10 and 11) map to transposon/retrotransposon-
  annotated genes. Repetitive regions like these are known to be more prone to
  genotyping and mapping artifacts, so these are reported as lower-confidence leads
  pending independent replication, not as credible candidate genes on their own.

---

## 📊 Visualizations

Manhattan and QQ plots for all four traits are saved in `results/figures/`:
- `manhattan_<trait>.png` / `qq_<trait>.png` for each of the four trait runs
- PCA scatter (PC1 vs PC2, colored by rice subpopulation)

All plots are 300-dpi, publication-ready exports.

---

## ⚠️ Limitations

- Only two of the four traits clear strict Bonferroni significance; the rest are
  suggestive-threshold results, appropriately weaker evidence than genome-wide
  significant hits.
- Gene mapping used the MSU7 annotation against MSU6-called SNP positions; a
  generous 200 kb window was used to absorb likely small coordinate discrepancies
  between builds, at some cost to mapping precision.
- Sample sizes per trait (193–267) are modest for GWAS by field standards, a
  consequence of incomplete phenotyping in the source dataset.

---

## ⚙️ Requirements

| Tool / Library | Version (recommended) |
|---|---|
| PLINK | 1.90b7.2 |
| GEMMA | 0.98.5 |
| Python | 3.8+ |
| pandas | latest |
| pdfplumber | latest |
| matplotlib | latest |
| scipy | latest |
| networkx | latest |
| Jupyter / Google Colab | either |

Install Python dependencies with:
```bash
pip install pandas pdfplumber matplotlib scipy networkx
```

---

## 📁 Repository structure

```
├── data/            # raw and intermediate data (genotype, phenotype, QC outputs)
├── scripts/         # analysis scripts / notebook
├── results/
│   ├── figures/     # Manhattan, QQ, and PCA plots
│   └── tables/      # association results, gene mapping tables
└── README.md
```

---

## ▶️ Reproducing the analysis

```bash
# Phase 2: QC
plink --file sativas413 --make-bed --out sativas413_binary
plink --bfile sativas413_binary --geno 0.1 --maf 0.05 --mind 0.1 --make-bed --out sativa_qc
plink --bfile sativa_qc --remove duplicates_to_remove.txt --make-bed --out sativa_qc_dedup
plink --bfile sativa_qc_dedup --pca 10 --out sativa_pca

# Phase 3: GWAS (GEMMA)
gemma -bfile sativa_qc_dedup -p gemma_phenotypes.txt -n 1 -gk 1 -o kinship
gemma -bfile sativa_qc_dedup -p gemma_phenotypes.txt -k kinship.cXX.txt \
      -c gemma_covariates.txt -n <trait_column> -lmm 4 -o <trait_name>
```

See `scripts/` for the full genotype-phenotype matching, QC, and gene-mapping code.

## 📚 Data sources

- RiceDiversity 44K genotype data: http://www.ricediversity.org/data/sets/44kgwas/
- Lin et al. 2018, *Botanical Studies* 59:32: https://as-botanicalstudies.springeropen.com/articles/10.1186/s40529-018-0248-4
