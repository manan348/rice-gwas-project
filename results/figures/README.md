# Figures

All figures are written by the DVC pipeline (`scripts/04_pca.py`, `07_plots.py`,
`08_gene_mapping.py`) at 300 dpi.

| File | Description |
|---|---|
| `pca_subpopulation.png` | PC1 vs PC2 of LD-pruned markers for all 360 accessions, coloured by subpopulation (unlabelled accessions in grey). |
| `manhattan_d41_2_LT.png` | Manhattan — D41-2 Lesion Type (suggestive only) |
| `manhattan_d41_2_DLA.png` | Manhattan — D41-2 Diseased Leaf Area (**clears Bonferroni**, chr11) |
| `manhattan_12yl_LT.png` | Manhattan — 12YL-DL-3-2 Lesion Type (**clears Bonferroni**, chr12 peak) |
| `manhattan_12yl_DLA.png` | Manhattan — 12YL-DL-3-2 Diseased Leaf Area (no association) |
| `qq_*.png` | QQ plot per trait with λ_GC; three traits are deflated (λ ≈ 0.89–0.91) |
| `manhattan_all_traits_panel.png` | Four-trait stacked Manhattan panel — chr12 signal in both LT traits |
| `hits_summary_by_chromosome.png` | Suggestive-or-better SNPs per chromosome and trait |
| `candidate_genes_chr12.png` | Chr12 region (MSU7 coordinates) with *Pita* and *Ptr* positions marked |

Dashed red line = Bonferroni (0.05/n); dotted blue = suggestive (1/n), n = SNPs tested
for that trait. Manhattan/QQ positions are MSU6 (the genotype build); the chr12 figure
uses SNPs lifted to MSU7 so they line up with the gene models.
