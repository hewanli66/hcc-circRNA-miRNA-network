# circRNA–miRNA–mRNA network in hepatocellular carcinoma — reproducible archive

Analysis archive for **"Cross-cohort circRNA–miRNA–mRNA network analysis identifies
NT5DC2 and HSD17B6 as prognostic markers in hepatocellular carcinoma"**.

**Author**: Wanli He — School of Medicine, Shangqiu Institute of Technology
**ORCID**: 0009-0002-0829-1909
**Version**: 1.0.0 · 2026-10-11
**License**: CC BY 4.0

## What this archive contains

| Directory | Contents |
|---|---|
| `results/` | All result tables underlying every number reported in the manuscript |
| `scripts/` | Complete analysis pipeline, in execution order |
| `figures/` | Figure source files (SVG) and 600 dpi submission renders (PNG) |
| `manuscript/` | Manuscript, cover letter and supplementary information |

## Public data are NOT redistributed here

Obtain the primary data from their original repositories using the accessions below.

| Dataset | Accession | Source |
|---|---|---|
| GSE122482 | GSE122482 | NCBI GEO |
| GSE14520 | GSE14520 | NCBI GEO |
| TCGA-LIHC | — | Genomic Data Commons |
| circBank | www.circbank.cn | circBank database |

## Execution order

```
 1. step1_map_filter.py           differential circRNA screening, host gene mapping
 2. step2_tcga_survival.py        COX / log-rank in TCGA-LIHC
 3. step3_validate14520.py        independent validation in GSE14520
 4. step4_immune.py               immune infiltration scoring and correlation
 5. step6_network.py              circRNA–miRNA–mRNA network assembly
 6. step7_miRNA_binding_v2.py     circBank retrieval, coordinate harmonisation, edges
 7. step8_network_qc.py           deduplication, quality scoring, tier assignment
 8. step9_nt5dc2_support.py       indirect support for NT5DC2
 9. step9b_fallback_support.py    indirect support for all fallback-matched genes
10. make_figs3.pyfigure generation
11. make_submission_figures.py    600 dpi render for submission
```

## Environment

Python 3.13.12. **All analyses use the standard library only** and therefore require
no third-party packages. `make_figs3.py` and `make_submission_figures.py` additionally
require `Pillow`; `make_submission_figures.py` uses a headless browser to rasterise SVG.

## Two technical issues encountered and resolved

1. **Chromosome nomenclature mismatch.** GSE122482 uses `chrN:start-end` while circBank
   uses `N:start-end`. Without harmonisation only 33 of 103 candidates matched (32%);
   after harmonisation 95 matched (92%).
2. **Duplicate edges from pagination overlap.** The circBank binding interface returns
   overlapping records at pagination boundaries, producing 766 duplicate edges (5.0%).
   Removed by composite key on (circRNA, host gene, miRNA, energy, score), leaving
   14,568 unique edges.

## Citation

If this archive is used, please cite the manuscript above.
