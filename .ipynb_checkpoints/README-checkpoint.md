# Overview
HOMED (Hierarchically Optimized Methylation Deconvolution) is a reference-based framework for estimating cell-type composition from bulk DNA methylation (DNAm) data. Unlike conventional methylation deconvolution approaches that treat cell types as flat entities, HOMED explicitly incorporates hierarchical relationships among related cell populations to improve resolution of closely related subtypes.

The framework combines purified FACS/MACS-derived methylation references with pseudo-ground-truth cell proportions generated through scRNA-seq-guided deconvolution of paired bulk RNA-seq samples. HOMED then performs hierarchical differential methylation analysis and iterative IDOL-based optimization to identify CpG libraries that maximize deconvolution accuracy. The resulting hierarchical reference matrices can be applied to new bulk methylation datasets to estimate both major cell lineages and finer cellular subtypes while maintaining biological consistency across hierarchical levels.

HOMED is designed for tissue-specific methylation deconvolution and can be applied to a wide range of tissues where purified methylation references and single-cell transcriptomic atlases are available.

A complete PBMC tutorial is available here: 

- [PBMC Vignette](https://github.com/yhdu36/HOMED/blob/main/vignettes/HOMED_PBMC.html)

# Installation

```r
library(devtools)

install_github("yhdu36/HOMED")
```

# References
Du, Yuheng, Paula A. Benny, Shayanki Lahiri, Fadhl M. AlAkwaa, Qianhui Huang, Yuansen Liu, Cameron B. Lassiter, Joshua Astern, Jonathan Riel, and Lana X. Garmire. 2026. Placental Molecular Subtypes of Severe Preeclampsia Reveal Divergent Aging Trajectories and Fetal Growth Outcomes. medRxiv. medRxiv. https://doi.org/10.64898/2026.06.02.26354756.


