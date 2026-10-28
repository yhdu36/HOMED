# Overview
HOMED (Hierarchically Optimized Methylation Deconvolution) is a reference-based framework for estimating cell-type composition from bulk DNA methylation (DNAm) data. Unlike conventional methylation deconvolution approaches that treat cell types as flat entities, HOMED explicitly incorporates hierarchical relationships among related cell populations to improve resolution of closely related subtypes.

The framework combines purified FACS/MACS-derived methylation references with pseudo-ground-truth cell proportions generated through scRNA-seq-guided deconvolution of paired bulk RNA-seq samples. HOMED then performs hierarchical differential methylation analysis and iterative IDOL-based optimization to identify CpG libraries that maximize deconvolution accuracy. The resulting hierarchical reference matrices can be applied to new bulk methylation datasets to estimate both major cell lineages and finer cellular subtypes while maintaining biological consistency across hierarchical levels.

HOMED is designed for tissue-specific methylation deconvolution and can be applied to a wide range of tissues where purified methylation references and single-cell transcriptomic atlases are available.

A complete PBMC tutorial is available here: 

- [PBMC Vignette](https://github.com/yhdu36/HOMED/blob/main/vignettes/HOMED_PBMC.html)

# Option 1: Install as an R package

```r
library(devtools)

install_github("yhdu36/HOMED")
```

### Option 2: Docker

HOMED is available as a prebuilt Docker image with R and all required dependencies installed. No local R or Bioconductor installation is needed.

**1. Install Docker**

Install [Docker Desktop](https://www.docker.com/products/docker-desktop/) (Windows/macOS) or Docker Engine (Linux).

**2. Download the HOMED image**

```bash
docker pull --platform linux/amd64 ghcr.io/yhdu36/homed:0.1.0
```

**3. Launch HOMED**

```bash
docker run --rm -it --platform linux/amd64 ghcr.io/yhdu36/homed:0.1.0
```

Inside the R console:

```r
library(HOMED)
packageVersion("HOMED")
```

**4. Work with local datasets (optional)**

Create `data` and `results` folders, then mount them into the container:

```bash
mkdir -p data results

docker run --rm -it --platform linux/amd64 \
  -v "$(pwd)/data:/data:ro" \
  -v "$(pwd)/results:/results" \
  ghcr.io/yhdu36/homed:0.1.0
```

Input files are accessible at `/data`, and results saved to `/results` will appear on your computer. The commands above use a Unix-style shell; Windows users can use WSL or Git Bash.

Refer to the [PBMC Vignette](https://github.com/yhdu36/HOMED/blob/main/vignettes/HOMED_PBMC.html) for example analyses.

# References
Liu, Y., Du, Y., Chen, Y., and Garmire, L. X. (2026). HOMED: Hierarchically Optimized Methylation Deconvolution for Accurate Cell-Type Estimation from Bulk DNA Methylation Data. bioRxiv. https://doi.org/10.64898/2026.06.08.730896

Du Y, Benny PA, Lahiri S, AlAkwaa FM, Huang Q, Liu Y, Lassiter CB, Astern J, Riel J, Garmire LX (2026). Placental Molecular Subtypes of Severe Preeclampsia Reveal Divergent Aging Trajectories and Fetal Growth Outcomes. medRxiv. https://doi.org/10.64898/2026.06.02.26354756.



