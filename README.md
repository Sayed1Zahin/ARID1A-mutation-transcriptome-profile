# Differential Gene Expression Analysis of GSE54979

**Targeting EZH2 in ARID1A-mutated Ovarian Clear Cell Carcinoma**

## 📌 Project Overview

This repository contains the complete workflow of data analysis and visualization for the microarray dataset [GSE54979](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54979), which explores gene expression changes resulting from ARID1A restoration or EZH2 inhibition (via GSK126) in ovarian clear cell carcinoma (OCCC) cells.

The project demonstrates a full pipeline using R for:

* Data preprocessing
* Differential expression analysis using the `limma` package
* Visualizations including volcano plots, MA plots, PCA, heatmaps, and more

## 🧬 Background

ARID1A, a chromatin remodeler, is frequently mutated in OCCC. Meanwhile, EZH2, a histone methyltransferase, is often upregulated. In this study, GSK126 (an EZH2 inhibitor) was used to investigate synthetic lethality in ARID1A-mutated cells.

**Experimental Design:**

* **Samples:** ARID1A-mutated cells, control cells, and cells treated with EZH2 inhibitor (GSK126)
* **Platform:** Expression profiling by array (Homo sapiens)
* **Goal:** Identify gene expression changes associated with ARID1A mutation and EZH2 inhibition, and uncover potential synthetic lethal interactions.

## 🛠️ Tools and Packages Used

* R (version ≥ 4.0)
* `limma`
* `Biobase`
* `affy`
* `ggplot2`
* `ggrepel`
* `pheatmap`
* `FactoMineR`, `factoextra` (for PCA)
* `EnhancedVolcano` (optional alternative)

## 📂 Repository Structure

```
├── Data/                   # Raw and processed expression data
│   ├── GSE54979_RAW/       # Raw CEL files or normalized expression matrix
│   ├── design_matrix.csv   # Sample metadata
│   └── *.csv               # Results from DE analysis
├── Figures/                # All generated plots
│   ├── volcano_ARID1A.png
│   ├── MA_plot.png
│   └── ...
├── Scripts/
│   ├── 01_preprocessing.R
│   ├── 02_differential_expression.R
│   └── 03_visualizations.R
├── README.md
└── sessionInfo.txt         # R session info for reproducibility
```

## 📊 Analyses Performed

* **Differential Expression:**

  * ARID1A vs Control
  * GSK126 vs Control
  * GSK126 vs ARID1A
* **Visualizations:**

  * Volcano plots
  * MA plots
  * Principal Component Analysis (PCA)
  * Heatmaps of top DE genes
  * Histogram of adjusted p-values

## 📈 Example: Volcano Plot

![Volcano Plot](Figures/volcano_ARID1A.png)

## ✅ How to Reproduce

1. Clone this repository:

   ```bash
   git clone https://github.com/yourusername/GSE54979_DE_analysis.git
   cd GSE54979_DE_analysis
   ```

2. Open R or RStudio and run the scripts in order:

   * `Scripts/01_preprocessing.R`
   * `Scripts/02_differential_expression.R`
   * `Scripts/03_visualizations.R`

3. All output tables and figures will be saved under `/Data` and `/Figures`.

## 📚 Citation

Original study:
Bitler, B. G., et al. (2015). *Targeting EZH2 methyltransferase activity in ARID1A mutated cells as a synthetic lethal therapeutic strategy*. **[GSE54979](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE54979)**.
