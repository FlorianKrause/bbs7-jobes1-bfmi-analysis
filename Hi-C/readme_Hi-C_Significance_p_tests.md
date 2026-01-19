# Hi-C Interaction Analysis & Differential Plotting

This directory contains R-based tools for analyzing and visualizing chromosomal interactions from **.hic** files. The pipeline is designed to perform differential interaction analysis between two mouse strains (**B6N** and **BFMI**) and generate high-resolution heatmaps showing statistical significance (p-values) of these contacts.

##  Contents

* **hic_analysis.r**: The primary R script for data extraction, statistical testing (t-tests), and visualization of Hi-C contact frequencies.
* **[Data Files]**: The script processes `.allValidPairs.hic` files (KR normalized) for replicates of B6N and BFMI strains.

---

## Analysis Workflow

### 1. Data Retrieval & Annotation

The script integrates genomic annotations to provide biological context to the interaction maps:

* **Ensembl biomaRt:** Retrieves protein-coding gene information for the target region (Chromosome 3).
* **AnnotationHub (CTCF):** Pulls mm10-specific **CTCF binding site** data to identify potential anchors of chromatin loops or Topologically Associating Domains (TADs).

### 2. Hi-C Data Extraction (`strawr`)

Using the **Knight-Ruiz (KR)** normalization via the `straw` API, the script extracts contact matrices for:

* **Resolution:** Adjustable bin sizes (25kb, 50kb, or 100kb).
* **Region:** Specifically focused on a critical region on **Chromosome 3** (approx. 30Mb to 40Mb).

### 3. Statistical Comparison

Instead of simple frequency plotting, this pipeline focuses on **differential interactions**:

* **Replicate Processing:** Aggregates data from 3 replicates per strain.
* **T-Tests:** Performs a bin-wise t-test to compare the interaction strength between B6N and BFMI replicates.
* **Transformation:** Converts p-values to a `-log10(p-value)` scale for heatmap visualization.

### 4. Visualization

The final output is a specialized interaction plot:

* **Heatmap:** Represents the significance of interaction differences (P-values).
* **Gene Overlays:** Automatically marks the positions of key genes such as *Bbs7*, *Exosc9*, *Ccna2*, *Trpc3*, and *4932438A13Rik*.
* **Structural Markers:** Highlights the **Bbs7 CTCF site** and specific interaction anchors in red.

---

## 📋 Requirements

### R Environment

The following packages are required for the analysis:

* `strawr`: Fast data extraction from .hic files.
* `biomaRt`: Genomic annotation retrieval.
* `AnnotationHub` & `GenomicRanges`: Managing chromosomal coordinates and external datasets.
* `RColorBrewer`: Color scales for heatmaps.

### Input Data

Ensure your `.hic` files are named according to the replicate structure (e.g., `FKR02_1bfmi_L001.allValidPairs.hic`) and located in the directory specified in the configuration.

---

## 🚀 Usage

1. **Configuration:** Open the script and update the `USER CONFIGURATION` section:
* `setwd()`: Set the path to your data.
* `bin_size`: Choose your desired resolution.
* `pS / pE`: Define the start and end positions (in bp) for the plot.


2. **Execution:**
Run the script within R or via terminal:
```bash
Rscript hic_analysis.r

```



---

## ⚖️ License & Attribution

### Scripts

Licensed under the **GNU General Public License v3.0 (GPL-3.0)**.

* **Original Author:** Danny Arends
* **Modified by:** Florian Krause (2021-2025)

### Dependencies & Data

* **strawr:** Aiden Lab (MIT)
* **Annotation Data:** Ensembl (Archive Nov 2020) and AnnotationHub (mm10).
