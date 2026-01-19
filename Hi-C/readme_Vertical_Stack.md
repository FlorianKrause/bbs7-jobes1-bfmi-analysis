# Hi-C Locus-Specific Heatmap Comparison

This R script is designed for the high-resolution visual comparison of Hi-C interaction matrices at a specific genomic locus across individual samples. It produces a vertical stack of heatmaps, allowing for the direct alignment of structural features across different replicates or strains.



##  Purpose

The script focuses on the **Bbs7 locus** (Chromosome 3) to investigate local chromatin architecture. By normalizing raw contact counts to the total library size, it enables a mathematically fair comparison between samples with different sequencing depths.

## 🛠 Features

* **Library Size Normalization:** Calculates the sum of all reads per `.hic` file and applies a scaling factor ($1.5 \times 10^7$) to ensure comparable signal intensities.
* **Vertical Stack Layout:** Generates a 3-panel vertical plot, facilitating the identification of boundary shifts or contact frequency changes along the genomic coordinates.
* **Locus Highlighting:** Overlays vertical reference lines (red, dashed) at the exact coordinates of the gene of interest (e.g., *Bbs7*).
* **Flexible Resolution:** Defaulted to **10kb resolution** for detailed inspection of local interactions.

##  Requirements

### R Packages
* `strawr`: Necessary for fast extraction of data from `.hic` files.

### Input Data
The script expects `.hic` files (using simplified naming conventions) located in the project directory:
* `BFMI1.hic`
* `B6N1.hic`
* `BFMI2.hic`

##  Usage

1.  **Configuration:** Update the `project_dir` path to point to your data directory.
2.  **Filename Adjustment:** If your files have different names, update the `files` vector in the configuration section.
3.  **Coordinate Setup:** The `val_start`, `val_end`, `sBBs`, and `eBBs` variables can be modified to target any other genomic region.
4.  **Run:**
    ```bash
    Rscript locus_comparison.r
    ```

## ️ License & Attribution

* **Original Author:** Danny Arends (2021)
* **Modified by:** Florian Krause (2025-2026)
* **License:** GNU General Public License v3.0 (GPL-3.0)