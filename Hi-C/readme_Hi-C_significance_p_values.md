# Hi-C Significance Heatmap (ggplot2)

This script provides a high-level visualization of pre-calculated significant chromosomal interactions using the `ggplot2` framework. It is specifically designed to render a clean, publication-ready heatmap of differential interaction p-values between **B6N** and **BFMI** mouse strains.

##  Purpose

While the primary analysis script handles raw data extraction from `.hic` files, this script is dedicated to:

* **Visualizing Significance:** Mapping negative log-transformed p-values () onto a genomic grid to highlight statistical differences.
* **Structural Highlighting:** Identifying and emphasizing specific "key structure" bins (marked with a blue border) that constitute critical architectural differences in the focus region.
* **Symmetry Rendering:** Displaying interactions across the diagonal to provide a standard Hi-C map perspective.

---

##  Features

* **Focused Resolution:** Optimized for a **25 kb** resolution grid.
* **Target Region:** Precision plotting for Chromosome 3 (Focus: 35.5 - 38.0 Mb).
* **Custom Styling:** Uses a specialized color gradient from light grey (low significance) to intense red (high significance).
* **Fixed Aspect Ratio:** Utilizes `coord_fixed()` to ensure the plot remains perfectly square, maintaining accurate spatial representation of genomic coordinates.

---

##  Requirements

### R Environment

The script requires the following CRAN packages:

* `ggplot2`: For advanced graphics and heatmap rendering.
* `dplyr`: For efficient data transformation and filtering.

If not present, the script will attempt to install them automatically:

```r
install.packages(c("ggplot2", "dplyr"))

```

---

##  Usage

1. **Modify Configuration:**
Open the script and adjust the `USER CONFIGURATION` section:
* `plot_min` / `plot_max`: Change the genomic axis limits (in Mb).
* `key_p_values`: Define the specific p-values used to trigger the highlight borders for structural bins.
* `highlight_color`: Change the color of the structural borders (default is blue).


2. **Execution:**
Run the script within your R IDE or via the terminal:
```bash
Rscript significance_plot.r

```



---

##  Output

The script generates a ggplot object that:

1. Plots the interactions as tiles.
2. Adds a dashed symmetry line.
3. Draws blue rectangles around bins identified as part of the "key structure" (e.g., specific interactions around the *Bbs7* locus).

---

## ⚖️ License & Attribution

### Scripts

Licensed under the **GNU General Public License v3.0 (GPL-3.0)**.

* **Author:** Florian Krause (2025)

### Dependencies

* **ggplot2:** H. Wickham (MIT)
* **dplyr:** Hadley Wickham et al. (MIT)
