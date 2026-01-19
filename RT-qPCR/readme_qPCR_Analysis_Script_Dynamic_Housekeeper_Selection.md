# qPCR Analysis Pipeline with Dynamic Normalization

This script provides a robust framework for analyzing Quantitative PCR (qPCR) data, featuring automated reference gene selection and relative expression visualization.

##  Purpose

The script is designed to handle multi-sample qPCR experiments, specifically comparing mouse strains (B6N vs. BFMI) across different tissues. It ensures high data integrity by selecting the most stable internal controls for normalization.

##  Workflow

1.  **Stability Check:** Evaluates candidate housekeepers (GAPDH, Rpl13a, rps29) and selects the two most stable genes per dataset based on Standard Deviation.
2.  **Expression Calculation:** Performs $dCt$ calculations for a list of target genes (CX3RC1, Aldh1l1, Tubb3, MBP, BBS7).
3.  **Normalization:** Scales $dCt$ values into relative expression levels for standardized visual comparison.
4.  **Visualization:** Generates barplots with error bars. The Y-axis is fixed at 12 to allow direct visual comparison between different gene plots.
5.  **Export:** Saves calculated expression values and statistics into a CSV file for downstream analysis.

##  Requirements

### R Environment
* **Base R:** No external libraries are strictly required as the script uses standard `stats` and `graphics` packages.
* **Input Files:** Requires tab-separated text files (`.txt`) containing qPCR results with columns: `Sample Name`, `Target Name`, and `Cт Mean`.

##  Usage

1.  **Paths:** Update the `project_dir` in the script to your data location.
2.  **Gene List:** If necessary, adjust `target_genes` or `candidate_housekeepers` in the configuration section.
3.  **Run:** Execute the script to generate the summary CSV and the expression plots.
    ```bash
    Rscript qpcr_analysis.r
    ```

## ⚖ License & Attribution

* **Original Author:** Danny Arends (2021)
* **Modified by:** Florian Krause (2021-2025)
* **License:** GNU General Public License v3.0 (GPL-3.0)