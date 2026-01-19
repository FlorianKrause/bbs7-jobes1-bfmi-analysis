# Hi-C Matrix Normalization & Visual Inspection

This script provides a specialized workflow for the high-resolution visualization and comparative analysis of Hi-C interaction matrices. It is optimized for the detailed inspection of the **Bbs7 locus** on Chromosome 3.



##  Purpose

The script facilitates a direct, normalized comparison between two mouse strains (**BFMI** vs. **B6N**) across multiple replicates. By calculating library sizes manually, it ensures that visual differences in interaction density are biologically significant and not artifacts of varying sequencing depths.

##  Workflow

1.  **Data Extraction:** Uses the `strawr` API to fetch raw interaction counts (no prior normalization) at **10kb resolution**.
2.  **Library Normalization:** * Calculates the total read sum for each `.hic` file.
    * Scales interaction frequencies based on a shared constant (`15,000,000`).
3.  **Region Filtering:** Focuses on a specific genomic window (default: 36.25 - 37.00 Mb on Chr 3).
4.  **Heatmap Generation:** * Constructs symmetric contact matrices.
    * Produces a 3x2 grid layout comparing BFMI and B6N replicates side-by-side.
    * Overlays the **Bbs7 gene locus** using red dashed reference lines.

##  Requirements

### R Packages
* `strawr`: Required for high-speed data access to `.hic` files.

### Input Files
The script expects `.hic` files named according to a simplified replicate structure in the project directory:
* `BFMI1.hic`, `BFMI2.hic`, `BFMI3.hic`
* `B6N1.hic`, `B6N2.hic`, `B6N3.hic`

##  Usage

1.  **Configuration:** Update the `project_dir` variable in the script to the location of your `.hic` files.
2.  **Locus Adjustment:**
    The coordinates for `val_start`, `val_end`, and the `sBBs`/`eBBs` markers can be adjusted to inspect other genomic regions of interest.
3.  **Run:**
    ```bash
    Rscript matrix_visualization.r
    ```

##  License & Attribution

* **Original Author:** Danny Arends (2021)
* **Modified by:** Florian Krause (2021-2025)
* **License:** GNU General Public License v3.0 (GPL-3.0)