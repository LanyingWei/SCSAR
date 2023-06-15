# SCSAR Package

SCSAR is an R package specifically designed to facilitate the process of 
cell-type annotation for single-cell RNA-seq data analysis. It takes differential gene expression information for each cell cluster as input. Leveraging the 
comprehensive cell type marker gene data available in the 
[CellMarker 2.0](http://bio-bigdata.hrbmu.edu.cn/CellMarker/) and 
[Cell Taxonomy](https://ngdc.cncb.ac.cn/celltaxonomy/) databases, SCSAR employs 
a scoring algorithm to accurately and reliably assign cell-type annotations.

The Python version of SCSAR is [SCSA](https://github.com/bioinfo-ibms-pumc/SCSA). Compared to SCSA, SCSAR further includes Cell Taxonomy
as the marker gene resource, and enhances the scoring scheme by introducing a 
penalty system that penalizes cell types containing marker genes that are not 
overexpressed in the given cluster. This refinement leads to improved accuracy 
in cell-type identification. In addition, SCSAR provides the user with the top 
scoring markers for each cluster, enabling them to verify the results. Finally, 
SCSAR is more user-friendly for researchers performing single-cell sequencing 
analysis in R. Users no longer need to switch to Python for cell-type 
annotation, simplifying the workflow.


## Installation

To install the SCSAR package, you can use the `remotes` package:

``` r
if (!requireNamespace("remotes", quietly=TRUE))
    install.packages("remotes")
remotes::install_github("LanyingWei/SCSAR")
```

## Vignette

For more details, please refer to the 
[online vignette](https://htmlpreview.github.io/?https://github.com/LanyingWei/SCSAR/blob/main/vignettes/SCSAR.html).
