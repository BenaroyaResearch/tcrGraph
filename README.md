# tcrGraph - Graph-based Analysis of T Cell Receptors

A package to facilitate building and analyzing networks of T cell receptors to help better understand clonality, antigen specificity, and other immunological questions.

### Installation

Open an R session, and execute the following code:

```R
if (!require("devtools")) install.packages("devtools")
devtools::install_github("BenaroyaResearch/tcrGraph")
```

## Requirements

The following R packages will be automatically installed if necessary: `stringr, dplyr, magrittr, rlang, igraph, httr`

## Vignettes

This package includes an introductory vignette demonstrating basic use of the package functions. To view this vignette, execute the following R code:

```R
vignette("tcrGraphIntro", "tcrGraph")
```

## Supported Functions and Classes

* `tcrGraph`
  * An S3 class, with implementations of methods `print`, `plot`, and `is.tcrGraph`. Generally you should create a `tcrGraph` object using the function `makeTcrGraph()`.

* `makeTcrGraph(tcrDf, link = "full_nt_sequence")`
  * Accepts a data frame of TCR information, with columns for libid, v_gene, j_gene, and the column indicated by the link argument. Returns a tcrGraph object.
  
* `getClonesFromTcrGraph(tcrData, maxA = 2, maxB = 2, maxD = -1, maxG = -1, format = "compressed", link = "full_nt_sequence")`
  * Accepts either a tcrGraph object or a data frame that can be used to construct a tcrGraph object. Generates summary information about the frequencies of the clones in the dataset. A 'clone' is defined as a maximal subgraph containing at most the number of alpha, beta, delta, and gamma chains set by `maxA`, `maxB`, `maxD`, and `maxG`.

