# DigestionSpecificity

An R package for analysis of cleavage site specificity (proteomics)

This package provides a set of functions that allows to plot cleavage site specificity based on identified peptide sequences. Both single amino acid level plot (barchart) and dipeptide-level specificity (heatmap) can be plotted.

To install this package, run the follwoings in R:

    library(devtools)
    install_github("ohgane/DigestionSpecificity")

If you don't have `devtools` yet, install it first by:

    install.packages("devtools")

After loading `DigestionSpecificity` package, you can check out a tutorial (package vignette) by:

    library(DigestionSpecificity)
    browseVignettes(package="DigestionSpecificity")
