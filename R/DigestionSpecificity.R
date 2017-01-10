#' Analyze Cleavage Site Specificity
#'
#' This package provides a set of functions that allows to plot
#' cleavage site specificity based on identified peptide sequences. Both single
#' amino acid level plot (barchart) and dipeptide-level specificity (heatmap) can
#' be plotted. For more detail, see the package vignette with \code{vignette(package="DigestionSpecificity")}.
#'
#' @seealso \code{\link{plotCleavage}}, \code{\link{plotCleavageMatrix}}
#'
#' @examples
#' data(LysCTryp)
#' plotCleavage(LysCTryp$Sequence, LysCTryp$AAs.After)
"_PACKAGE"
