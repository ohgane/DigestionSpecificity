#' Peptide identification results of HeLa cells digested with Elastase
#'
#' Identification results from a bottom-up proteomics experiment. HeLa cell lysate was digested with Elastase using TFE digestion protocol on an AmiconUltra 10K, and the digest was analyzed on a Thermo LTQ-Orbitrap XL in data-dependent acquisition mode (orbitrap detector for MS1 at 60000 resolution, and ion trap for MS2). The obtained data was searched against a subset of human protein database (from uniprot, ~1500 sequences) by using SearchGUI (X!Tandem and OMSSA, unspecific search, precursor tolerance 4 ppm, and fragment ion tolerance 0.4 Da). The identification results were processed with PeptideShaker (FDR 1%) and exported to xls file (as "Default Peptide Report").
#'
#' @docType data
#'
#' @usage data(Elastase)
#'
#' @format An object of class \code{data.frame}, with its rows corresponding to each peptide identifications. Its columns include \code{Sequence}, \code{AAs.After} etc.
#'
#' @keywords datasets
#'
#' @examples
#' data(Elastase)
#' head(Elastase)
"Elastase"
