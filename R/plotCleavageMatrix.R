#' Plot cleavage site matrix from identified peptide sequences
#'
#' This function allows to examine specificity of enzymatic digestion (or chemical cleavage) at dipeptide sequence level, based on MS/MS identification result (search without any enzyme specification, or "unspecific" search). Can be useful to check cleavage selectivity & efficiency of each enzyme used in multiple digestion or of less-specific enzymes.
#' @param seq A character vector of identified peptide sequences without any modifications. Should only made up of 20 amino acids one letter characters.
#' @param after A character vector of peptide sequences after cleavage site. Only the first amino acid is used, but can be of longer length.
#' @param hide.plot Logical, default FALSE. If set to TRUE, no plot will be drawn. Can be used just to obtain a data.frame of cleavage site frequency for manual plotting.
#' @param normalize Logical, default FALSE (stacked barchart of cleaved and not-cleaved).
#' @param col.regions Color palette function like \code{heat.colors} etc.
#'
#' @return If assigned to a variable, this function returns a cleavage site matrix, with its row being the first amino acid and its column being the second amino acid.
#'
#' @seealso \code{\link{plotCleavage}}, \code{\link{extractDCasMat}}
#'
#' @export
#'
#' @examples
#' data(LysCTryp)
#' plotCleavageMatrix(LysCTryp$Sequence, LysCTryp$AAs.After, normalize=TRUE)

plotCleavageMatrix=function(seq, after, hide.plot=FALSE, normalize=FALSE, col.regions=NULL){
  # Define two utility functions
  extractAAsAfter=function(char, n=1){
    f=function(x){
      res=substr(x, 1, n)
      if (nchar(res)==0){
        res=paste(rep(" ", n), collapse="")
      }
      return(res)
    }
    res=as.vector(sapply(char, f)) # Extract first 2 amino acids (from "CC; CC; CC" etc)
    return(res)
  }
  extractSeqLast=function(char, n=1){
    f=function(x) substr(x, nchar(x)-(n-1),nchar(x)) # extract last two
    res=as.vector(sapply(char, f)) # Vectorize
    return(res)
  }

  # Extract cleavage site dipeptide
  clsite=paste(extractSeqLast(seq, n=1), extractAAsAfter(after, n=1), sep="")
  # remove entries containing spaces (C-terminal fragment like "R " causes error in extractDC)
  clmat=extractDCasMat(clsite[!grepl(" ", clsite)], count=TRUE)

  # Extract dipeptide frequency from sequences
  firstafter=extractAAsAfter(after, n=1)
  firstafter=gsub(" ", "", firstafter, fixed=TRUE)
  theomat=extractDCasMat(paste(seq, firstafter, sep=""), count=TRUE)

  # Normalize cleavage site count
  if(normalize){
    clmat=clmat/theomat
    clmat[is.nan(clmat)]=0
  }

  # Define a default palette

  if (is.null(col.regions)){
    pal_magma=grDevices::colorRampPalette(c("white",
                                            rev(viridis::magma(12))), bias = 1)
    col.regions=pal_magma(256)
  }

  # Plot cleavage site matrix
  if (!hide.plot){
    p=lattice::levelplot(clmat, col.regions=col.regions,
                         scales=list(tck=c(1,0)), cuts=256,
                         xlab="N-term of cleavage site",
                         ylab="C-term of cleavage site")
    print(p)
  }
  invisible(clmat)
}
