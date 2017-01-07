#' Plot cleavage selectivity as a barchart from identified peptide sequences
#'
#' This function allows to examine specificity of enzymatic digestion (or chemical cleavage), based on MS/MS identification result (search without any enzyme specification, or "unspecific" search). Can be useful to check cleavage selectivity & efficiency of each enzyme used in multiple digestion or of less-specific enzymes.
#' @param seq A character vector of identified peptide sequences without any modifications. Should only made up of 20 amino acids one letter characters.
#' @param after A character vector of peptide sequences after cleavage site. Only the first amino acid is used, but can be of longer length.
#' @param hide.plot Logical, default FALSE. If set to TRUE, no plot will be drawn. Can be used just to obtain a data.frame of cleavage site frequency for manual plotting.
#' @param normalize Logical, default FALSE (stacked barchart of cleaved and not-cleaved).
#'
#' @return If (and only if) assigned to a variable, this function returns a data.frame with columns (\code{AA}, \code{Freq}, \code{terminal}, \code{cleavage}) for custom plotting.
#'
#' @seealso \code{\link{plotCleavageMatrix}}
#'
#' @export
#'
#' @examples
#' data(Elastase)
#' plotCleavage(Elastase$Sequence, Elastase$AAs.After)
#' plotCleavage(Elastase$Sequence, Elastase$AAs.After, normalize=TRUE)

plotCleavage=function(seq, after, hide.plot=FALSE, normalize=FALSE){
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
  # Concatenate identified peptide sequences as a character
  seqall=paste(seq, collapse="")
  # Extract last AA of identified sequences, and concatenate as a character like "KKRRFR...."
  left=paste(extractSeqLast(seq, n=1), collapse="")
  right=paste(extractAAsAfter(after, n=1), collapse="")
  right=gsub(" ", "", right, fixed=TRUE) # remove spaces (C-terminal fragment)

  # Extract AA composition
  aacLeft=protr::extractAAC(left)*nchar(left)
  aacRight=protr::extractAAC(right)*nchar(right)
  aacAll=protr::extractAAC(seqall)*nchar(seqall)

  ### Stacked barplot
  dfAll=data.frame(AA=names(aacAll), Freq=aacAll)
  # Left, cleaved
  dfLeft=data.frame(AA=names(aacLeft), Freq=aacLeft, cleavage="")
  dfLeft$terminal="Cleavage at C-term of"
  dfLeft$cleavage="cleaved"
  # Left, not cleaved
  dfLeftNot=dfLeft
  dfLeftNot$Freq=dfAll$Freq-dfLeft$Freq
  dfLeftNot$cleavage="not cleaved"
  # Right, cleaved
  dfRight=data.frame(AA=names(aacRight), Freq=aacRight, cleavage="")
  dfRight$terminal="Cleavage at N-term of"
  dfRight$cleavage="cleaved"
  # Right, not cleaved
  dfRightNot=dfRight
  dfRightNot$Freq=dfAll$Freq-dfRight$Freq
  dfRightNot$cleavage="not cleaved"
  # Combine
  AAsummary=rbind(dfLeft, dfLeftNot, dfRight, dfRightNot)
  rownames(AAsummary)=NULL
  if (normalize){
    AAsummary$Freq=AAsummary$Freq/rep(aacAll, 4)
  }

  if (!hide.plot & !normalize){
    p=lattice::barchart(Freq~AA | terminal, groups=AAsummary$cleavage, data=AAsummary,
               stack=TRUE, auto.key=list(columns=2),
               scales=list(tck=c(1,0), alternating=FALSE), origin=0,
               par.settings=lattice::standard.theme(col=FALSE),
               xlab="Amino acids", ylab="Frequency (counts)", ylim=c(0,NA), layout=c(2,1))
    print(p)
  } else if (!hide.plot & normalize){
    p=lattice::barchart(I(Freq*100)~AA | terminal,
                        data=subset(AAsummary, AAsummary$cleavage=="cleaved"),
                        stack=TRUE, auto.key=list(title="Cleavage frequency", columns=2),
                        scales=list(tck=c(1,0), alternating=FALSE), origin=0,
                        par.settings=lattice::standard.theme(col=FALSE),
                        xlab="Amino acids", ylab="Frequency (%)", ylim=c(-5,105), layout=c(2,1))
    print(p)
  }
  invisible(AAsummary)
}
