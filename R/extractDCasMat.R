#' Calculate dipeptide frequency matrix
#'
#' This function calculates a dipeptide frequency (or count) matrix from a vector of sequences.
#' @param seqvec A character vector of peptide/protein sequences.
#' @param count Logical, default TRUE. If TRUE, returns total count of dipeptide. If FALSE, returns frequency.
#'
#' @return A matrix of dipeptide frequency/count, with its row being the first amino acid and its column being the second amino acid.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # An artificial, very short example
#' seq=c("AHGTYIR", "VCSDWEQLPMK", "FGGNTIPCAWMR", "QTYILMNK", "ASTEMNDNK")
#' dipepmat=extractDCasMat(seq)
#' print(dipepmat) # Note that the presence of two "NK" motifs
#' lattice::levelplot(dipepmat, xlab="First", ylab="Second", col.regions=heat.colors)
#' }

extractDCasMat=function(seqvec, count=TRUE){
  #For each sequences, count AAs (calculate frequency, and convert to count)
  countlist=round(t(sapply(seqvec, protr::extractDC))*sapply(seqvec, nchar))
  #Combine all the sequences by summing
  theo=apply(t(countlist), 1, sum)
  if (count){#return count matrix
    mat=matrix(theo, 20, 20)
  } else {#return frequency matrix
    mat=matrix(theo/sum(theo), 20, 20)
  }
  aalist=c("A","R","N","D","C","E","Q","G","H","I","L","K","M","F","P","S","T","W","Y","V")
  rownames(mat)=aalist # First AA
  colnames(mat)=aalist # Seconod AA
  return(mat)
}
