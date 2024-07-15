
#' Title Obtain DNA methylation beta value matrix for all miRNAs
#'
#' @param beta The DNA methylation beta value matrix with rownames annotated to
#' Illumina array probe IDs.
#' @param array array A parameter specifying the type of input data matrix, it should be
#' one of the "450k" for Illumina 450k array, or "EPIC" for Illumina EPIC (850k)
#' array, or "EPICv2" for Illumina EPICv2 (935k) array, default for "450k".
#' @param IlmnID A logical to indicate whether the IlmnID (which comprise the 'Name'
#' and a four-character suffix, e.g. 'cg09617579_BC12') was used for EPICv2 data, default
#' for FALSE.
#' @param type A parameter specifying the miRNA type used. it should be either
#' "mature" for mature miRNA, or "pre" for precursor miRNA, default for "mature".
#' @param upstream An integer specifying upstream boundary of the promoter region
#' from transcription start site, the maximum is 5000, default for 1500.
#' @param downstream An integer specifying downstream boundary of the promoter region
#' from transcription start site, the maximum is 5000, default for 500.
#' @param method A parameter specifying how the methylation level for a particular
#' miRNA be calculated. "mean" for that average beta value of all CpGs within
#' promoter region, "nearest" for that the beta value of the nearest CpG from
#' the transcription start site as representative, "loess" for taht the loess
#' fitted average beta value of all CpGs within promoter region.
#'
#' @return DNA methylation matrix with rows are miRNAs columns are samples
#' @export
#'
#' @examples
#' data(beta450k.m)
#' miRBeta.m <- probe2miR(beta=beta450k.m, array="450k", type="Mature", upstream=2000, downstream=2000, method="mean")
#'
probe2miR <- function(beta, array="450k", IlmnID=FALSE, type="Mature", upstream=1500, downstream=500, method="mean"){
  up <- upstream
  down <- downstream
  method <- method
  if(array=="450k"){
    if(type=="Mature"){
      betaGene <- lapply(probe450kPosAnno2Mature, function(x) getGeneBeta(x, beta, upstream=up, downstream=down, method=method ) )
    }else if(type=="Pre"){
      betaGene <- lapply(probe450kPosAnno2Pre, function(x) getGeneBeta(x, beta, upstream=up, downstream=down, method=method) )
    }
  }else if(arrayType=="EPIC"){
    if(type=="Mature"){
      betaGene <- lapply(probeEPICPosAnno2Mature, function(x) getGeneBeta(x, beta, upstream=up, downstream=down, method=method) )
    }else if(type=="Pre"){
      betaGene <- lapply(probeEPICPosAnno2Pre, function(x) getGeneBeta(x, beta, upstream=up, downstream=down, method=method) )
    }
  }else if(arrayType=="EPICv2"){
    if(type=="Mature"){
      if(IlmnID == FALSE){
        betaGene <- lapply(probeEPICv2PosAnno2Mature, function(x) getGeneBeta(x, beta, upstream=up, downstream=down, method=method) )
      }else if(IlmnID == TRUE){
        betaGene <- lapply(probeEPICv2PosAnno2Mature_IlmnID, function(x) getGeneBeta(x, beta, upstream=up, downstream=down, method=method) )
      }
    }else if(type=="Pre"){
      if(IlmnID == FALSE){
        betaGene <- lapply(probeEPICv2PosAnno2Pre, function(x) getGeneBeta(x, beta, upstream=up, downstream=down, method=method) )
      }else if(IlmnID == TRUE){
        betaGene <- lapply(probeEPICv2PosAnno2Pre_IlmnID, function(x) getGeneBeta(x, beta, upstream=up, downstream=down, method=method) )
      }
    }
  }
  data.m <- do.call(rbind, betaGene)
}
