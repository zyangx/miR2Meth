#' Title Obtain beta value for each individual miRNA
#'
#' @param pos A vector consisting of the relative positions of the CpGs within
#' promoter region of a miRNA to its transcription start site, the names of
#' the elements in the vector should be the CpG ID of the Illumina Infinium
#' array, for example: "cg26257082".
#' @param beta The DNA methylation beta value matrix with rownames annotated to
#' Illumina array probe ID.
#' @param upstream An integer specifying upstream boundary of the promoter region
#' from transcription start site, the maximum is 5000, default for 1500.
#' @param downstream An integer specifying downstream boundary of the promoter region
#' from transcription start site, the maximum is 5000, default for 500.
#' @param method A parameter specifying how the methylation level for a particular
#' miRNA be calculated. "mean" for that average beta value of all CpGs within
#' promoter region, "nearest" for that the beta value of the nearest CpG from
#' the transcription start site as representative, "loess" for that the loess
#' fitted average beta value of all CpGs within promoter region.
#'
#'
#' @return a vector of the methylation level for particular gene across all samples
#' @export
#'
#' @examples
#' data(beta450k.m)
#' cpgPos181a5p.v <- c(-1105, -960, 63, 2362)
#' names(cpgPos181a5p.v) <- c("cg19639967", "cg07837852", "cg26257082", "cg09517106")
#' beta181a5p.v <- getGeneBeta(cpgPos181a5p.v, beta450k.m,  upstream=2000, downstream=2000, method="mean" )
#'
getGeneBeta <- function(pos, beta, upstream=1500, downstream=500, method="mean"){
  if(upstream > 5000){upstream <- 5000}
  if(downstream > 5000){downstream <- 5000}

  upstream <- as.numeric(upstream)  * -1
  downstream <- as.numeric(downstream)

  pos <- pos[intersect(which(pos >= upstream) , which(pos <= downstream) )]
  tmp.idx <- match(names(pos), rownames(beta))
  pos <- pos[which(!is.na(tmp.idx))]
  tmp.idx <- tmp.idx[which(!is.na(tmp.idx))]

  if(length(pos) > 0){
    if(length(pos) == 1){
      geneBeta <- beta[tmp.idx,]
      return(geneBeta)
    }else{
      geneBeta.m <- beta[tmp.idx,]
      if(method == "mean"){
        geneBeta <- colMeans(geneBeta.m)
      }else if(method == "nearest"){
        nearest <- which(abs(pos) == min(abs(pos)))[1]
        geneBeta <- geneBeta.m[nearest,]
      }else if(method == "loess"){
        geneBeta <- apply(geneBeta.m, 2, function(sampleBeta) {
          tmpBeta <- data.frame(x=pos, y=sampleBeta)
          mean(loess(y ~ x, span=10, tmpBeta)$fitted)
        })
      }
      return(geneBeta)
    }
  }
}

