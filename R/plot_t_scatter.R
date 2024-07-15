

#' Title Plot t-value scatter plot from differential promoter DNAm methylation
#' and differential expression.
#'
#' @param statM matrix of differential methylation analysis for the miRNAs
#' @param statE matrix of differential expression analysis for the miRNAs
#' @param FDR_threshold FDR threshold for the significantly differentially methylated
#' miRNAs and differentially expressed miRNAs
#'
#' @return Graphical output of the t-values
#' @export
#'
#' @examples
#' data(beta450k.m)
#' data(miRExpRPM.m)
#' miRBeta.m <- probe2miR(beta=beta450k.m, array="450k", type="Mature", upstream=2000, downstream=2000, method="mean")
#' pheno.v  <- c(rep("Normal", 5), rep("Tumor", 5))
#' statM.df <- DoWilcox(miRBeta.m, pheno.v, case = "Tumor", control = "Normal")
#' logExp.m <- log2(miRExpRPM.m + 1)
#' statE.df <- diffmiR(logExp.m, pheno.v, case = "Tumor", control = "Normal", method="limma")
#' plot_t_scatter(statM.df, statE.df, FDR_threshold = 0.05)
#'
#' @import ggplot2
plot_t_scatter <- function(statM, statE, FDR_threshold = 0.05){
  interMIR.v <- intersect(rownames(statM), rownames(statE))
  statEM <- cbind(statM[interMIR.v, c(3,5)] , statE[interMIR.v, c(3,5)])
  colnames(statEM) <- c("t_M", "FDR_M", "t_E", "FDR_E")
  statEM$Pattern <- "NULL"

  HyperUp.idx <- intersect( intersect(which(statEM[,1] > 0), which(statEM[,2] < FDR_threshold)), intersect( which(statEM[,3] > 0) , which(statEM[,4] < FDR_threshold)) )
  statEM[HyperUp.idx, 5] <- "HyperUp"
  HypoUp.idx <- intersect( intersect(which(statEM[,1] < 0), which(statEM[,2] < FDR_threshold)), intersect( which(statEM[,3] > 0) , which(statEM[,4] < FDR_threshold)) )
  statEM[HypoUp.idx, 5] <- "HypoUp"
  HypoDown.idx <- intersect( intersect(which(statEM[,1] < 0), which(statEM[,2] < FDR_threshold)), intersect( which(statEM[,3] < 0) , which(statEM[,4] < FDR_threshold)) )
  statEM[HypoDown.idx, 5] <- "HypoDown"
  HyperDown.idx <- intersect( intersect(which(statEM[,1] > 0), which(statEM[,2] < FDR_threshold)), intersect( which(statEM[,3] < 0) , which(statEM[,4] < FDR_threshold)) )
  statEM[HyperDown.idx, 5] <- "HyperDown"
  statEM <- within(statEM, Pattern <- factor(Pattern, levels = c("NULL", "HyperUp", "HypoUp", "HypoDown", "HyperDown" )))

  tmp.idx <- which(statEM[,2] < FDR_threshold)
  xmin <- min(abs(statEM[tmp.idx,1]))

  tmp.idx <- which(statEM[,4] < FDR_threshold)
  ymin <- min(abs(statEM[tmp.idx,3]))

  xmax <- max(abs(statEM[,1]))
  ymax <- max(abs(statEM[,3]))

  p <- ggplot(statEM , aes(x=t_M, y=t_E, color=Pattern)) + geom_point() + xlim(-xmax-1, xmax+1) + ylim(-ymax-1, ymax+1) + xlab("t(M)") + ylab("t(E)")
  p <- p + scale_color_manual(values=c("NULL" = "gray", "HyperUp" = "red", "HypoUp" = "blue", "HypoDown" = "orange", "HyperDown" = "yellowgreen" ))
  p <- p + geom_vline(xintercept= c(-xmin, xmin), colour="black", linetype="dashed") + geom_hline(yintercept= c(-ymin, ymin), colour="black", linetype="dashed") +  geom_vline(xintercept=0) +  geom_hline(yintercept=0) + theme_light() + theme(legend.position="none")
  p <- p + annotate(geom="text", x=xmax*0.5,  y=xmax*0.5,  label="Hyper-Up", color="black")
  p <- p + annotate(geom="text", x=-xmax*0.5, y=xmax*0.5,  label="Hypo-Up", color="black")
  p <- p + annotate(geom="text", x=-xmax*0.5, y=-xmax*0.5, label="Hypo-Down", color="black")
  p <- p + annotate(geom="text", x=xmax*0.5,  y=-xmax*0.5, label="Hyper-Down", color="black")
  return(p)

}

