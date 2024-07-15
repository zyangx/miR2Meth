#' Title Get differential analysis statistics
#'
#' @param data.m data matrix for average DNA methylation probes mapped miRNA
#' promoters or miRNA expression.
#' @param pheno.v The phenotype vector.
#' @param case Case group indicated in the phenotype vector.
#' @param control Control group indicated in the phenotype vector.
#' @param method A parameter specifying which method was used for differential
#' analysis, should be one of the c("limma", "ttest", "wilcox"), default for "limma".
#'
#' @return data matrix for differential analysis statistics.
#' @export
#'
#' @examples
#' data(beta450k.m)
#' miRBeta.m <- probe2miR(beta=beta450k.m, array="450k", type="Mature", upstream=2000, downstream=2000, method="mean")
#' pheno.v  <- c(rep("Normal", 5), rep("Tumor", 5))
#' statM.df <- diffmiR(miRBeta.m, pheno.v, case = "Tumor", control = "Normal", method="limma")
#'
diffmiR <- function(data.m, pheno.v, case="", control="", method="limma"){
  if(method == "limma"){
    stat.df <- Dolimma(data.m, pheno.v, case, control)
  }else if(method == "ttest"){
    stat.df <- Dottest(data.m, pheno.v, case, control)
  }else if(method == "wilcox"){
    stat.df <- DoWilcox(data.m, pheno.v, case, control)
  }
  return(stat.df)
}

#' Title Get differential analysis statistics
#'
#' @param data.m data matrix for average DNA methylation probes mapped miRNA
#' promoters or miRNA expression.
#' @param pheno.v The phenotype vector.
#' @param case Case group indicated in the phenotype vector.
#' @param control Control group indicated in the phenotype vector.
#'
#' @return data matrix for differential analysis statistics.
#' @export
#'
#' @examples
#' data(beta450k.m)
#' miRBeta.m <- probe2miR(beta=beta450k.m, array="450k", type="Mature", upstream=2000, downstream=2000, method="mean")
#' pheno.v  <- c(rep("Normal", 5), rep("Tumor", 5))
#' statM.df <- Dolimma(miRBeta.m, pheno.v, case = "Tumor", control = "Normal")
#'
#' @importFrom limma lmFit contrasts.fit eBayes topTable
Dolimma <- function(data.m, pheno.v, case="", control=""){
  ### construct model matrix
  sampletype.f <- as.factor(pheno.v);
  design.matrix <- model.matrix(~0 + sampletype.f);
  colnames(design.matrix) <- levels(sampletype.f);

  fit  <- lmFit(data.m , design.matrix)
  contrast.matrix <- matrix(0, nrow=ncol(design.matrix) , ncol=1)
  rownames(contrast.matrix) <- colnames(design.matrix);
  colnames(contrast.matrix) <- c(paste(case, "--", control, sep=""));
  contrast.matrix[which(rownames(contrast.matrix) == case), ] <- 1;
  contrast.matrix[which(rownames(contrast.matrix) == control), ] <- -1;
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  DEG.m   <- topTable(fit2, adjust = "fdr", coef = 1,  number=nrow(data.m), p.value = 1)
  return(DEG.m);
}

#' Title Get differential analysis statistics
#'
#' @param data.m data matrix for average DNA methylation probes mapped miRNA
#' promoters or miRNA expression.
#' @param pheno.v The phenotype vector.
#' @param case Case group indicated in the phenotype vector.
#' @param control Control group indicated in the phenotype vector.
#'
#' @return data matrix for differential analysis statistics.
#' @export
#'
#' @examples
#' data(beta450k.m)
#' miRBeta.m <- probe2miR(beta=beta450k.m, array="450k", type="Mature", upstream=2000, downstream=2000, method="mean")
#' pheno.v  <- c(rep("Normal", 5), rep("Tumor", 5))
#' statM.df <- Dottest(miRBeta.m, pheno.v, case = "Tumor", control = "Normal")
#'
Dottest <- function(data.m, pheno.v, case="", control=""){
  lst <- apply(data.m, 1, function(x) t.test(x[which(pheno.v == control)], x[which(pheno.v == case)] ) )
  df <- do.call(rbind, lapply(lst, function(x)  c(logFC = unname(diff(x$estimate)), AveExpr = NA, t = -unname(x$statistic), P.Value = unname(x$p.value) ) ) )
  df[,2] <-  rowMeans(data.m)
  df <- transform(df, adj.P.Val=p.adjust(P.Value, method="BH"))
  DEG.m <- df[order(df[,5]),]
  return(DEG.m);
}

#' Title Get differential analysis statistics
#'
#' @param data.m data matrix for average DNA methylation probes mapped miRNA
#' promoters or miRNA expression.
#' @param pheno.v The phenotype vector.
#' @param case Case group indicated in the phenotype vector.
#' @param control Control group indicated in the phenotype vector.
#'
#' @return data matrix for differential analysis statistics.
#' @export
#'
#' @examples
#' data(beta450k.m)
#' miRBeta.m <- probe2miR(beta=beta450k.m, array="450k", type="Mature", upstream=2000, downstream=2000, method="mean")
#' pheno.v  <- c(rep("Normal", 5), rep("Tumor", 5))
#' statM.df <- DoWilcox(miRBeta.m, pheno.v, case = "Tumor", control = "Normal")
#'
DoWilcox <- function(data.m, pheno.v, case="", control=""){
  DEG.m  <- as.data.frame(matrix(NA, nrow=nrow(data.m), ncol=5) )
  rownames(DEG.m)  <- rownames(data.m)
  colnames(DEG.m)  <- c("logFC", "AveExpr", "t", "P.Value", "adj.P.Val")
  DEG.m[,1] <-  rowMeans(data.m[,which(pheno.v == case)]) - rowMeans(data.m[,which(pheno.v == control)])
  DEG.m[,2] <-  rowMeans(data.m)
  for(i in 1:nrow(data.m)) {
    wilcox.o <- wilcox.test(data.m[i ,which(pheno.v == control)] , data.m[i, which(pheno.v == case)])
    DEG.m[i, 3] <- -wilcox.o$statistic
    DEG.m[i, 4] <- wilcox.o$p.value
  }
  DEG.m[,5] <-  p.adjust(DEG.m[, 4], method="BH")
  DEG.m <- DEG.m[order(DEG.m[,5]),]
  return(DEG.m)
}

