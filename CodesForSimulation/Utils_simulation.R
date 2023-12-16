Plot_ROC <- function(flag, Pvals, models,ltypes, cols = NULL, leg){
  library(ROCR)
  ROC.DE <- function(DE.gs, pval) {
    pred = prediction(1-pval, DE.gs)
    perf = performance(pred,"tpr","fpr")
    perf
  }
  
  AUC.DE <- function(pval, DE.gs) {
    ii = !is.na(pval)
    library(ROCR)
    pred = prediction(1-pval[ii], DE.gs[ii])
    auc.tmp = performance(pred,"auc")
    auc = as.numeric(auc.tmp@y.values)
    auc
  }
  
  nmodel = length(Pvals)
  if(is.null(cols)){
    cols = 1:nmodel
  }
  # all.flag = rep(flag, ncol(Pvals[[1]]))
  
  all.flag = as.vector(flag)
  auc = rep(NA, nmodel)
  for(i in 1:nmodel){
    this.pval = as.vector(Pvals[[i]])
    idx = is.na(this.pval)
    roc = ROC.DE(all.flag[!idx], this.pval[!idx])
    ###
    auc[i] = round(AUC.DE(this.pval, all.flag), 3)
    ###
    par(cex.axis=1.1)
    if (i == 1){
      plot(roc, xlim = c(0,1), lty = ltypes[i], col = cols[i], ylim = c(0,1),
           lwd = 2)
      #axis(2,cex.axis=1.5)
    }else{
      plot(roc, add = TRUE, col = cols[i], lty = ltypes[i],lwd = 2)
      # axis(2,cex.axis=1.5)
      
    }
  }
  if(leg)
    legend("bottomright", legend = paste0(models, ": ", auc), 
           col = cols, lty=ltypes, lwd = 1.5,cex = 0.8, bty = "n")
}


Plot_PrecTop <- function(flag, Pvals, Ntop = 1000,
                         models, cols = NULL, ltypes,
                         yylim,
                         leg,
                         leg.pos = "bottomleft"){
  ### calculate the proportion of true DM sites among the top ranked sites
  Ranks = seq(100, Ntop, 100)
  cum.prop = matrix(NA, nrow = length(Ranks), ncol = length(Pvals))
  for (m in 1:length(Pvals)) {
    thisPval = as.matrix(Pvals[[m]])
    this.cum.prop = matrix(NA, nrow = length(Ranks), ncol = ncol(thisPval))
    for (isim in 1:ncol(thisPval)) {
      pval.sim = thisPval[, isim]
      flag.sim = as.matrix(flag)[, isim]
      tmp.flag = flag.sim[order(pval.sim)]
      
      # tmp.flag = flag[order(pval.sim)] ## commented on May 28, 2022
      
      ### claculate cumulative proportion 
      this.cum.prop[, isim] = cumsum(tmp.flag)[Ranks]/ Ranks
    }
    if(ncol(this.cum.prop) > 1){
      cum.prop[, m] = rowMeans(this.cum.prop, na.rm = TRUE)
    }else{
      cum.prop[, m] = this.cum.prop
    }
  }
  
  
  matplot(cum.prop, type = "l", xaxt = "n", 
          xlab = "Top ranked regions", ylab = "True DM m6A regions",
          col = cols, lty = ltypes, lwd = 2.2,
          cex.axis = 1.1,
          ylim = yylim)
  if(leg){
    legend(leg.pos, legend = models, cex = 0.9, col= cols,  lty = ltypes,
           bty = "n")
  }
  
  
  axis(1, at = seq(length(Ranks)/5, length(Ranks), by = length(Ranks)/5),
       labels = paste0(seq(length(Ranks)/5, length(Ranks), by = length(Ranks)/5), "00"),
       cex.axis = 1.1) 
  
}




### FDR and type I error under threshold t0
FDR_t <- function(t0, pval, fdr = NULL, flag){
  if(length(fdr) == 0)
    fdr = p.adjust(pval, method = "fdr")
  R = sum(fdr < t0, na.rm = TRUE)
  V = sum(fdr < t0 & (!flag), na.rm = TRUE)
  return(V/R)
}

TypeIEr_t <- function(p0, pval, flag){
  V = sum(pval < p0 & (!flag), na.rm = TRUE)
  return(V/sum(!flag))
}

logit <- function(x){
  x[x==0] = 0.0001
  x[x==1] = 0.999
  log(x/(1-x))
}

stra.measure <- function(flag, padj, fdr0){
  
  #### calculate different measures: precision, recall, F1 score, ...
  padj[is.na(padj)] = 1
  num.TP = sum(flag  & padj < fdr0, na.rm = TRUE)
  num.FP = sum(!flag & padj < fdr0, na.rm = TRUE)
  num.TN = sum(!flag & padj >= fdr0, na.rm = TRUE)
  num.FN =  sum(flag & padj >= fdr0, na.rm = TRUE)
  
  # if(sum(flag) ==0 ){
  #   Power = 0
  # }else{
  #   Power = num.TP/(num.TP + num.FN)
  # }
  Power = num.TP/(num.TP + num.FN)
  
  # if(sum(padj < fdr0, na.rm = TRUE) == 0 ){
  #   FDR = 0
  # }else{
  #   FDR = num.FP/(num.TP + num.FP)
  # }
  # Precision = 1 - FDR
  
  Precision = num.TP/(num.TP + num.FP)
  
  F1score = 2*(Power*Precision)/(Power+Precision)
  
  TypeI = num.FP/(num.FP + num.TN)
  res = list(power = Power, precision = Precision, f1score= F1score, typeI = TypeI)
  return(res)
}



stra.measure.efsz <- function(flag, padj, fdr0, beta, beta0){
  
  #### calculate different measures: precision, recall, F1 score, ...
  z = as.numeric(flag)
  z.s = as.numeric(abs(beta)>beta0)
  
  G0 = sum(!flag)
  G1a = sum( abs(beta) > 0 & abs(beta) <= beta0 )
  G1b = sum(abs(beta)>beta0)
  
  ######
  padj[is.na(padj)] = 1
  D = padj < fdr0
  V = sum(!flag & padj < fdr0)
  R = sum(padj < fdr0)
  
  
  #### summarize
  TypeI = V/G0
  FDR = ifelse(R==0, 0, V/R)
  Power = sum( padj < fdr0 & flag, na.rm = TRUE)/sum(flag) 
  Precision = 1-FDR
  F1score = 2*(Power*Precision)/(Power+Precision)
  
  #t.Power = sum(D*z.s)/sum(z.s)
 # FDC = sum(D*(1-z.s))/sum(D*z.s)
  
  t.Power = sum(padj < fdr0 & abs(beta) > beta0)/sum(abs(beta) > beta0)
  Sb = sum(abs(beta) > beta0 & padj < fdr0)
  FDC = V/Sb
  
  r10 = sum(flag)/sum(!flag)
  ####
  res = list(TypeI = TypeI, FDR = FDR, 
             power = Power,
             t.power = t.Power, 
             FDC = FDC,
             f1score = F1score,
             r10 = r10)
  return(res)
}


Plot_ROC_OR <- function(flag, beta, beta0, Pvals,ltypes, cols = NULL, leg){
  library(ROCR)
  ROC.DE <- function(DE.gs, pval) {
    pred = prediction(1-pval, DE.gs)
    perf = performance(pred,"tpr","fpr")
    perf
  }
  
  AUC.DE <- function(pval, DE.gs) {
    ii = !is.na(pval)
    library(ROCR)
    pred = prediction(1-pval[ii], DE.gs[ii])
    auc.tmp = performance(pred,"auc")
    auc = as.numeric(auc.tmp@y.values)
    auc
  }
  
  nOR = length(beta0)
  if(is.null(cols)){
    cols = 1:nOR
  }
  # all.flag = rep(flag, ncol(Pvals[[1]]))
  
  all.flag = as.vector(flag)
  auc = rep(NA, nOR)
  for(i in 1:nOR){
    this.pval = as.vector(Pvals)
    flag_roc = all.flag
    flag_roc[flag_roc & (beta < log(beta0[i]))] = FALSE
    idx = is.na(this.pval)
    roc = ROC.DE(flag_roc[!idx], this.pval[!idx])
    ###
    auc[i] = round(AUC.DE(this.pval, flag_roc), 3)
    ###
    par(xaxt = "n", yaxt = "n")
    if (i == 1){
      plot(roc, ylab = "", xlab = " ", xlim = c(0,1), col = cols[i], ylim = c(0,1),
           lwd = 2)
      par(xaxt = "s", yaxt = "s")
      axis(1,padj= 0.3,
           cex.axis = 2.5, lwd = 2, tck = -0.015)
      axis(2,
           cex.axis = 2.5, lwd = 2, tck = -0.015)
    }else{
      plot(roc, ylab = "", xaxt = "n", yaxt = "n", xlab = " ", add = TRUE, col = cols[i],lwd = 2)
      # axis(2,cex.axis=1.5)
      
    }
  }
  if(leg)
    legend("bottomright", legend = paste0("OR-",beta0, ": ", auc), 
           col = cols, lwd = 2,cex = 1, bty = "n")
}

Plot_PR_OR <- function(flag, beta, beta0, Pvals,ltypes, cols = NULL, leg){
  library(ROCR)
  PR.DE <- function(DE.gs, pval) {
    pred = prediction(1-pval, DE.gs)
    perf = performance(pred,"prec","rec")
    perf
  }
  
  AUC.DE <- function(pval, DE.gs) {
    ii = !is.na(pval)
    library(ROCR)
    pred = prediction(1-pval[ii], DE.gs[ii])
    auc.tmp = performance(pred,"aucpr")
    auc = as.numeric(auc.tmp@y.values)
    auc
  }
  
  nOR = length(beta0)
  if(is.null(cols)){
    cols = 1:nOR
  }
  # all.flag = rep(flag, ncol(Pvals[[1]]))
  
  all.flag = as.vector(flag)
  auc = rep(NA, nOR)
  for(i in 1:nOR){
    this.pval = as.vector(Pvals)
    flag_pr = all.flag
    flag_pr[flag_pr & (beta < log(beta0[i]))] = FALSE
    idx = is.na(this.pval)
    pr = PR.DE(flag_pr[!idx], this.pval[!idx])
    ###
    auc[i] = round(AUC.DE(this.pval, flag_pr), 3)
    ###
    par(xaxt = "n", yaxt = "n")
    if (i == 1){
      plot(pr, ylab = "", xlab = " ", xlim = c(0,1), col = cols[i], ylim = c(0,1),
           lwd = 3, cex = 1.5)
      par(xaxt = "s", yaxt = "s")
      axis(1,padj= 0.3,
           cex.axis = 2.5, lwd = 2, tck = -0.015)
      axis(2,
           cex.axis = 2.5, lwd = 2, tck = -0.015)
    }else{
      plot(pr, ylab = "", xaxt = "n", yaxt = "n", xlab = " ", add = TRUE, col = cols[i],lwd = 3, cex = 1.5)
      #axis(2,cex.axis=1.5)
      #axis(1,cex.axis=1.5)
    }
  }
  if(leg)
    legend("topright", legend = paste0("OR-",beta0, ": ", auc), 
           col = cols, lwd = 2,cex = 1, bty = "n")
}



