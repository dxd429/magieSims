##################################################  fun for DESeq2
DiffPeak.DESeq2 <- function(counts, nreps, 
                            model,
                            design,
                            sf = NULL) {
  library(DESeq2)
  # design = data.frame(Reps = c(rep(paste0("Rep", 1:nreps[1]),each = 2),
  #                              rep(paste0("Rep", 1:nreps[2]),each = 2)),
  #                     IP = rep(c("Input", "IP"), sum(nreps)), 
  #                     Trt = rep(c("Ctrl", "Trt"), 2*nreps ))
  # model = ~ Reps + IP + Trt + IP*Trt
  # model.matrix(model, design)
  # 
  dds <- DESeqDataSetFromMatrix(countData = counts,
                                colData = design,
                                design = model)
  
  if(length(sf) ==  0){
    ## compute size factor from all counts
    sf = colSums(counts)
    sf = sf/median(sf)
  }
  counts.norm = sweep(counts, 2, sf, FUN = "/")
  sizeFactors(dds) = sf
  dds <- DESeq(dds, test = "Wald")
  
  ###
  #dds <- nbinomWaldTest(dds)
  
  ###
  
  resultsNames(dds) # name of each variable in the design matrix used for DESeq
  res = DESeq2::results(dds, name = "IPIP.TrtT2D") 
  dat = data.frame(pval = res$pvalue, stat = res$stat)
  return(dat)
}