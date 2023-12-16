source("../simDat.TRESS.R")
source("../simPara.TRESS.R")
source("../EstiFromRealDat.R")
source("../Utils_simulation.R")
source("../Utils_realdata.R")
library(Matrix)
library(matrixStats)
library(miceadds)
library(parallel)
source("../DiffPeak.DESeq2.R")
source.all("../ALLFUNS/TRESS_funs/")
source("../DiffPeak_RADAR.R")


#### dataset + strategy: GSE114150+BB, GSE47217+BB, GSE46705+NB, GSE48037+NB, GSE120024+NB


dataset = "GSE114150"
strategy = "NB"

load(paste0("C:/Users/dd284/OneDrive/harry/m6a_power/Code_forDaoyu/EstimatedPara_", dataset, ".rda"))
N.sites = min(10000, nrow(allPara$Coef))

SampleInfo = table(unlist(lapply(strsplit(colnames(allPara$Candidates$Counts), 
                                          split = "_"  ), 
                                 function(x) x[1])))
SampleName = names(SampleInfo)

N.reps.Ctrl = c(2, 3, 5, 7, 10)
N.reps.Trt = c(2, 3, 5, 7, 10)
nsims = 10
SumInfo = vector("list", length = length(N.reps.Ctrl))
for (ireps in seq_along(N.reps.Ctrl)) {
  this.reps = c(N.reps.Ctrl[ireps], N.reps.Trt[ireps])
  New.design = data.frame(predictor = rep(SampleName, this.reps))
  New.model = ~1+predictor
  model.matrix(New.model, New.design)
  
  res.para = SimPara.TRESS(nreps = this.reps,
                           nsites = N.sites, 
                           p.prop = 0.1,
                           EstiPara = allPara,
                           design = New.design,
                           model = New.model,
                           adjTheta = TRUE,
                           PhiStrategy = strategy)
  FLAG = matrix(NA, nrow = N.sites, ncol = nsims)
  pval.TRESS = pval.exomePeak2 =  pval.RADAR = 
    fdr.TRESS = fdr.exomePeak2 = fdr.RADAR = 
    meanInput = matrix(NA, nrow = N.sites, ncol = nsims)
  for (isim in 1:nsims) {
    #####
    res.sim = SimDat.TRESS.Coverage(nreps = this.reps, nsites = N.sites,
                                    mu = res.para$mu, 
                                    phi = res.para$phi, 
                                    theta = res.para$theta, 
                                    Coverage = 1,
                                    seed = 12345+isim)
    meanInput[ ,isim] = rowMeans(res.sim$counts[, seq(1, ncol(res.sim$counts), 2)])
    FLAG[, isim] = res.para$flag
    ####################################################
    TRESS.DMR = CallDMRs.paramEsti(counts = res.sim$counts,
                                   sf = res.sim$sf,
                                   variable = New.design,
                                   model = New.model)
    TRESS.test = TRESS_DMRtest(DMR = TRESS.DMR, contrast = c(0,1))
    pval.TRESS[, isim] = TRESS.test$pvalue
    fdr.TRESS[, isim]  = TRESS.test$padj
    
    ### exomePeak2
    design.exome2 = data.frame(Reps = c(rep(paste0("Rep", 1:this.reps[1]),each = 2),
                                        rep(paste0("Rep", 1:this.reps[2]),each = 2)),
                               IP = rep(c("Input", "IP"), sum(this.reps)),
                               Trt = rep(c("Ctrl", "T2D"), 2*this.reps ))
    model.exome2 = ~ Reps + IP + Trt + IP*Trt
    res.exome2 = DiffPeak.DESeq2(counts = res.sim$counts, nreps = this.reps,
                                 model = model.exome2, design = design.exome2)
    pval.exomePeak2[, isim] = res.exome2$pval
    fdr.exomePeak2[, isim] = p.adjust(res.exome2$pval, method = "fdr")
    
    
    # RADAR
    radar.idx = c(seq(1, ncol(res.sim$counts), 2), seq(2, ncol(res.sim$counts), 2))
    radar.counts = res.sim$counts[, radar.idx]
    radar.var <- New.design
    res.RADAR = DiffPeak.RADAR(counts = radar.counts, variable = radar.var)
    pval.RADAR[, isim] = res.RADAR$p_value
    fdr.RADAR[, isim]  = p.adjust(res.RADAR$p_value, method = "fdr")
  }
  
  # ###
  SumInfo[[ireps]] = list(flag = FLAG,
                          pval.TRESS = pval.TRESS,
                          fdr.TRESS = fdr.TRESS,
                          pval.exomePeak2 = pval.exomePeak2,
                          fdr.exomePeak2 = fdr.exomePeak2,
                          pval.RADAR = pval.RADAR,
                          fdr.RADAR = fdr.RADAR,
                          meanInput = meanInput)
  
}

save(SumInfo, file = "sim_bat1.rda")







