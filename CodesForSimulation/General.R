source("simDat.TRESS.R")
source("simPara.TRESS.R")
source("Utils_simulation.R")
source("Utils_realdata.R")
library(TRESS)

# library(Matrix)
# library(matrixStats)
#library(miceadds)
#source.all("../TRESS_funs/")


#### dataset+strategy: GSE114150+BB, GSE47217+BB, GSE46705+NB, GSE48037+NB, GSE120024+NB

dataset = "GSE114150"

load(paste0("./Results/EstimatedPara_", dataset, ".rda"))



SampleInfo = table(unlist(lapply(strsplit(colnames(allPara$Candidates$Counts), split = "_"  ),
                               function(x) x[1])))
SampleName = names(SampleInfo)

New.reps = c(3,3)
New.design = data.frame(predictor = rep(SampleName, New.reps))
New.model = ~1 + predictor
model.matrix(New.model, New.design)

strategy = SelectSimStrategyByKL(allPara, nreps = New.reps, dataset = dataset,
                                 nsites = 10000, p.prop = 0.1,PLOT = TRUE)

#####

#####

res.para = SimPara.TRESS(nreps = New.reps,
                         nsites = 10000, 
                         p.prop = 0.1,
                         EstiPara = allPara,
                         design = New.design,
                         model = New.model,
                         adjTheta = TRUE,
                         PhiStrategy = "BB")
res.sim = SimDat.TRESS.Coverage(nreps = New.reps,
                                 nsites = length(res.para$flag), 
                                 mu = res.para$mu, 
                                 phi = res.para$phi, 
                                 theta = res.para$theta, 
                                 Coverage = 1,
                                 seed = 12345)
## QQ plot
KL = CompareSimRealCount(Sim.Candidates = res.sim, 
                         Real.Candidates = allPara$Candidates,
                         PLOT = TRUE)


########## compare simulated data with real data
# mean of input and IP
set.seed(seed = 12345)
nsites = min(nrow(allPara$Coef), 10000)
idx = sample(1:nrow(allPara$Coef), nsites)
rcount = sweep(allPara$Candidates$Counts, 2, allPara$Candidates$sf, FUN = "/" )
mean.real.X = rowMeans(rcount[, seq(1, ncol(rcount), 2)],
                       na.rm = TRUE)[idx]
mean.real.Y = rowMeans(rcount[, seq(2, ncol(rcount), 2)], 
                       na.rm = TRUE)[idx]

scount = sweep(res.sim$counts, 2, res.sim$sf, FUN = "/" )
mean.sim.X = rowMeans(scount[, seq(1, ncol(scount), 2)], na.rm = TRUE)
mean.sim.Y = rowMeans(scount[, seq(2, ncol(scount), 2)], na.rm = TRUE)


### scatter plot
par(mfrow = c(1,2))
plot(mean.real.X, mean.sim.X, log = "xy")
abline(0,1, col = "red")

plot(mean.real.Y, mean.sim.Y, log = "xy")
abline(0,1, col = "red")


Xdat = data.frame(real.X = log(mean.real.X+0.5), sim.X = log(mean.sim.X+0.5))
Ydat = data.frame(real.Y = log(mean.real.Y+0.5), sim.Y = log(mean.sim.Y+0.5))

par(mfrow = c(2,2))
boxplot(Xdat[!res.para$flag, ], main = "non-DMRs", ylab = "mean count")
boxplot(Ydat[!res.para$flag, ], main = "non-DMRs", ylab = "mean count")
boxplot(Xdat[res.para$flag, ], main = "DMRs", ylab = "mean count")
boxplot(Ydat[res.para$flag, ], main = "DMRs", ylab = "mean count")




####################################################
library(doParallel)
TRESS.DMR = CallDMRs.paramEsti(counts = res.sim$counts,
                               sf = res.sim$sf,
                               variable = New.design,
                               model = New.model)
TRESS.test = TRESS_DMRtest(DMR = TRESS.DMR, contrast = c(0,1))
pval.TRESS = TRESS.test$pvalue
V = length(which(TRESS.test$padj<0.05&!res.para$flag))
R = length(which(TRESS.test$padj < 0.05))
V/R
power = length(which(TRESS.test$padj<0.05&res.para$flag))/sum(res.para$flag)
c(power,V/R)


#### RADAR
source("DiffPeak_RADAR.R")
radar.idx = c(seq(1, ncol(res.sim$counts), 2), seq(2, ncol(res.sim$counts), 2))
radar.counts = res.sim$counts[, radar.idx]
radar.var <- New.design
res.RADAR = DiffPeak.RADAR(counts = radar.counts, variable = radar.var)
pval.RADAR = res.RADAR$p_value
V = length(which(res.RADAR$fdr<0.05&!res.para$flag))
R = length(which(res.RADAR$fdr < 0.05))
V/R
power = length(which(res.RADAR$fdr<0.05&res.para$flag))/sum(res.para$flag)
c(power,V/R)


### exomePeak2
source("DiffPeak.DESeq2.R")
design.exome2 = data.frame(Reps = c(rep(paste0("Rep", 1:New.reps[1]),each = 2),
                                    rep(paste0("Rep", 1:New.reps[2]),each = 2)),
                           IP = rep(c("Input", "IP"), sum(New.reps)),
                           Trt = rep(c("Ctrl", "T2D"), 2*New.reps ))
model.exome2 = ~ Reps + IP + Trt + IP*Trt
res.exome2 = DiffPeak.DESeq2(counts = res.sim$counts, nreps = New.reps, 
                             model = model.exome2, design = design.exome2) 
pval.exome2 = res.exome2$pval
padj.exome2 = p.adjust(pval.exome2, method = "fdr")
V = length(which( padj.exome2< 0.05&!res.para$flag))
R = length(which(padj.exome2< 0.05))
V/R
power = length(which(padj.exome2<0.05&res.para$flag))/sum(res.para$flag)
c(power,V/R)


#### Statified presentation: obtain strata according to the Input count of DMR
Stratums = Strata_ByDMRInput(res.sim = res.sim, P.adj = TRESS.test$padj)
# save(Stratums, file = paste0("./Results/", dataset, "/Stratums_", dataset, "_",
#                              strategy, ".rda"))


meanInput = rowMeans(res.sim$counts[, seq(1, ncol(res.sim$counts), 2)])


######
OR = 2
rTRESS = rexomePeak2 = rRADAR = matrix(NA, nrow = 7, ncol = ncol(Stratums))
rownames(rTRESS) = rownames(rexomePeak2) = rownames(rRADAR) = 
  c("StratPower", "StratTargetP", "StratTypeI",
    "StratFDR", "StratFDC","StratF1score","StratR10")
for (istra in seq_len(ncol(Stratums))) {
  print(paste0("For regions with mean count in [",Stratums[1, istra], 
               ", ",Stratums[2, istra]  ,") ..."))
  this.id = which(meanInput >= Stratums["Int.L", istra] & meanInput < Stratums["Int.U", istra])
  
  ##### flag density in each strata
  res1 = stra.measure.efsz(flag = res.para$flag[this.id], 
                           padj = TRESS.test$padj[this.id], 
                           fdr0 = 0.05,
                           beta = res.para$beta[this.id], 
                           beta0 = log(OR))
  rTRESS["StratPower", istra] = res1$power
  rTRESS["StratTargetP", istra] = res1$t.power
  rTRESS["StratTypeI", istra] = res1$TypeI 
  rTRESS["StratFDR", istra] = res1$FDR 
  rTRESS["StratFDC", istra] = res1$FDC
  rTRESS["StratF1score", istra] = res1$f1score
  rTRESS["StratR10", istra] = res1$r10
  
  

  res2 = stra.measure.efsz(flag = res.para$flag[this.id], 
                           padj = padj.exome2[this.id], 
                           fdr0 = 0.05,
                           beta = res.para$beta[this.id], 
                           beta0 = log(OR))
  rexomePeak2["StratPower", istra] = res2$power
  rexomePeak2["StratTargetP", istra] = res2$t.power
  rexomePeak2["StratTypeI", istra] = res2$TypeI 
  rexomePeak2["StratFDR", istra] = res2$FDR 
  rexomePeak2["StratFDC", istra] = res2$FDC
  rexomePeak2["StratF1score", istra] = res2$f1score
  rexomePeak2["StratR10", istra] = res2$r10
  
  
  res3 = stra.measure.efsz(flag = res.para$flag[this.id], 
                           padj = res.RADAR$fdr[this.id], 
                           fdr0 = 0.05,
                           beta = res.para$beta[this.id], 
                           beta0 = log(OR))
  rRADAR["StratPower", istra] = res3$power
  rRADAR["StratTargetP", istra] = res3$t.power
  rRADAR["StratTypeI", istra] = res3$TypeI 
  rRADAR["StratFDR", istra] = res3$FDR 
  rRADAR["StratFDC", istra] = res3$FDC
  rRADAR["StratF1score", istra] = res3$f1score
  rRADAR["StratR10", istra] = res3$r10
}
rTRESS
rexomePeak2
rRADAR

matplot(t(rTRESS[1:2, ]), col = c(1,2))

matplot(t(rTRESS[c("StratPower", "StratTargetP"), ]), col = c(1,2))
matplot(t(rexomePeak2[c("StratPower", "StratTargetP"), ]), col = c(1,2))
matplot(t(rRADAR[c("StratPower", "StratTargetP"), ]), col = c(1,2))

matplot(cbind(rTRESS["StratTargetP", ], rexomePeak2["StratTargetP", ],
              rRADAR["StratTargetP", ]))
matplot(cbind(rTRESS["StratFDC", ], rexomePeak2["StratFDC", ],
              rRADAR["StratFDC", ]))



####
method.col = c("#3C5488B2", "#DC0000B2", "gold3")
names(method.col) = c("exomePeak2", "TRESS", "RADAR")
flag = res.para$flag
par(mfrow = c(1,2))
Plot_PrecTop(flag =  matrix(flag, ncol = 1),
             Pvals = list(
               exomePeak2 = matrix(pval.exome2, ncol = 1),
               TRESS = matrix(TRESS.test$pvalue, ncol = 1),
               RADAR = matrix(res.RADAR$fdr, ncol = 1)
             ),
             Ntop = 1000,
             models = c("exomePeak2","TRESS", "RADAR"),
             cols = method.col[c("exomePeak2","TRESS", "RADAR")],
             ltypes = c(4, 1, 2),
             yylim = c(0.2, 1),
             leg = TRUE)
title(main = paste0(New.reps[1], "Reps"))


Plot_ROC(flag = matrix(flag, ncol = 1),
         Pvals = list(
           exomePeak2 = matrix(pval.exome2, ncol = 1),
           TRESS = matrix(TRESS.test$pvalue, ncol = 1),
           RADAR = matrix(res.RADAR$fdr, ncol = 1)
         ),
         models = c("exomePeak2", "TRESS", "RADAR"),
         cols = method.col[c("exomePeak2","TRESS", "RADAR")],
         ltypes = c(4,1, 2),
         leg = TRUE)
title(main = paste0(New.reps[1], " Reps"))


## ROC
method.col = c("#3C5488B2", "#DC0000B2", "gold3", "darkgreen", "blue", "purple")
names(method.col) = c("1.5", "2", "4", "6", "8", "10")
flag = rep(SumInfo$`n = 5`$flag[,1],100)
pdf(file = paste0("./UpdatedFigures/", "ROC_3.pdf"))
Plot_ROC_OR(flag = matrix(flag, ncol = 1),
         Pvals = matrix(unlist(SumInfo$`n = 3`$pval.TRESS), ncol = 1),
         beta = rep(res.para$`n=3`$beta, 100),
         beta0 = c(1.5, 2, 4, 6, 8, 10),
         cols = method.col[c("1.5", "2", "4", "6", "8", "10")],
         leg = TRUE)
title(main = paste0(New.reps[1], " Reps"))
dev.off()
## Dispersion vs. MeanInput
MeanInput = apply(allPara$Candidates$Counts[,seq(1, ncol(allPara$Candidates$Counts), by = 2)], 1, mean)
smoothingSpline = smooth.spline(MeanInput, allPara$Dispersion, spar=0.35)
# plot(x = MeanInput, y = allPara$Dispersion)
# lines(smoothingSpline, col = "red")
pdf(file = paste0("./UpdatedFigures/", "DispersionInput.pdf"))
smoothScatter(allPara$Dispersion ~ MeanInput, xlab = "MeanInput", ylab = "Dispersion", nrpoints = Inf,
              bandwidth = 30, pch = 1, cex = 0.1, main = "GSE114150")
lines(smoothingSpline, col = "red")
dev.off()
boxplot(allPara$Dispersion, main = "GSE114150", ylab = "Dispersion")
library(ggplot2)
dispersion_df = data.frame(Dispersion = allPara$Dispersion)
pdf(file = paste0("./UpdatedFigures/", "Dispersion_violin.pdf"))
ggplot(dispersion_df, aes(x=factor(1), y=Dispersion)) + 
  geom_violin(fill="skyblue", color="darkblue") +  # Fill and border color
  labs(title = "GSE114150", y = "Dispersion", x = "") +
  theme_classic() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())
dev.off()
###############ORs Targeted Power
CUTOFFs = seq(0.05, 0.2, by = 0.05)
OverallMeasure_OR = vector("list", length = length(SumInfo))
names(OverallMeasure_OR) = names(SumInfo)
for (ireps in 1:length(SumInfo)) {
  tmp = SumInfo[[ireps]]
  tmp.Beta = res.para[[ireps]]$beta
  
  rTRESS = vector("list", length = 6)
  names(rTRESS) = c("Power", "TargetP.2","TargetP.4" ,
                                         "TargetP.6","TargetP.8", "TargetP.10")
  rTRESS$Power = rTRESS$TargetP.2 = rTRESS$TargetP.4 = rTRESS$TargetP.6 = rTRESS$TargetP.8 = 
    rTRESS$TargetP.10 = matrix(NA, nrow = ncol(tmp$flag), ncol = length(CUTOFFs))
  
  for (ic in 1:length(CUTOFFs)) {
    cutoff = CUTOFFs[ic]
    for (isim in 1:ncol(SumInfo$`n = 2`$flag)) {
      print(paste0("In the ", isim, "-th run: "), sep = "\n")
      
      ###### TRESS
      res1 = stra.measure.efsz(flag = tmp$flag[, isim], padj = tmp$fdr.TRESS[, isim], 
                               fdr0 = cutoff, beta = tmp.Beta, beta0 = log(2))
      rTRESS$Power[isim, ic] = res1$power
      rTRESS$TargetP.2[isim, ic] = res1$t.power
      res1 = stra.measure.efsz(flag = tmp$flag[, isim], padj = tmp$fdr.TRESS[, isim], 
                               fdr0 = cutoff, beta = tmp.Beta, beta0 = log(4))
      rTRESS$TargetP.4[isim, ic] = res1$t.power 
      res1 = stra.measure.efsz(flag = tmp$flag[, isim], padj = tmp$fdr.TRESS[, isim], 
                               fdr0 = cutoff, beta = tmp.Beta, beta0 = log(6))
      rTRESS$TargetP.6[isim, ic] = res1$t.power 
      res1 = stra.measure.efsz(flag = tmp$flag[, isim], padj = tmp$fdr.TRESS[, isim], 
                               fdr0 = cutoff, beta = tmp.Beta, beta0 = log(8))
      rTRESS$TargetP.8[isim, ic] = res1$t.power
      res1 = stra.measure.efsz(flag = tmp$flag[, isim], padj = tmp$fdr.TRESS[, isim], 
                               fdr0 = cutoff, beta = tmp.Beta, beta0 = log(10))
      rTRESS$TargetP.10[isim, ic] = res1$t.power
    }
    
  }
  OverallMeasure_OR[[ireps]] = list(TRESS = rTRESS)
}
####
for (i in 1:4){
  
  CUTOFFs = seq(0.05, 0.2, 0.05)
  # YLIM = c(1, 1, 1, 1, 1, 1)
  # names(YLIM) = c("Power", "TargetP.2","TargetP.4" ,
  #                 "TargetP.6","TargetP.8", "TargetP.10")
  # 
  # YLIM.min = c(0.4, 0.4, 0.4, 0.4, 0.4, 0.4)
  # names(YLIM.min) = c("Power", "TargetP.2","TargetP.4" ,
  #                     "TargetP.6","TargetP.8", "TargetP.10")
  
  GROUP.COL = c("#3C5488B2", "#DC0000B2", "gold3", "darkgreen", "blue", "purple")
  # names(GROUP.COL$Power) = names(GROUP.COL$TargetP.2)  = 
  #   names(GROUP.COL$TargetP.4)  = names(GROUP.COL$TargetP.6)  = 
  #   names(GROUP.COL$TargetP.8)  = names(GROUP.COL$TargetP.10) = paste0("FDR = ", CUTOFFs)
  
  #group.lty = c(6, 6, 6, 6)
  #names(group.lty) = paste0("FDR = ", CUTOFFs)
  
  m.TRESS = matrix(NA, nrow = length(OverallMeasure_OR), ncol = 6)
  for (ireps in 1:length(OverallMeasure_OR)) {
    m.TRESS[ireps, ] = c(colMeans(OverallMeasure_OR[[ireps]]$TRESS[[1]], na.rm = TRUE)[i],
                         colMeans(OverallMeasure_OR[[ireps]]$TRESS[[2]], na.rm = TRUE)[i],
                         colMeans(OverallMeasure_OR[[ireps]]$TRESS[[3]], na.rm = TRUE)[i],
                         colMeans(OverallMeasure_OR[[ireps]]$TRESS[[4]], na.rm = TRUE)[i],
                         colMeans(OverallMeasure_OR[[ireps]]$TRESS[[5]], na.rm = TRUE)[i],
                         colMeans(OverallMeasure_OR[[ireps]]$TRESS[[6]], na.rm = TRUE)[i])
  }
  colnames(m.TRESS) = c("Power", paste0("OR = ", c(2,4,6,8,10)))
  rownames(m.TRESS) = names(OverallMeasure_OR)
  
  pdf(file = paste0("./UpdatedFigures/", dataset, "_NoStrata_OR", paste0("_FDR", c(0.05, 0.1, 0.15, 0.2))[i],"_TRESS.pdf"))
  par(bty="l",lwd=2,mgp = c(3,1.5,0))
  matplot(m.TRESS,  xaxt = "n", yaxt = "n", 
          xlab = "", 
          ylab = "",
          col = GROUP.COL,
          lty = 6,
          pch = c(2, 20,3,8,12,17),
          type = "o",
          lwd = 3.5,
          bty = 'l', 
          cex.axis = 1.1,
          ylim = c(0.4, 1)
  )
  #title(paste0("TRESS"), line = -0.01, cex.main = 1.3, adj = 0.1)
  
  axis(1, at = 1:nrow(m.TRESS), labels = sub("n = ", "", rownames(m.TRESS)), vjust= 0.3, cex.axis = 2.9, lwd = 2, tck = -0.015) 
  axis(2, cex.axis = 2.9, lwd = 2, tck = -0.015)
  
  dev.off()
  
}

pdf(file = paste0("./UpdatedFigures/", dataset, "_NoStrata_OR_TRESS_legend.pdf"))
plot.new()
legend("bottomright", legend = c("Power", paste0("OR = ", c(2,4,6,8,10))),
       col = GROUP.COL,
       lty = rep(6, 6), pch = c(2, 20,3,8,12,17), bty = "l")
dev.off()

#########FDC over 1
for (i in 1:length(OverallMeasure)){
  fdc_m = OverallMeasure[[i]][["TRESS"]]$FDC
  fdc_perc = apply(fdc_m, 2, function(x) mean(x>1))
  print(paste0(names(OverallMeasure)[i],": ", fdc_perc))
}
#######
method.col = c("#3C5488B2", "#DC0000B2", "gold3", "darkgreen", "blue", "purple")
names(method.col) = c("1.5", "2", "4", "6", "8", "10")
flag = rep(SumInfo$`n = 5`$flag[,1],100)
pdf(file = paste0("./UpdatedFigures/", "PR_3.pdf"))
Plot_PR_OR(flag = matrix(flag, ncol = 1),
            Pvals = matrix(unlist(SumInfo$`n = 3`$pval.TRESS), ncol = 1),
            beta = rep(res.para$`n=3`$beta, 100),
            beta0 = c(1.5, 2, 4, 6, 8, 10),
            cols = method.col[c("1.5", "2", "4", "6", "8", "10")],
            leg = TRUE)
title(main = paste0(New.reps[1], " Reps"))
dev.off()
########Sequencing depth
library(Rsamtools)

bamfs = list.files("~/Library/CloudStorage/OneDrive-Personal/harry/m6a_power/writing/GSE114150")
for (bams in bamfs){
  bamFile <- paste0("~/Library/CloudStorage/OneDrive-Personal/harry/m6a_power/writing/GSE114150/", bamfs[1])
  
  bam <- scanBam(bamFile)
  save(bam, "ki1.rda")
  depth <- Rsamtools::pileup(bamFile)
  
  averageDepth <- mean(depth$pileup)
  print(averageDepth)
}

ll = 1
for (bams in bamfs){
  bamFile <- paste0("~/Library/CloudStorage/OneDrive-Personal/harry/m6a_power/writing/GSE114150/", bams)
  
  bam <- scanBam(bamFile)
  save(bam, paste0("ki", ll, ".rda"))
  depth <- Rsamtools::pileup(bamFile)
  
  averageDepth <- mean(depth$pileup)
  ll = ll + 1
  print(averageDepth)
}

####violin
library(ggplot2)

input_mean = unlist(data.frame(res.sim$counts[,seq(1,11,2)]))
ip_mean = unlist(data.frame(res.sim$counts[,seq(2,12,2)]))
labelled_vector = ifelse(res.para$flag, "DMR", "non-DMR")
data <- data.frame(
  value = c(input_mean, ip_mean),
  half = factor(rep(labelled_vector, 12), levels = c("non-DMR", "DMR")),
  condition = factor(rep(rep(c("Control", "Case"), each = 30000),2), levels = c("Control", "Case")),
  group = factor(rep(c("Input", "IP"), each = 60000))
)
GeomSplitViolin <- ggproto("GeomSplitViolin", GeomViolin, 
                           draw_group = function(self, data, ..., draw_quantiles = NULL) {
                             data <- transform(data, xminv = x - violinwidth * (x - xmin), xmaxv = x + violinwidth * (xmax - x))
                             grp <- data[1, "group"]
                             newdata <- plyr::arrange(transform(data, x = if (grp %% 2 == 1) xminv else xmaxv), if (grp %% 2 == 1) y else -y)
                             newdata <- rbind(newdata[1, ], newdata, newdata[nrow(newdata), ], newdata[1, ])
                             newdata[c(1, nrow(newdata) - 1, nrow(newdata)), "x"] <- round(newdata[1, "x"])
                             
                             if (length(draw_quantiles) > 0 & !scales::zero_range(range(data$y))) {
                               stopifnot(all(draw_quantiles >= 0), all(draw_quantiles <=
                                                                         1))
                               quantiles <- ggplot2:::create_quantile_segment_frame(data, draw_quantiles)
                               aesthetics <- data[rep(1, nrow(quantiles)), setdiff(names(data), c("x", "y")), drop = FALSE]
                               aesthetics$alpha <- rep(1, nrow(quantiles))
                               both <- cbind(quantiles, aesthetics)
                               quantile_grob <- GeomPath$draw_panel(both, ...)
                               ggplot2:::ggname("geom_split_violin", grid::grobTree(GeomPolygon$draw_panel(newdata, ...), quantile_grob))
                             }
                             else {
                               ggplot2:::ggname("geom_split_violin", GeomPolygon$draw_panel(newdata, ...))
                             }
                           })

geom_split_violin <- function(mapping = NULL, data = NULL, stat = "ydensity", position = "identity", ..., 
                              draw_quantiles = NULL, trim = TRUE, scale = "area", na.rm = FALSE, 
                              show.legend = NA, inherit.aes = TRUE) {
  layer(data = data, mapping = mapping, stat = stat, geom = GeomSplitViolin, 
        position = position, show.legend = show.legend, inherit.aes = inherit.aes, 
        params = list(trim = trim, scale = scale, draw_quantiles = draw_quantiles, na.rm = na.rm, ...))
}
# Create the violin plot
pdf(file = paste0("./UpdatedFigures/", "Violin_Input.pdf"))
ggplot(subset(data, group == "Input"), aes(x = half, y = value, fill = condition)) +
  geom_split_violin(trim = T) +
  scale_fill_manual(values = c("Case" = "green4", "Control" = "#3C5488B2")) +
  theme_classic() +
  labs(fill = "condition")+
  ylim(c(0,500))
dev.off()
pdf(file = paste0("./UpdatedFigures/", "Violin_IP.pdf"))
ggplot(subset(data, group == "IP"), aes(x = half, y = value, fill = condition)) +
  geom_split_violin(trim = T) +
  scale_fill_manual(values = c("Case" = "green4", "Control" = "#3C5488B2")) +
  theme_classic() +
  labs(fill = "condition")+
  ylim(c(0,500))
dev.off()
########################
#######Dispersion Plots
for (datas in c("GSE114150", "GSE46705", "GSE48037", "GSE47217", "GSE120024")){
  load(paste0("~/Library/CloudStorage/OneDrive-Personal/harry/m6a_power/Code_forDaoyu/EstimatedPara_", datas, ".rda"))
  MeanInput = apply(allPara$Candidates$Counts[,seq(1, ncol(allPara$Candidates$Counts), by = 2)], 1, mean)
  smoothingSpline = smooth.spline(MeanInput, allPara$Dispersion, spar=0.35)
  # plot(x = MeanInput, y = allPara$Dispersion)
  # lines(smoothingSpline, col = "red")
  pdf(file = paste0("./UpdatedFigures/DispersionInput_", datas, ".pdf"))
  smoothScatter(allPara$Dispersion ~ MeanInput, xaxt = "n", yaxt = "n", xlab = "", ylab = "", nrpoints = Inf,
                bandwidth = 30, pch = 1, cex = 0.1, main = "", cex.axis = 2.9)
  axis(1, padj= 0.6, cex.axis = 2.9, lwd = 2, tck = -0.015) 
  axis(2, cex.axis = 2.9, lwd = 2, tck = -0.015)
  lines(smoothingSpline, col = "red")
  dev.off()
  
  library(ggplot2)
  dispersion_df = data.frame(Dispersion = allPara$Dispersion)
  pdf(file = paste0("./UpdatedFigures/Dispersion_violin_", datas, ".pdf"))
  print(ggplot(dispersion_df, aes(x=factor(1), y=Dispersion)) + 
    geom_violin(fill="skyblue", color="darkblue") +  # Fill and border color
    ylab(" ") +
    xlab(" ") +
    theme_bw() +
    theme_classic() +
    theme(axis.text.x = element_text(size = 33, vjust = 0.85, color = "black"),
          axis.text.y = element_text(hjust = 1, size = 33, color = "black"),
          legend.position = "none",
          axis.line = element_line(colour = "black", size = 0.8), 
          axis.ticks = element_line(colour = "black", size = 0.8),
          axis.ticks.length = unit(0.2, "cm")
    ))
  dev.off()
  rm(allPara)
}

##########distribution of alpha
for (datas in c("GSE114150", "GSE46705", "GSE48037", "GSE47217", "GSE120024")){
  load(paste0("C:/Users/dd284/OneDrive/harry/m6a_power/Code_forDaoyu/EstimatedPara_", datas, ".rda"))
 # MeanInput = apply(allPara$Candidates$Counts[,seq(1, ncol(allPara$Candidates$Counts), by = 2)], 1, mean)
  #smoothingSpline = smooth.spline(MeanInput, allPara$Dispersion, spar=0.35)
  # plot(x = MeanInput, y = allPara$Dispersion)
  # lines(smoothingSpline, col = "red")
 # pdf(file = paste0("./UpdatedFigures/DispersionInput_", datas, ".pdf"))
 # smoothScatter(allPara$Dispersion ~ MeanInput, xaxt = "n", yaxt = "n", xlab = "", ylab = "", nrpoints = Inf,
#              bandwidth = 30, pch = 1, cex = 0.1, main = "", cex.axis = 2.9)
 # axis(1, padj= 0.6, cex.axis = 2.9, lwd = 2, tck = -0.015) 
 # axis(2, cex.axis = 2.9, lwd = 2, tck = -0.015)
 # lines(smoothingSpline, col = "red")
 # dev.off()
  
  library(ggplot2)
  Alpha_df = data.frame(Alpha = allPara$Coef[,1])
  pdf(file = paste0("./UpdatedFigures/Alpha_violin_", datas, ".pdf"))
  print(ggplot(Alpha_df, aes(x=factor(1), y=Alpha)) + 
          geom_violin(fill="skyblue", color="darkblue") +  # Fill and border color
          ylab(expression(alpha)) +
          xlab(" ") +
          theme_bw() +
          theme_classic() +
          theme(axis.text.x = element_text(size = 33, vjust = 0.85, color = "black"),
                axis.text.y = element_text(hjust = 1, size = 33, color = "black"),
                legend.position = "none",
                axis.line = element_line(colour = "black", size = 0.8), 
                axis.ticks = element_line(colour = "black", size = 0.8),
                axis.ticks.length = unit(0.2, "cm"),
                axis.title.y = element_text(size = rel(3.5))
          ) +
          scale_y_continuous(limits = c(-5,10)) +
          geom_hline(yintercept = 0, color = "red", linetype = "dashed")
          
        )
  dev.off()
  rm(allPara)
}




