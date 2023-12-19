source("./ALLFUNS/simDat.TRESS.R")
source("./ALLFUNS/simPara.TRESS.R")
source("./ALLFUNS/Utils_realdata.R")
library(Matrix)
library(matrixStats)
library(miceadds)
library(ComplexHeatmap)
library(RColorBrewer)
library(jcolors)
library(colorspace)
source.all("./ALLFUNS/TRESS_funs/")

#### dataset+strategy: GSE114150+BB, GSE47217+BB, GSE46705+NB, GSE48037+NB, GSE120024+NB
load("./RawResults/EstimatedPara_GSE114150.rda")
SampleInfo = table(unlist(lapply(strsplit(colnames(allPara$Candidates$Counts), split = "_"  ),
                                 function(x) x[1])))
SampleName = names(SampleInfo)

New.reps = c(2,2)
New.design = data.frame(predictor = rep(SampleName, New.reps))
New.model = ~1 + predictor
model.matrix(New.model, New.design)

res.para = SimPara.TRESS(nreps = New.reps,
                         nsites = 10000, 
                         p.prop = 0.1,
                         EstiPara = allPara,
                         design = New.design,
                         model = New.model,
                         adjTheta = TRUE,
                         PhiStrategy = "BB")
res.sim = SimDat.TRESS(nreps = New.reps,
                       nsites = length(res.para$flag),
                       mu = res.para$mu,
                       phi = res.para$phi,
                       theta = res.para$theta,
                       seed = 12345)

############# qqplot
# mean of input and IP
set.seed(seed = 12345)
nsites = min(nrow(allPara$Coef), 10000)
idx = sample(1:nrow(allPara$Coef), nsites)
##divided by sf
rcount = sweep(allPara$Candidates$Counts, 2, allPara$Candidates$sf, FUN = "/" )
mean.real.X = rowMeans(rcount[, seq(1, ncol(rcount), 2)],
                       na.rm = TRUE)
mean.real.Y = rowMeans(rcount[, seq(2, ncol(rcount), 2)], 
                       na.rm = TRUE)

scount = sweep(res.sim$counts, 2, res.sim$sf, FUN = "/" )
mean.sim.X = rowMeans(scount[, seq(1, ncol(scount), 2)], na.rm = TRUE)
mean.sim.Y = rowMeans(scount[, seq(2, ncol(scount), 2)], na.rm = TRUE)

limits <- c(min(range(mean.real.X)[1], range(mean.sim.X)[1]), max(range(mean.real.X)[2], range(mean.sim.X)[2]))

pdf("./UpdatedFigures/qqplot_meanInput_dividedBysf.pdf")
par(cex.axis = 2.5, lwd = 2, tck = -0.015, mgp = c(3,1.5,0))
qqplot(mean.real.X, mean.sim.X, main = "", xlab = "", ylab = "",
       xlim = limits, ylim = limits)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 3)
title("")
box(lwd=2)
dev.off()


limits <- c(min(range(mean.real.Y)[1], range(mean.sim.Y)[1]), max(range(mean.real.Y)[2], range(mean.sim.Y)[2]))
pdf("./UpdatedFigures/qqplot_meanIP_dividedBysf.pdf")
par(cex.axis = 2.5, lwd = 2, tck = -0.015, mgp = c(3,1.5,0))
qqplot(mean.real.Y, mean.sim.Y, main = "", xlab = "", ylab = "",
       xlim = limits, ylim = limits)
abline(a = 0, b = 1, col = "red", lty = 2, lwd = 3)
title("")
box(lwd = 2)
dev.off()
##Not divided by sf
# rcount = allPara$Candidates$Counts
# mean.real.X = rowMeans(rcount[, seq(1, ncol(rcount), 2)],
#                        na.rm = TRUE)
# mean.real.Y = rowMeans(rcount[, seq(2, ncol(rcount), 2)], c
#                        na.rm = TRUE)
# 
# scount = res.sim$counts
# mean.sim.X = rowMeans(scount[, seq(1, ncol(scount), 2)], na.rm = TRUE)
# mean.sim.Y = rowMeans(scount[, seq(2, ncol(scount), 2)], na.rm = TRUE)
# 
# limits <- c(min(range(mean.real.X)[1], range(mean.sim.X)[1]), max(range(mean.real.X)[2], range(mean.sim.X)[2]))
# 
# pdf("./UpdatedFigures/qqplot_meanInput_NodividedBysf.pdf")
# qqplot(mean.real.X, mean.sim.X, main = "", xlab = "Real Input", ylab = "Simulated Input",
#        xlim = limits, ylim = limits)
# abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
# title("Q-Q Plot")
# dev.off()
# 
# 
# limits <- c(min(range(mean.real.Y)[1], range(mean.sim.Y)[1]), max(range(mean.real.Y)[2], range(mean.sim.Y)[2]))
# pdf("./UpdatedFigures/qqplot_meanIP_NodividedBysf.pdf")
# qqplot(mean.real.Y, mean.sim.Y, main = "", xlab = "Real IP", ylab = "Simulated IP",
#        xlim = limits, ylim = limits)
# abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
# title("Q-Q Plot")
# dev.off()
##divided by sf, not mean
# rcount = sweep(allPara$Candidates$Counts, 2, allPara$Candidates$sf, FUN = "/" )
# real.X = unlist(rcount[, seq(1, ncol(rcount), 2)][idx,])
# real.Y = unlist(rcount[, seq(2, ncol(rcount), 2)][idx,])
# 
# scount = sweep(res.sim$counts, 2, res.sim$sf, FUN = "/" )
# sim.X = unlist(scount[, seq(1, ncol(scount), 2)])
# sim.Y = unlist(scount[, seq(2, ncol(scount), 2)])
# 
# limits <- c(min(range(real.X)[1], range(sim.X)[1]), max(range(real.X)[2], range(sim.X)[2]))
# 
# pdf("./UpdatedFigures/qqplot_Input_dividedBysf.pdf")
# qqplot(real.X, sim.X, main = "", xlab = "Real Input", ylab = "Simulated Input",
#        xlim = limits, ylim = limits)
# abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
# title("Q-Q Plot")
# dev.off()
# 
# 
# limits <- c(min(range(real.Y)[1], range(sim.Y)[1]), max(range(real.Y)[2], range(sim.Y)[2]))
# pdf("./UpdatedFigures/qqplot_IP_dividedBysf.pdf")
# qqplot(real.Y, sim.Y, main = "", xlab = "Real IP", ylab = "Simulated IP",
#        xlim = limits, ylim = limits)
# abline(a = 0, b = 1, col = "red", lty = 2, lwd = 2)
# title("Q-Q Plot")
# dev.off()
####################################sim mu and real m
source("./ALLFUNS/EstiFromRealDat.R")
library(parallel)
dataset = "GSE114150"

load("./ALLFUNS/GSE114150/LivervsKidney_Candidates.rda")
nreps = c(3,3)
colnames(Candidates$Counts) =  c(paste0("Ctrl_input_rep", 1:nreps[1]), paste0("Ctrl_IP_rep", 1:nreps[1]),
                                 paste0("Trt_input_rep", 1:nreps[2]), paste0("Trt_IP_rep", 1:nreps[2]))

allPara = EstiParaFromReal(Candidates = Candidates)
save(allPara, file = "esti_para_mu.rda")
##sim
New.reps = c(3,3)
New.design = data.frame(predictor = rep(SampleName, New.reps))
New.model = ~1 + predictor
mu_sims = NULL
for (i in 1:100){
  res.para = SimPara.TRESS(nreps = New.reps,
                           nsites = 10000,
                           p.prop = 0.1,
                           EstiPara = allPara,
                           design = New.design,
                           model = New.model,
                           adjTheta = TRUE,
                           PhiStrategy = "BB",
                           seed = i)
  mu_sims = cbind(mu_sims, res.para$mu[,c(1,6)])
  idx = res.para$idx
}

mu_real = allPara$mu[,c(1,6)][idx,]

mu_sim_control = mu_sims[,seq(1,ncol(mu_sims), 2)]
mu_sim_case = mu_sims[,seq(2,ncol(mu_sims), 2)]
mu_real_control = mu_real[,1]
mu_real_case = mu_real[,2]


pdf("./UpdatedFigures/mu_scatter.pdf")
par(cex.axis = 2.5, lwd = 2, tck = -0.015, mgp = c(3,1.5,0))
plot(mu_sim_case[324,], main = "", ylab = "", xlab = "")
points(mu_real_case[324],col="blue",pch=19)
text(6, mu_real_case[324],  "Real mu",
     cex=1.5, pos=3,col="blue") 
title("")
box(lwd = 2)
dev.off()
plot(mu_sim_case[324,])
points(mu_real_case[324],col="blue",pch=19)
#####################sequencing depth
rm(list = ls())
load("./ProcessedResults/Coverage_GSE114150.rda" )
Model = "TRESS"
if(Model == "TRESS"){
  Coverages = Coverage.TRESS
} else{
  Coverages = Coverage.exomePeak2
}

colnames(Coverages)[c(3,12,14)]
colnames(Coverages)[c(3,12,14)] = c("Target Power", "Type I error", "FDC")

## reshape data
dat = data.frame(Values = c(Coverages$Power, Coverages$`Target Power`,
                            Coverages$`Type I error`, Coverages$FDC,
                            Coverages$FDR),
                 Replicates = rep(Coverages$Replicates, 5),
                 Depth = rep(Coverages$Depth, 5),
                 Measures = rep(c("Power", "Target Power", "Type I error", "FDC", "FDR"),
                                each = nrow(Coverages) )
)

dat$Replicates = sub("n = ", "", dat$Replicates)
dat$Replicates = factor(dat$Replicates,
                        levels = c("2",  "3",  "5",  "7", "10"))

dat$Depth = factor(dat$Depth,
                   levels = c("0.3",  "0.5",  "0.7",  "1", "3", "5", "7"))
dat$Measures = factor(dat$Measures,
                      levels = c("Power", "Target Power", "Type I error", "FDC", "FDR"))


library(ggplot2)
#Depth.col = c("grey", "skyblue", "yellowgreen", "red", "darksalmon", "wheat2", "gold3" )
Depth.col = rev(viridis::viridis(7))

pdf(file = paste0("./UpdatedFigures/", "TP_box_", Model, ".pdf"))
print(ggplot(dat[dat$Measures == "Target Power", ],
             aes(x= Replicates, y=Values, color = Depth )) +
        geom_boxplot( width = 1, outlier.color = "grey", outlier.size = 0.5)+
        scale_color_manual(values = Depth.col)+
        ylab(" ") +
        xlab(" ") +
        theme_bw()+
        theme_classic() +
        scale_y_continuous(limits = c(0,1)) +
        theme(axis.text.x = element_text(size = 33, vjust = 0.85, color = "black"),
              axis.text.y = element_text(hjust = 1, size = 33, color = "black"),
              legend.position = "none",
              axis.line = element_line(colour = "black", size = 0.8), 
              axis.ticks = element_line(colour = "black", size = 0.8),
              axis.ticks.length = unit(0.2, "cm")
#              legend.margin=margin(-30,1,1,1),
 #             legend.box.margin=margin(-5,-5,-5,-5)
        )
)
dev.off()


pdf(file = paste0("./UpdatedFigures/", "FDC_box_", Model, ".pdf"))
print(ggplot(dat[dat$Measures == "FDC", ],
             aes(x= Replicates, y=Values, color = Depth )) +
        geom_boxplot( width = 1, outlier.color = "grey", outlier.size = 0.5)+
        scale_color_manual(values = Depth.col)+
        ylab(" ") +
        xlab(" ") +
        theme_bw()+
        theme_classic() +
        scale_y_continuous(limits = c(0, ceiling(max(dat[dat$Measures == "FDC", ]$Values)/0.2)*0.2)) +
        theme(axis.text.x = element_text(size = 33, vjust = 0.85, color = "black"),
              axis.text.y = element_text(hjust = 1, size = 33, color = "black"),
              legend.position = "none",
              axis.line = element_line(colour = "black", size = 0.8), 
              axis.ticks = element_line(colour = "black", size = 0.8),
              axis.ticks.length = unit(0.2, "cm")
              #              legend.margin=margin(-30,1,1,1),
              #             legend.box.margin=margin(-5,-5,-5,-5)
        )
)
dev.off()
###########FDR
pdf(file = paste0("./UpdatedFigures/", "FDR_box_", Model, ".pdf"))
print(ggplot(dat[dat$Measures == "FDR", ],
             aes(x= Replicates, y=Values, color = Depth )) +
        geom_boxplot( width = 1, outlier.color = "grey", outlier.size = 0.5)+
        scale_color_manual(values = Depth.col)+
        ylab(" ") +
        xlab(" ") +
        theme_bw()+
        theme_classic() +
        scale_y_continuous(limits = c(0,1)) +
        theme(axis.text.x = element_text(size = 33, vjust = 0.85, color = "black"),
              axis.text.y = element_text(hjust = 1, size = 33, color = "black"),
              legend.position = "none",
              axis.line = element_line(colour = "black", size = 0.8), 
              axis.ticks = element_line(colour = "black", size = 0.8),
              axis.ticks.length = unit(0.2, "cm")
              #              legend.margin=margin(-30,1,1,1),
              #             legend.box.margin=margin(-5,-5,-5,-5)
        )
)
dev.off()


##########legend
library(ggpubr)
pdf(file = paste0("./UpdatedFigures/", "depth_box_legend", Model, ".pdf"))
full_plot = ggplot(dat[dat$Measures == "FDC", ],
             aes(x= Replicates, y=Values, color = Depth)) +
        geom_boxplot( width = 1, outlier.color = "grey", outlier.size = 0.5)+
        scale_color_manual(values = Depth.col)+
        ylab(" ") +
        xlab(" ") +
        theme_bw()+
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45,hjust = 0.8,
                                         size = 20, vjust = 0.85, color = "black"),
              axis.text.y = element_text(hjust = 1, size = 20, color = "black"),
              legend.position = "bottom",
              axis.line = element_line(colour = "black", size = 0.8), 
              axis.ticks = element_line(colour = "black", size = 0.8),
              axis.ticks.length = unit(0.2, "cm"),
              legend.margin=margin(-30,1,1,1),
              legend.box.margin=margin(-5,-5,-5,-5)
        ) +
  guides(color = guide_legend(nrow = 1))

leg = get_legend(full_plot)
print(as_ggplot(leg))
dev.off()

pdf(file = paste0("./UpdatedFigures/depth_box_legend_spectrum", Model, ".pdf"))
legend_image <- as.raster(matrix(Depth.col, ncol=1))
plot(c(0.8,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=1.14, y = seq(0.03,0.97,length.out = 7), labels = rev(c(0.3,0.5,0.7,1,3,5,7)), cex = 3.5)
rasterImage(legend_image, 0, 0, 1,1)
dev.off()

#######TP FDC in one
library(dplyr)

Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(FDC), list(name = mean)) -> FDC
FDC$Replicates = sub("n = ", "", FDC$Replicates)
FDC$Replicates = factor(FDC$Replicates, levels = c("2",  "3",  "5",  "7", "10"))

Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(`Target Power`), list(name = mean)) -> TP
TP$Replicates = sub("n = ", "", TP$Replicates)
TP$Replicates = factor(TP$Replicates, levels = c("2",  "3",  "5",  "7", "10"))

FDC.0.3 = FDC[FDC$Depth == 0.3, ]
TP.0.3 = TP[TP$Depth == 0.3, ]
FDC.0.5 = FDC[FDC$Depth == 0.5, ]
TP.0.5 = TP[TP$Depth == 0.5, ]
FDC.0.7 = FDC[FDC$Depth == 0.7, ]
TP.0.7 = TP[TP$Depth == 0.7, ]
FDC.1 = FDC[FDC$Depth == 1, ]
TP.1 = TP[TP$Depth == 1, ]
FDC.3 = FDC[FDC$Depth == 3, ]
TP.3 = TP[TP$Depth == 3, ]
FDC.5 = FDC[FDC$Depth == 5, ]
TP.5 = TP[TP$Depth == 5, ]
FDC.7 = FDC[FDC$Depth == 7, ]
TP.7 = TP[TP$Depth == 7, ]

pdf(file = paste0("./UpdatedFigures/TPFDC_Line_", Model, ".pdf"))
par(mfrow = c(1,1), mai = c(0.8, 0.8, 0.8, 0.8), mgp = c(3,1.5,0))
plot(NA, xlim=c(1,5), ylim=c(0.05,0.9), ylab = "", xaxt = "n", yaxt = "n", xlab = " ")
axis(1, at = 1:5,
     labels = paste0(c("2",  "3",  "5",  "7", "10")),
     cex.axis = 2.5, lwd = 2, tck = -0.015)

axis(2,
     cex.axis = 2.5, lwd = 2, tck = -0.015)

points(1:5, TP.0.3$name[order(TP.0.3$Replicates)], type="b", pch=17, col = Depth.col[1], lwd=2)
points(1:5, TP.0.5$name[order(TP.0.5$Replicates)], type="b", pch=17, col = Depth.col[2], lwd=2)
points(1:5, TP.0.7$name[order(TP.0.7$Replicates)], type="b", pch=17, col = Depth.col[3], lwd=2)
points(1:5, TP.1$name[order(TP.1$Replicates)], type="b", pch=17, col = Depth.col[4], lwd=2)
points(1:5, TP.3$name[order(TP.3$Replicates)], type="b", pch=17, col = Depth.col[5], lwd=2)
points(1:5, TP.5$name[order(TP.5$Replicates)], type="b", pch=17, col = Depth.col[6], lwd=2)
points(1:5, TP.7$name[order(TP.7$Replicates)], type="b", pch=17, col = Depth.col[7], lwd=2)


par(new = TRUE)                             # Add new plot
plot(1:5, FDC.0.3$name[order(FDC.0.3$Replicates) ], type = "b", pch = 17, lty = "longdash",lwd=2,
     col = Depth.col[1], bg = Depth.col[1], ylab = " ", xlab = " ", ylim = c(0.05,0.9),
     axes = FALSE)
#mtext("FDC", side = 4, line = 2)             # Add second axis label
axis(4, cex.axis = 2.5, lwd = 2, tck = -0.015)


points(1:5, FDC.0.3$name[order(FDC.0.3$Replicates) ], type = "b", pch = 17, lty = "longdash",lwd=2,
       col = Depth.col[1], bg=Depth.col[1])
points(1:5, FDC.0.5$name[order(FDC.0.5$Replicates) ], type = "b", pch = 17,  lty = "longdash",lwd=2,
       col = Depth.col[2], bg=Depth.col[2])
points(1:5, FDC.0.7$name[order(FDC.0.7$Replicates) ], type = "b", pch = 17,  lty = "longdash",lwd=2,
       col = Depth.col[3], bg=Depth.col[3])
points(1:5, FDC.1$name[order(FDC.1$Replicates) ], type = "b", pch = 17,  lty = "longdash",lwd=2,
       col = Depth.col[4], bg=Depth.col[4])
points(1:5, FDC.3$name[order(FDC.3$Replicates) ], type = "b", pch = 17,  lty = "longdash",lwd=2,
       col = Depth.col[5], bg=Depth.col[5])
points(1:5, FDC.5$name[order(FDC.5$Replicates) ], type = "b", pch = 17,  lty = "longdash",lwd=2,
       col = Depth.col[6], bg=Depth.col[6])
points(1:5, FDC.7$name[order(FDC.7$Replicates) ], type = "b", pch = 17,  lty = "longdash",lwd=2,
       col = Depth.col[7], bg=Depth.col[7])
box(lwd = 2)
dev.off()

pdf(file = paste0("./UpdatedFigures/TPFDC_Line_legend_", Model, ".pdf"))
plot.new()
legend("bottomright", legend = paste0("Depth = ", c(0.3, 0.5, 0.7, 1, 3, 5, 7)),
       col = Depth.col,
       lty = rep(2, 7), pch = 17, bty = "n")
dev.off()

pdf(file = paste0("./UpdatedFigures/TPFDC_Line_legend_linetype_", Model, ".pdf"))
plot.new()
legend("bottomright", legend = c("Targeted Power", "FDC"),
       lty = c(1, 2), pch = 17, bty = "n")
dev.off()
#######OR effect
load("./ProcessedResults/Coverage_GSE114150.rda")
#OR_color = c("#F9CF78", "#F39C12", "#F27020", "#E74C3C", "#C0392B", "#9B111B")
#magma_colors = viridis::viridis(6, option = "A")
OR_color = rev(jcolors("pal4"))
if(Model == "TRESS"){
  Coverages = Coverage.TRESS
} else{
  Coverages = Coverage.exomePeak2
}

#colnames(Coverages)[c(3,12,14)]
colnames(Coverages)[c(3,12,14)] = c("Target Power", "Type I error", "FDC")


library(dplyr)
Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(FDC.1.5), list(name = mean)) -> FDC.1.5
FDC.1.5$Replicates = sub("n = ", "", FDC.1.5$Replicates)
FDC.1.5$Replicates = factor(FDC.1.5$Replicates, levels = c("2",  "3",  "5",  "7", "10"))
Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(TargetP.1.5), list(name = mean)) -> TP.1.5
TP.1.5$Replicates = sub("n = ", "", TP.1.5$Replicates)
TP.1.5$Replicates = factor(TP.1.5$Replicates, levels = c("2",  "3",  "5",  "7", "10"))

Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(FDC), list(name = mean)) -> FDC.2
FDC.2$Replicates = sub("n = ", "", FDC.2$Replicates)
FDC.2$Replicates = factor(FDC.2$Replicates, levels = c("2",  "3",  "5",  "7", "10"))
Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(`Target Power`), list(name = mean)) -> TP.2
TP.2$Replicates = sub("n = ", "", TP.2$Replicates)
TP.2$Replicates = factor(TP.2$Replicates, levels = c("2",  "3",  "5",  "7", "10"))


Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(FDC.3), list(name = mean)) -> FDC.3
FDC.3$Replicates = sub("n = ", "", FDC.3$Replicates)
FDC.3$Replicates = factor(FDC.3$Replicates, levels = c("2",  "3",  "5",  "7", "10"))
Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(TargetP.3), list(name = mean)) -> TP.3
TP.3$Replicates = sub("n = ", "", TP.3$Replicates)
TP.3$Replicates = factor(TP.3$Replicates, levels = c("2",  "3",  "5",  "7", "10"))


Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(FDC.4), list(name = mean)) -> FDC.4
FDC.4$Replicates = sub("n = ", "", FDC.4$Replicates)
FDC.4$Replicates = factor(FDC.4$Replicates, levels = c("2",  "3",  "5",  "7", "10"))
Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(TargetP.4), list(name = mean)) -> TP.4
TP.4$Replicates = sub("n = ", "", TP.4$Replicates)
TP.4$Replicates = factor(TP.4$Replicates, levels = c("2",  "3",  "5",  "7", "10"))



Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(FDC.5), list(name = mean)) -> FDC.5
FDC.5$Replicates = sub("n = ", "", FDC.5$Replicates)
FDC.5$Replicates = factor(FDC.5$Replicates, levels = c("2",  "3",  "5",  "7", "10"))
Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(TargetP.5), list(name = mean)) -> TP.5
TP.5$Replicates = sub("n = ", "", TP.5$Replicates)
TP.5$Replicates = factor(TP.5$Replicates, levels = c("2",  "3",  "5",  "7", "10"))



Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(FDC.6), list(name = mean)) -> FDC.6
FDC.6$Replicates = sub("n = ", "", FDC.6$Replicates)
FDC.6$Replicates = factor(FDC.6$Replicates, levels = c("2",  "3",  "5",  "7", "10"))
Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(TargetP.6), list(name = mean)) -> TP.6
TP.6$Replicates = sub("n = ", "", TP.6$Replicates)
TP.6$Replicates = factor(TP.6$Replicates, levels = c("2",  "3",  "5",  "7", "10"))



Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(FDC.7), list(name = mean)) -> FDC.7
FDC.7$Replicates = sub("n = ", "", FDC.7$Replicates)
FDC.7$Replicates = factor(FDC.7$Replicates, levels = c("2",  "3",  "5",  "7", "10"))
Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(TargetP.7), list(name = mean)) -> TP.7
TP.7$Replicates = sub("n = ", "", TP.7$Replicates)
TP.7$Replicates = factor(TP.7$Replicates, levels = c("2",  "3",  "5",  "7", "10"))



Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(FDC.8), list(name = mean)) -> FDC.8
FDC.8$Replicates = sub("n = ", "", FDC.8$Replicates)
FDC.8$Replicates = factor(FDC.8$Replicates, levels = c("2",  "3",  "5",  "7", "10"))
Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(TargetP.8), list(name = mean)) -> TP.8
TP.8$Replicates = sub("n = ", "", TP.8$Replicates)
TP.8$Replicates = factor(TP.8$Replicates, levels = c("2",  "3",  "5",  "7", "10"))



Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(FDC.9), list(name = mean)) -> FDC.9
FDC.9$Replicates = sub("n = ", "", FDC.9$Replicates)
FDC.9$Replicates = factor(FDC.9$Replicates, levels = c("2",  "3",  "5",  "7", "10"))
Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(TargetP.9), list(name = mean)) -> TP.9
TP.9$Replicates = sub("n = ", "", TP.9$Replicates)
TP.9$Replicates = factor(TP.9$Replicates, levels = c("2",  "3",  "5",  "7", "10"))


Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(FDC.10), list(name = mean)) -> FDC.10
FDC.10$Replicates = sub("n = ", "", FDC.10$Replicates)
FDC.10$Replicates = factor(FDC.10$Replicates, levels = c("2",  "3",  "5",  "7", "10"))
Coverages %>%
  group_by(Replicates, Depth) %>%
  summarise_at(vars(TargetP.10), list(name = mean)) -> TP.10
TP.10$Replicates = sub("n = ", "", TP.10$Replicates)
TP.10$Replicates = factor(TP.10$Replicates, levels = c("2",  "3",  "5",  "7", "10"))



#############################  Effect of OR on FDC, and target power when depth = 1
d = 1
pdf(file = paste0("./UpdatedFigures/EffectofOR_FDC_", Model, ".pdf"))
par(mgp = c(3,1.5,0))
plot(NA, xlim=c(1,5), ylim=c(0,10), ylab = "", xaxt = "n", yaxt = "n", xlab = " ")
axis(1, at = 1:5,
     labels = c("2",  "3",  "5",  "7", "10"),
     cex.axis = 2.5, lwd = 2, tck = -0.015)
axis(2,
     cex.axis = 2.5, lwd = 2, tck = -0.015)
tmp = FDC.1.5[FDC.1.5$Depth == d, ]
points(1:5, tmp$name[order(tmp$Replicates)], type="b", pch=20, col = OR_color[1], lwd=2)

tmp = FDC.2[FDC.2$Depth == d, ]
points(1:5, tmp$name[order(tmp$Replicates)], type="b", pch=20, col = OR_color[2], lwd=2)

tmp = FDC.4[FDC.3$Depth == d, ]
points(1:5, tmp$name[order(tmp$Replicates)], type="b", pch=20, col = OR_color[3], lwd=2)

tmp = FDC.6[FDC.3$Depth == d, ]
points(1:5, tmp$name[order(tmp$Replicates)], type="b", pch=20, col = OR_color[4], lwd=2)

tmp = FDC.8[FDC.3$Depth == d, ]
points(1:5, tmp$name[order(tmp$Replicates)], type="b", pch=20, col = OR_color[5], lwd=2)

tmp = FDC.10[FDC.10$Depth == d, ]
points(1:5, tmp$name[order(tmp$Replicates)], type="b", pch=20, col = OR_color[6], lwd=2)

abline(h = 1, lty="dashed", col = "grey")


box(lwd = 2)
#title(paste0(""), line = 0.5)
dev.off()


pdf(file = paste0("./UpdatedFigures/EffectofOR_legend_", Model, ".pdf"))
legend_image <- as.raster(matrix(OR_color, ncol=1))
plot(c(0.8,2),c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main = '')
text(x=1.14, y = seq(0.03,0.97,length.out = 6), labels = rev(c(1.5,2,4,6,8,10)), cex = 3.5)

rasterImage(legend_image, 0, 0, 1,1)
dev.off()

#
d = 1
pdf(file = paste0("./UpdatedFigures/EffectofOR_TP_", Model, ".pdf"))
par(mgp = c(3,1.5,0))
plot(NA, xlim=c(1,5), ylim=c(0,1), ylab = "", xaxt = "n", yaxt = "n", xlab = " ")
axis(1, at = 1:5,
     labels = c("2",  "3",  "5",  "7", "10"),
     cex.axis = 2.5, lwd = 2, tck = -0.015)
axis(2,
     cex.axis = 2.5, lwd = 2, tck = -0.015)
tmp = TP.1.5[TP.1.5$Depth == d, ]
points(1:5, tmp$name[order(tmp$Replicates)], type="b", pch=20, col = OR_color[1], lwd=2)

tmp = TP.2[TP.2$Depth == d, ]
points(1:5, tmp$name[order(tmp$Replicates)], type="b", pch=20, col = OR_color[2], lwd=2)

tmp = TP.4[TP.4$Depth == d, ]
points(1:5, tmp$name[order(tmp$Replicates)], type="b", pch=20, col = OR_color[3], lwd=2)

tmp = TP.6[TP.6$Depth == d, ]
points(1:5, tmp$name[order(tmp$Replicates)], type="b", pch=20, col = OR_color[4], lwd=2)

tmp = TP.8[TP.8$Depth == d, ]
points(1:5, tmp$name[order(tmp$Replicates)], type="b", pch=20, col = OR_color[5], lwd=2)

tmp = TP.10[TP.10$Depth == d, ]
points(1:5, tmp$name[order(tmp$Replicates)], type="b", pch=20, col = OR_color[6], lwd=2)


box(lwd = 2)
dev.off()
##############by FDR threshold
rm(list = ls())
dataset = "GSE114150"
strategy = "BB"
leng.pos = c("bottomright", "bottomright", "topright", "topright", "topright" )
names(leng.pos) = c("Power", "TargetP", "TypeI", "FDC", "FDR")
for (Measures in c("Power", "TargetP", "TypeI", "FDC", "FDR")){
  
  ######## load un-stratified criteria
  ### each criteria matrix has 4 columns, corresponding to 4 different 
  ###    FDR thresholds: [0.05, 0.1, 0.15, 0.2]
  load("./ProcessedResults/OverallMeasure_GSE114150_BB.rda") 
  
  
  CUTOFFs = seq(0.05, 0.2, 0.05)
  YLIM = c(1, 1, 0.1, 1.5, 1, 1, 0.2)
  names(YLIM) = c("Power",  "TargetP", "TypeI", "FDC", "FDR", "F1score", "R10")
  
  YLIM.min = c(0.4, 0.4, 0, 0, 0, 0.2, 0)
  names(YLIM.min) = c("Power",  "TargetP", "TypeI", "FDC", "FDR", "F1score", "R10")
  
  GROUP.COL = list(Power = c("lightpink", "lightpink1", "lightpink2", "lightpink3"),
                   TargetP = c("lightpink", "lightpink1", "lightpink2", "lightpink3"),
                   FDC = c("cadetblue2", "cadetblue3", "cadetblue", "cadetblue4"), 
                   TypeI = c("cadetblue2", "cadetblue3", "cadetblue", "cadetblue4"),
                   FDR = c("cadetblue2", "cadetblue3", "cadetblue", "cadetblue4")
  )
  names(GROUP.COL$Power) = names(GROUP.COL$TargetP)  = 
    names(GROUP.COL$TypeI)  = names(GROUP.COL$FDR)  = 
    names(GROUP.COL$FDC)  = paste0("FDR = ", CUTOFFs)
  
  group.lty = c(6, 6, 6, 6)
  names(group.lty) = paste0("FDR = ", CUTOFFs)
  
  m.TRESS = m.exomePeak2 = matrix(NA, nrow = length(OverallMeasure), ncol = length(CUTOFFs))
  for (ireps in 1:length(OverallMeasure)) {
    id = which(grepl( Measures, names(OverallMeasure[[ireps]]$TRESS) ))
    m.TRESS[ireps, ] = colMeans(OverallMeasure[[ireps]]$TRESS[[id]], na.rm = TRUE)
    
    id = which(grepl( Measures, names(OverallMeasure[[ireps]]$exomePeak2) ))
    m.exomePeak2[ireps, ] = colMeans(OverallMeasure[[ireps]]$exomePeak2[[id]], na.rm = TRUE)
  }
  colnames(m.TRESS) = colnames(m.exomePeak2) = paste0("FDR = ", seq(0.05, 0.2, 0.05))
  rownames(m.TRESS) = rownames(m.exomePeak2) = names(OverallMeasure)
  
  pdf(file = paste0("./UpdatedFigures/", dataset, "_NoStrata_", Measures,"_TRESS.pdf"))
  par(bty="l",lwd=2,mgp = c(3,1.5,0))
  matplot(m.TRESS,  xaxt = "n", yaxt = "n", 
          xlab = "", 
          ylab = "",
          col = GROUP.COL[[Measures]],
          lty = group.lty[colnames(m.TRESS)],
          pch = c(2, 20,3,8),
          type = "o",
          lwd = 3.5,
          bty = 'l', 
          cex.axis = 1.1,
          ylim = c(YLIM.min[Measures], YLIM[Measures])
  )
  #title(paste0("TRESS"), line = -0.01, cex.main = 1.3, adj = 0.1)
  
  axis(1, at = 1:nrow(m.TRESS), labels = sub("n = ", "", rownames(m.TRESS)), vjust= 0.3, cex.axis = 2.9, lwd = 2, tck = -0.015) 
  axis(2, cex.axis = 2.9, lwd = 2, tck = -0.015)
  
  dev.off()
  
  pdf(file = paste0("./UpdatedFigures/", dataset, "_NoStrata_", Measures,"_exomePeak2.pdf"))
  par(bty="l",lwd=2, mgp = c(3,1.5,0))
  matplot(m.exomePeak2,  xaxt = "n", yaxt = "n", 
          xlab = "", 
          ylab = "",
          col = GROUP.COL[[Measures]],
          lty = group.lty[colnames(m.exomePeak2)],
          lwd = 3.5,
          pch = c( 2,20,3,8),
          type = "o",
          bty = 'l', 
          cex.axis = 1.1,
          ylim = c(YLIM.min[Measures], YLIM[Measures])
  )
  #title(paste0("exomePeak2"), line = -0.01, cex.main = 1.3, adj = 0.1)

  axis(1, at = 1:nrow(m.exomePeak2), labels = sub("n = ", "", rownames(m.exomePeak2)), vjust= 0.3, cex.axis = 2.9, lwd = 2, tck = -0.015) 
  axis(2, cex.axis = 2.9, lwd = 2, tck = -0.015)
  dev.off()
  
  pdf(file = paste0("./UpdatedFigures/", dataset, "_NoStrata_", Measures,"_legend.pdf"))
  plot.new()
  legend(leng.pos[Measures], legend = colnames(m.exomePeak2), 
         cex = 0.8, 
         col = GROUP.COL[[Measures]],
         bty = "n", 
         pch = c( 2, 20,3,8),
         lty = group.lty[colnames(m.exomePeak2)],
         lwd = rep(2, 1))
  dev.off()
}
########strata
rm(list = ls())
dataset = "GSE114150"
strategy = "BB"
### stratified value for each criteria: 
### each criteria matrix has six columns, corresponding to six stratum:
###  (0, 10%] (10%, 30%] (30%, 50%] (50%, 70%] (70%, 90%] (90%, 100%]
### The exact count range for each stratum is saved in Stratums_GSE114150_BB.rda
load("./ProcessedResults/StratMeasure_GSE114150_BB.rda") 
load("./RawResults/Stratums_GSE114150_BB.rda") 


Stratums = round(Stratums)
CUTOFFs = seq(0.05, 0.2, 0.05)
YLIM = c(1, 1, 0.05, 1.0, 0.7, 1, 0.2)
names(YLIM) = c("Power",  "TargetP", "TypeI", "FDC", "FDR", "F1score", "R10")
group.col = rev(colorspace::heat_hcl(5, c = c(80,30), l = c(30,90), power = c(1/5,1.5)))
group.col[1] = "#E6D895"
group.lty = c( 3, 6, 4, 2, 5)
names(group.col) = c("n = 2", "n = 3", "n = 5", "n = 7", "n = 10")
names(group.lty) = c("n = 2", "n = 3", "n = 5", "n = 7", "n = 10")
for (Measures in c("Power", "TargetP", "TypeI", "FDC", "FDR")){
  
  leng.pos = c("bottomright", "bottomright", "topleft", "topleft", "topleft" )
  names(leng.pos) = c("Power", "TargetP", "TypeI", "FDC", "FDR")
  
  txt.dist = -YLIM/5 
  
  for (ic in 1:length(CUTOFFs)) {
    ########
    m.TRESS = m.exomePeak2 = matrix(NA, nrow = ncol(Stratums), ncol = length(StratMeasure))
    for (ireps in 1:length(StratMeasure)) {
      id = which(grepl( Measures, names(StratMeasure[[ireps]][[ic]]$TRESS) ))
      m.TRESS[, ireps] = colMeans(StratMeasure[[ireps]][[ic]]$TRESS[[id]], na.rm = TRUE)
      
      id = which(grepl( Measures, names(StratMeasure[[ireps]][[ic]]$exomePeak2) ))
      m.exomePeak2[, ireps] = colMeans(StratMeasure[[ireps]][[ic]]$exomePeak2[[id]], na.rm = TRUE)
    }
    colnames(m.TRESS) = colnames(m.exomePeak2) = names(StratMeasure)
    rownames(m.TRESS) = rownames(m.exomePeak2) = c(paste0("(", Stratums[,1][1], ", ", Stratums[,1][2], "]"),
                                                   paste0("(", Stratums[,2][1], ", ", Stratums[,2][2], "]"),
                                                   paste0("(", Stratums[,3][1], ", ", Stratums[,3][2], "]"),
                                                   paste0("(", Stratums[,4][1], ", ", Stratums[,4][2], "]"),
                                                   paste0("(", Stratums[,5][1], ", ", Stratums[,5][2], "]"),
                                                   paste0("(", Stratums[,6][1], ", ", "Inf", "]"))
    
    
    ####
    pdf(file = paste0("./UpdatedFigures/", dataset, "_ForEachCutoff_", Measures,"_TRESS_",
                      CUTOFFs[ic],
                      ".pdf"))
    par(bty="l",lwd=2, mai = c(1.9, 0.82, 0.6, 0.42))
    ### TRESS
    matplot(m.TRESS, 
            xaxt = "n", 
            yaxt = "n",
            xlab = " ", # "Stratum", 
            ylab = "",
            col = group.col[colnames(m.TRESS)], 
            lty = group.lty[colnames(m.TRESS)],
            lwd = 3.5,
            pch = 0:4,
            type = "o",
            bty = "l",
            cex.axis = 1.1,
            ylim = c(0, YLIM[Measures])
    )
    # title(paste0("TRESS: ", names(StratMeasure$`n = 2`)[ic]), 
    #       line = -0.01, cex.main = 1.3, adj = 0.1)

    axis(1, at = 1:ncol(Stratums), labels = FALSE, 
         cex.axis = 2.5,
        lwd = 2, tck = -0.015) 
    axis(2,
         cex.axis = 2.5, lwd = 2, tck = -0.015)
    text(x = 1:ncol(Stratums),
         y = txt.dist[Measures],
         labels = c(paste0("(", Stratums[,1][1], ", ", Stratums[,1][2], "]"),
                    paste0("(", Stratums[,2][1], ", ", Stratums[,2][2], "]"),
                    paste0("(", Stratums[,3][1], ", ", Stratums[,3][2], "]"),
                    paste0("(", Stratums[,4][1], ", ", Stratums[,4][2], "]"),
                    paste0("(", Stratums[,5][1], ", ", Stratums[,5][2], "]"),
                    paste0("(", Stratums[,6][1], ", ", "Inf", "]")),
         xpd = NA,
         srt = 45, ### rotate by 45
         cex = 2.5,
         adj = 0.7)
    
     dev.off()
     

     if(Measures == "FDC"){
       heat_col = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "Greens")))(100))
     }else if (Measures == "FDR"){
       heat_col = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "Blues")))(100))
     }else {
       heat_col = rev(colorRampPalette(rev(brewer.pal(n = 7, name = "Reds")))(100))
     }
     
     pdf(paste0("./UpdatedFigures/", dataset, "_HeatmapStrata_", Measures,"_TRESS_",
                CUTOFFs[ic],
                ".pdf"))
     draw(Heatmap(t(m.TRESS[,rev(colnames(m.TRESS))]),
             cluster_rows = F,
             cluster_columns = F,
             # width = ncol(mean_cor1)*unit(25, "mm"), 
             # height = nrow(mean_cor1)*unit(25, "mm"),
             rect_gp = gpar(col = "grey", lwd = 2, type = "none"),
             cell_fun = function(j, i, x, y, width, height, fill) {
               
               
               grid.rect(x, y, width, height, gp = gpar(fill = fill, col = fill))
               grid.text(sprintf("%.2f", t(m.TRESS[,rev(colnames(m.TRESS))])[i, j]), x, y, gp = gpar(fontsize = 22))
               
             },
             column_names_rot = 45,
             column_names_centered = T,
             
             col = heat_col,
             heatmap_legend_param = list(
               legend_gp = gpar(fontsize = 22),
               title = ""
             )
             ,row_names_gp = gpar(fontsize = 25),
             column_names_gp = gpar(fontsize = 25)
             ,show_heatmap_legend = T
     ))
    dev.off()
    ####
    pdf(file = paste0("./UpdatedFigures/", dataset, "_ForEachCutoff_", Measures,"_exomePeak2_",
                      CUTOFFs[ic],
                      ".pdf"))
    par(bty="l",lwd=2, mai = c(1.9, 0.82, 0.6, 0.42))
    matplot(m.exomePeak2, xaxt = "n", yaxt = "n",
            xlab = " ", # "Stratum", 
            ylab = "",
            col = group.col[colnames(m.TRESS)], 
            lty = group.lty[colnames(m.TRESS)],
            lwd = 3.5,
            pch = 0:4,
            bty = "l",
            type = "o",
            cex.axis = 1.1,
            ylim = c(0, YLIM[Measures])
    )
    # title(paste0("exomePeak2: ", names(StratMeasure$`n = 2`)[ic]), 
    #       line = -0.01, cex.main = 1.3, adj = 0.1)

    axis(1, at = 1:ncol(Stratums), labels = FALSE,
         cex.axis = 2.5,
         lwd = 2, tck = -0.015) 
    axis(2,
         cex.axis = 2.5, lwd = 2, tck = -0.015)
    text(x = 1:ncol(Stratums),
         y = txt.dist[Measures], 
         labels = c(paste0("(", Stratums[,1][1], ", ", Stratums[,1][2], "]"),
                    paste0("(", Stratums[,2][1], ", ", Stratums[,2][2], "]"),
                    paste0("(", Stratums[,3][1], ", ", Stratums[,3][2], "]"),
                    paste0("(", Stratums[,4][1], ", ", Stratums[,4][2], "]"),
                    paste0("(", Stratums[,5][1], ", ", Stratums[,5][2], "]"),
                    paste0("(", Stratums[,6][1], ", ", "Inf", "]")),
         xpd = NA,
         srt = 45, 
         cex = 2.5,
         adj = 0.7)
     dev.off()
     
     pdf(paste0("./UpdatedFigures/", dataset, "_HeatmapStrata_", Measures,"_exomePeak2_",
                CUTOFFs[ic],
                ".pdf"))
     draw(Heatmap(t(m.exomePeak2[,rev(colnames(m.exomePeak2))]),
             cluster_rows = F,
             cluster_columns = F,
             # width = ncol(mean_cor1)*unit(25, "mm"), 
             # height = nrow(mean_cor1)*unit(25, "mm"),
             #rect_gp = gpar(col = "grey", lwd = 2, type = "none"),
             cell_fun = function(j, i, x, y, width, height, fill) {
               
               
               grid.rect(x, y, width, height, gp = gpar(fill = fill, col = fill))
               grid.text(sprintf("%.2f", t(m.exomePeak2[,rev(colnames(m.exomePeak2))])[i, j]), x, y, gp = gpar(fontsize = 22))
               
             },
             column_names_rot = 45,
             column_names_centered = T,
             
             col = heat_col,
             heatmap_legend_param = list(
               legend_gp = gpar(fontsize = 22),
               title = ""
             )
             ,row_names_gp = gpar(fontsize = 25),
             column_names_gp = gpar(fontsize = 25)
             ,show_heatmap_legend = T
     ))
     dev.off()
     
     pdf(file = paste0("./UpdatedFigures/", dataset, "_ForEachCutoff_legend",
                       ".pdf"))
     plot.new()
     legend(leng.pos[Measures], legend = colnames(m.TRESS), 
            cex = 0.8, 
            col = group.col[colnames(m.TRESS)], 
            bty = "n", 
            pch = 0:4,
            lty = group.lty[colnames(m.TRESS)],
            lwd = rep(2, 1))
     dev.off()
  }
}
######heatmap for strata
# library(ComplexHeatmap)
# pdf("heatmap.pdf")
# Heatmap(t(m.TRESS[,rev(colnames(m.TRESS))]),
#         cluster_rows = F,
#         cluster_columns = F,
#         # width = ncol(mean_cor1)*unit(25, "mm"), 
#         # height = nrow(mean_cor1)*unit(25, "mm"),
#         rect_gp = gpar(col = "grey", lwd = 2, type = "none"),
#         cell_fun = function(j, i, x, y, width, height, fill) {
#           
#         
#             grid.rect(x, y, width, height, gp = gpar(fill = fill, col = fill))
#             grid.text(sprintf("%.2f", t(m.TRESS[,rev(colnames(m.TRESS))])[i, j]), x, y, gp = gpar(fontsize = 22))
#          
#         },
#         column_names_rot = 315,
#         column_names_centered = T,
#         
#         col = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
#         heatmap_legend_param = list(
#           legend_gp = gpar(fontsize = 12),
#           title = ""
#         )
#         ,row_names_gp = gpar(fontsize = 12),
#         column_names_gp = gpar(fontsize = 12)
#         ,show_heatmap_legend = T
# )
# dev.off()
######################violin plot of stratified results
strata.col = rev(viridis::viridis(6))
Model = "TRESS"

load("./ProcessedResults/StratMeasure_GSE114150_BB.rda") 
load("./RawResults/Stratums_GSE114150_BB.rda") 

Stratums = round(Stratums)

strata_names = c(paste0("(", Stratums[,1][1], ", ", Stratums[,1][2], "]"),
                 paste0("(", Stratums[,2][1], ", ", Stratums[,2][2], "]"),
                 paste0("(", Stratums[,3][1], ", ", Stratums[,3][2], "]"),
                 paste0("(", Stratums[,4][1], ", ", Stratums[,4][2], "]"),
                 paste0("(", Stratums[,5][1], ", ", Stratums[,5][2], "]"),
                 paste0("(", Stratums[,6][1], ", ", "Inf", "]"))
library(reshape2)
tempsave = StratMeasure$`n = 5`$`FDR = 0.05`$TRESS$StratTargetP
colnames(tempsave) = strata_names
tp_5_0.05 = melt(tempsave)
names(tp_5_0.05)[2:3] = c("Strata", "TP")
tp_5_0.05 = tp_5_0.05[,-1]

pdf(file = paste0("./UpdatedFigures/", "TP_violin_", Model, ".pdf"))
print(ggplot(tp_5_0.05,
             aes(x= Strata, y = TP, color = Strata )) +
        geom_violin(trim=T, size = 1, width = 1) +
        geom_point(aes(fill = Strata, color = Strata), position = position_jitter(seed = 1, width = 0.2))+
        scale_color_manual(values = strata.col)+
        ylab(" ") +
        xlab(" ") +
        theme_bw()+
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45,hjust = 0.8,
                                         size = 29, vjust = 0.85, color = "black"),
              axis.text.y = element_text(hjust = 1, size = 33, color = "black"),
              legend.position = "none",
              axis.line = element_line(colour = "black", size = 0.8), 
              axis.ticks = element_line(colour = "black", size = 0.8),
              axis.ticks.length = unit(0.2, "cm")
              #              legend.margin=margin(-30,1,1,1),
              #             legend.box.margin=margin(-5,-5,-5,-5)
        )
)
dev.off()

tempsave = StratMeasure$`n = 5`$`FDR = 0.05`$TRESS$StratFDC
colnames(tempsave) = strata_names
fdc_5_0.05 = melt(tempsave)
names(fdc_5_0.05)[2:3] = c("Strata", "FDC")
fdc_5_0.05 = fdc_5_0.05[,-1]

pdf(file = paste0("./UpdatedFigures/", "FDC_violin_", Model, ".pdf"))
print(ggplot(fdc_5_0.05,
             aes(x= Strata, y = FDC, color = Strata )) +
        geom_violin(trim=T, size = 1, width = 1) +
        geom_point(aes(fill = Strata, color = Strata), position = position_jitter(seed = 1, width = 0.2))+
        scale_color_manual(values = strata.col)+
        ylab(" ") +
        xlab(" ") +
        theme_bw()+
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45,hjust = 0.8,
                                         size = 29, vjust = 0.85, color = "black"),
              axis.text.y = element_text(hjust = 1, size = 33, color = "black"),
              legend.position = "none",
              axis.line = element_line(colour = "black", size = 0.8), 
              axis.ticks = element_line(colour = "black", size = 0.8),
              axis.ticks.length = unit(0.2, "cm")
              #              legend.margin=margin(-30,1,1,1),
              #             legend.box.margin=margin(-5,-5,-5,-5)
        ) + 
        ylim(c(0,0.5))
)
dev.off()
##########
tempsave = cbind(StratMeasure$`n = 2`$`FDR = 0.05`$TRESS$StratTargetP[,3],
                 StratMeasure$`n = 3`$`FDR = 0.05`$TRESS$StratTargetP[,3],
                 StratMeasure$`n = 5`$`FDR = 0.05`$TRESS$StratTargetP[,3],
                 StratMeasure$`n = 7`$`FDR = 0.05`$TRESS$StratTargetP[,3],
                 StratMeasure$`n = 10`$`FDR = 0.05`$TRESS$StratTargetP[,3])
colnames(tempsave) = c("2","3","5","7","10")
tp_s3_0.05 = melt(tempsave)
names(tp_s3_0.05)[2:3] = c("Sample_Size", "TP")
tp_s3_0.05 = tp_s3_0.05[,-1]
tp_s3_0.05$Sample_Size = factor(tp_s3_0.05$Sample_Size, levels = c("2","3","5","7","10"))
pdf(file = paste0("./UpdatedFigures/", "TP_S3_violin_", Model, ".pdf"), height = 5.64)
par(bty="l",lwd=2, mai = c(1.9, 0.82, 0.6, 0.42))
print(ggplot(tp_s3_0.05,
             aes(x= Sample_Size, y = TP, color = Sample_Size)) +
        geom_violin(trim=T, size = 1, width = 1) +
        geom_point(aes(fill = Sample_Size, color = Sample_Size), position = position_jitter(seed = 1, width = 0.2))+
        scale_color_manual(values = group.col)+
        ylab(" ") +
        xlab(" ") +
        theme_bw()+
        theme_classic() +
        theme(axis.text.x = element_text(size = 33, vjust = 0.85, color = "black"),
              axis.text.y = element_text(hjust = 1, size = 33, color = "black"),
              legend.position = "none",
              axis.line = element_line(colour = "black", size = 0.8), 
              axis.ticks = element_line(colour = "black", size = 0.8),
              axis.ticks.length = unit(0.2, "cm")
              #              legend.margin=margin(-30,1,1,1),
              #             legend.box.margin=margin(-5,-5,-5,-5)
        ) +
        ylim(c(0.4, 1))
)
dev.off()

tempsave = cbind(StratMeasure$`n = 2`$`FDR = 0.05`$TRESS$StratFDC[,3],
                 StratMeasure$`n = 3`$`FDR = 0.05`$TRESS$StratFDC[,3],
                 StratMeasure$`n = 5`$`FDR = 0.05`$TRESS$StratFDC[,3],
                 StratMeasure$`n = 7`$`FDR = 0.05`$TRESS$StratFDC[,3],
                 StratMeasure$`n = 10`$`FDR = 0.05`$TRESS$StratFDC[,3])
colnames(tempsave) = c("2","3","5","7","10")
fdc_s3_0.05 = melt(tempsave)
names(fdc_s3_0.05)[2:3] = c("Sample_Size", "FDC")
fdc_s3_0.05 = fdc_s3_0.05[,-1]
fdc_s3_0.05$Sample_Size = factor(fdc_s3_0.05$Sample_Size, levels = c("2","3","5","7","10"))
pdf(file = paste0("./UpdatedFigures/", "FDC_S3_violin_", Model, ".pdf"), height = 5.64)
print(ggplot(fdc_s3_0.05,
             aes(x= Sample_Size, y = FDC, color = Sample_Size)) +
        geom_violin(trim=T, size = 1, width = 1) +
        geom_point(aes(fill = Sample_Size, color = Sample_Size), position = position_jitter(seed = 1, width = 0.2))+
        scale_color_manual(values = group.col)+
        ylab(" ") +
        xlab(" ") +
        theme_bw()+
        theme_classic() +
        theme(axis.text.x = element_text(size = 33, vjust = 0.85, color = "black"),
              axis.text.y = element_text(hjust = 1, size = 33, color = "black"),
              legend.position = "none",
              axis.line = element_line(colour = "black", size = 0.8), 
              axis.ticks = element_line(colour = "black", size = 0.8),
              axis.ticks.length = unit(0.2, "cm")
              #              legend.margin=margin(-30,1,1,1),
              #             legend.box.margin=margin(-5,-5,-5,-5)
        ) +
        ylim(c(0, 0.8))
)
dev.off()
