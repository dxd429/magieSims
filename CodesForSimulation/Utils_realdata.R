#### Utils_Realdata
FilterLowCountRegion <-function(Candidates, Quant = 0.99){
  #### filtering step: rm regions with extremely small input, or extremely high input count, 
  #### or extremely high IP count
  library(matrixStats)
  counts.norm = sweep(Candidates$Counts, 2, Candidates$sf, FUN = "/")
  X.norm = round(counts.norm[, seq(1, ncol(Candidates$Counts), 2)])+1
  Y.norm = round(counts.norm[, seq(2, ncol(Candidates$Counts), 2)])+1
  
  i.0 = which(rowSums(counts.norm[, seq(1, ncol(counts.norm), 2)]) > 20 )
  Q.x = quantile(X.norm, prob = Quant); i.x = which( rowSums(X.norm <= Q.x)== ncol(X.norm) )
  Q.y = quantile(Y.norm, prob = Quant); i.y = which( rowSums(Y.norm <= Q.y)== ncol(Y.norm) )
  idx = intersect(i.0, intersect(i.x, i.y))
  
  Candidates$Counts = Candidates$Counts[idx, ]; rownames(Candidates$Counts) = 1:nrow(Candidates$Counts)
  Candidates$Regions = Candidates$Regions[idx, ];rownames(Candidates$Regions) = 1:nrow(Candidates$Regions)
  
  return(Candidates)
}

CompareSimRealCount <- function(Sim.Candidates, Real.Candidates,
                                input.lim = 1000, ip.lim = 1500,
                                PLOT = TRUE, outputPDF = FALSE){
  #### compare simulated and real candidate count density
  sim.X = sweep(Sim.Candidates$counts[, seq(1, ncol(Sim.Candidates$counts), by = 2)], 2, 
                Sim.Candidates$sf[seq(1, ncol(Sim.Candidates$counts), by = 2)], FUN = "/" )
  sim.Y = sweep(Sim.Candidates$counts[, seq(2, ncol(Sim.Candidates$counts), by = 2)], 2, 
                Sim.Candidates$sf[seq(2, ncol(Sim.Candidates$counts), by = 2)], FUN = "/" )
  Real.X = sweep(Real.Candidates$Counts[, seq(1, ncol(Real.Candidates$Counts), 2)], 2,
                 Real.Candidates$sf[seq(1, ncol(Real.Candidates$Counts), 2)], FUN = "/")
  Real.Y = sweep(Real.Candidates$Counts[, seq(2, ncol(Real.Candidates$Counts), 2)], 2,
                 Real.Candidates$sf[seq(2, ncol(Real.Candidates$Counts), 2)], FUN = "/")
  
  if(PLOT){
   # par(mfrow = c(3, 2))
  #  hist(log(sim.X), 40, probability = TRUE, xlim = c(0, 8), ylim = c(0, 0.6))
  #  hist(log(sim.Y), 40, probability = TRUE, xlim = c(0, 8), ylim = c(0, 0.5))
  #  hist(log(Real.X), 40, probability = TRUE, xlim = c(0, 8), ylim = c(0, 0.6))
  #  hist(log(Real.Y), 30, probability = TRUE, xlim = c(0, 8), ylim = c(0, 0.5))
    par(mfrow = c(1, 2))
    qqplot(x = rowMeans(Real.X), y = rowMeans(sim.X), pch = 16, cex = 1,
           ylab = "", xaxt = "n", yaxt = "n", xlab = " ", main = "Q-Q Plot",
           xlim=c(0,input.lim), ylim = c(0,input.lim))
    axis(1,
         cex.axis = 2.5, lwd = 2, tck = -0.015)
    axis(2,
         cex.axis = 2.5, lwd = 2, tck = -0.015)
    abline(0, 1, col = "red", lwd = 1.5)
    
    qqplot(x = rowMeans(Real.Y), y = rowMeans(sim.Y), pch = 16, cex = 1,
           ylab = "", xaxt = "n", yaxt = "n", xlab = " ", main = "Q-Q Plot",
           xlim=c(0,ip.lim), ylim = c(0,ip.lim))
    axis(1,
         cex.axis = 2.5, lwd = 2, tck = -0.015)
    axis(2,
         cex.axis = 2.5, lwd = 2, tck = -0.015)
    abline(0,1,col = "red", lwd = 1.5)
  }
  
  #### KL divergence
 KL =  (KL.qp(Real.dist = Real.X, Sim.dist = sim.X) + KL.qp(Real.dist = Real.Y, Sim.dist = sim.Y))/2
 cat("The KL divergence of simulated count to real count is ", as.numeric(KL), sep = "\n")
 return(KL)
}


KL.qp <- function(Real.dist, Sim.dist){
  Real.range = range(log(Real.dist[Real.dist > 1]));
  print(Real.range);
  int = seq(Real.range[1], Real.range[2], by = 0.1)
  if(max(int) > Real.range[2]){
    xbrks = int
  } else{
    xbrks = c(int, max(int) + 0.1)
  }
  p.real = hist(log(Real.dist[Real.dist > 1]), breaks = xbrks, 
                plot = FALSE)$density*0.1
  sum(p.real)
  p.real[p.real == 0] = min(p.real[p.real > 0])/2
  p.real = p.real/sum(p.real)
  min(p.real)
  
  #####
  Sim.dist = Sim.dist[Sim.dist > 1]
  Sim.dist = Sim.dist[ log(Sim.dist)>= Real.range[1] & log(Sim.dist) <= Real.range[2] ]
  q.sim = hist(log(Sim.dist), breaks = xbrks,  plot = FALSE)$density*0.1
  
  q.sim[q.sim == 0] = min(q.sim[q.sim>0])/2
  q.sim = q.sim/sum(q.sim)
  min(q.sim)
  KL.sTOr = sum(p.real*log(p.real/q.sim))
  
  par(mfrow = c(1, 1))
  par(mgp = c(3,1.5,0))
  plot(x = xbrks[-length(xbrks)], y = p.real, col = "pink3", pch =1, xlim = c(0, max(xbrks[-length(xbrks)])), ylab = "", xaxt = "n", yaxt = "n", xlab = " ", lwd=2)
  axis(1,at = c(0,1,2,3,4,5,6,7,8),
       labels = c("0",  "1",  "2",  "3", "4", "5",  "6",  "7",  "8"),
       cex.axis = 2.5, lwd = 2, tck = -0.015)
  axis(2,
       cex.axis = 2.5, lwd = 2, tck = -0.015)
  
  points(x = xbrks[-length(xbrks)], y = q.sim, col = "skyblue3", #"skyblue3",
         pch = 16, cex = 1)
  legend("topleft", legend = c("Real", "Simulated"), 
         col = c("pink3", "skyblue3"), 
         pch = c(1, 16), pt.lwd = c(2, 1),bty = "n")
  
  return(KL.sTOr)
}


# ####### theta by relationship with phi.NB in MLE
# coef.lm = coef(lm(log(allPara$theta) ~ log(allPara$Dispersion)))
# set.seed(12345)
# theta = exp(rnorm(length(phi), mean = coef.lm[1] + coef.lm[2]*log(res.BB$PostPhi.BB), 
#                   sd = 0.1))
# plot(log(theta),log(res.theta$theta), cex = 0.5, col = "#00000030")
# abline(0,1,col = "red")
# ########

SelectSimStrategyByKL <- function(allPara,dataset,
                                  nreps, nsites, p.prop, 
                                  PLOT = FALSE,
                                  nsim = 100){
  #####
  SampleInfo = table(unlist(lapply(strsplit(colnames(allPara$Candidates$Counts), split = "_"  ),
                                   function(x) x[1])))
  SampleName = names(SampleInfo)
  
  design = data.frame(predictor = rep(SampleName, nreps))
  model = ~1 + predictor
  
  strategy = c("NB", "BB")
  
  KL = matrix(NA, nrow = nsim, ncol = length(strategy))
  colnames(KL) = strategy
  for (ist in 1:length(strategy)) {
     cat("Strategy ", strategy[ist], sep = "\n")
    res.para = SimPara.TRESS(nreps = nreps,nsites = nsites, p.prop = p.prop,
                             EstiPara = allPara, design = design,
                             model = model, adjTheta = TRUE,
                             PhiStrategy = strategy[ist])
    for (isim in 1:nsim) {
      res.sim = SimDat.TRESS(nreps = nreps,  nsites = length(res.para$flag), 
                             mu = res.para$mu, phi = res.para$phi, 
                             theta = res.para$theta, 
                             seed = 12344+isim)
      KL[isim, ist] = CompareSimRealCount(Sim.Candidates = res.sim,
                                          Real.Candidates = allPara$Candidates,
                                          PLOT = FALSE)
    }
  }
  
  p.value = wilcox.test(x = KL[, "NB"], y = KL[, "BB"], alternative = "g")$p.value
  
  if(PLOT){
   # pdf(file = paste0(getwd(), "/Results/", dataset, "_KL.pdf"), width = 6, height = 5)
    boxplot(KL, ylab = "KL divergence")
    if(p.value < 0.0001){
      legend("bottomleft", legend = paste0("P-value < 0.0001"), bty = "n")
    }else{
      legend("bottomleft", legend = paste0("P-value = ", round(p.value, 4)), bty = "n")
    }
    title(dataset, line = 0.5)
  #  dev.off()
  }
  write.csv(KL, file = paste0("./Results/KL_", dataset, ".csv"))
  
  if(p.value < 0.05){
    return("BB")
  }else{
    return("NB")
  }
}






adjSimTheta <- function(theta.ini, Sim.Candidates, Real.Candidates,
                        Inflate.Input = 1, 
                        Inflate.IP = 1){
  
  SimsInfo = table(unlist(lapply(strsplit(colnames(Sim.Candidates$counts), split = "_"  ), 
                                 function(x) x[1])))
  SimsName = names(SimsInfo)
  nreps = as.numeric(SimsInfo)/2
  
  #####
  sim.X = sweep(Sim.Candidates$counts[, seq(1, ncol(Sim.Candidates$counts), by = 2)], 2, 
                Sim.Candidates$sf[seq(1, ncol(Sim.Candidates$counts), by = 2)], FUN = "/" )
  sim.X.Ctrl = sim.X[, grepl(SimsName[1], colnames(sim.X))]
  sim.X.Trt = sim.X[, grepl(SimsName[2], colnames(sim.X))]
  
  sim.Y = sweep(Sim.Candidates$counts[, seq(2, ncol(Sim.Candidates$counts), by = 2)], 2, 
                Sim.Candidates$sf[seq(2, ncol(Sim.Candidates$counts), by = 2)], FUN = "/" )
  sim.Y.Ctrl = sim.Y[, grepl(SimsName[1], colnames(sim.Y))]
  sim.Y.Trt = sim.Y[, grepl(SimsName[2], colnames(sim.Y))]
  
  ## Real data
  SampleInfo = table(unlist(lapply(strsplit(colnames(Real.Candidates$Counts), split = "_"  ), 
                                   function(x) x[1])))
  SampleName = names(SampleInfo)
  
  Real.X = sweep(Real.Candidates$Counts[, seq(1, ncol(Real.Candidates$Counts), 2)], 2,
                 Real.Candidates$sf[seq(1, ncol(Real.Candidates$Counts), 2)], FUN = "/")
  Real.X.Ctrl = Real.X[, grepl(SampleName[1], colnames(Real.X))]
  Real.X.Trt = Real.X[, grepl(SampleName[2], colnames(Real.X))]
  
  Real.Y = sweep(Real.Candidates$Counts[, seq(2, ncol(Real.Candidates$Counts), 2)], 2,
                 Real.Candidates$sf[seq(2, ncol(Real.Candidates$Counts), 2)], FUN = "/")
  Real.Y.Ctrl = Real.Y[, grepl(SampleName[1], colnames(Real.Y))]
  Real.Y.Trt = Real.Y[, grepl(SampleName[2], colnames(Real.Y))]
  

  ######################################################
  ##########
  myfun = mean
  s.Input = myfun(log(rowMeans(Real.X + 1)))/myfun(log(rowMeans(sim.X + 1)))
  s.IP = myfun(log(rowMeans(Real.Y+1)))/myfun(log(rowMeans(sim.Y+1)))
  
  s.Input = rep(s.Input*Inflate.Input, sum(nreps))
  s.IP = rep(s.IP*Inflate.IP, sum(nreps))
  ######################################################
  
  
  adj.theta = theta.ini
  adj.theta$Input = sweep(theta.ini$Input, 2, s.Input, FUN = "*")
  adj.theta$IP = sweep(theta.ini$IP, 2, s.IP, FUN = "*")
  
  return(adj.theta)
}



Strata_ByDMRInput <- function(res.sim, P.adj){
  #### Statified presentation: obtain strata according to the Input count of DMR
  meanInput = rowMeans(res.sim$counts[, seq(1, ncol(res.sim$counts), 2)])
  sigInput = rowMeans(res.sim$counts[which(P.adj< 0.05),
                                     seq(1, ncol(res.sim$counts), 2)])
  
  strat = quantile(sigInput, prob = seq(0.1, 1, by = 0.2))
  Int.L = c(0, strat)
  Int.U = c(strat, max(res.sim$counts[, seq(1, ncol(res.sim$counts), 2)]))
  Stratums = rbind(Int.L, Int.U)
  colnames(Stratums) = c("(0, 10%]", "(10%, 30%]", "(30%, 50%]", "(50%, 70%]",
                         "(70%, 90%]","(90%, 100%]")
  return(Stratums)
}
