SimPara.TRESS <- function(nreps, nsites, p.prop,
                          EstiPara,
                          design, model,
                          PhiStrategy = "NB",
                          adjTheta = TRUE,
                          extraInf.theta.Input = 1,
                          extraInf.theta.IP = 1,
                          seed = 12345){
  
  
  ###
  set.seed(seed)
  nsites = min(nrow(allPara$Coef), nsites)
  idx = sample(1:nrow(allPara$Coef), nsites)
  sim.Coef = allPara$Coef[idx, ]
  
  if(PhiStrategy == "NB"){
    sim.phi = allPara$Dispersion[idx]
    sim.theta = allPara$theta[idx]
    if(mean(log(allPara$Dispersion[idx]), na.rm = TRUE) > -3){
      sim.phi = exp(log(allPara$Dispersion[idx])-2.5)
      sim.theta = allPara$theta[idx]*exp(-2.5)
    }
  }else if(PhiStrategy == "BB"){
    sim.phi = EstiPara$Dispersion.BB[idx]
    sim.theta = EstiPara$theta.BB[idx]
  }

  
  ##### generate new mu
  set.seed(seed = seed)
  flag = rep(FALSE, nsites)
  idx.TP = sample(1:nsites, size = round(nsites*p.prop))
  flag[idx.TP] = TRUE
  
  #### alpha and beta
  set.seed(seed = seed)
  sim.alpha = sim.Coef[, 1]
  
  ### 1. simulate beta from uniform distribution
  # set.seed(seed = seed)
  # sim.beta = rep(0, nsites)
  # sim.beta[idx.TP[1:(nsites*p.prop/2)]] = runif(nsites*p.prop/2, min = 1, max = 2)
  # sim.beta[idx.TP[(nsites*p.prop/2 + 1):(nsites*p.prop)]] = runif(nsites*p.prop/2,
  #                                                                 min = -2, max = -1)
 
  ### 2. sample from originally estimated beta which are relatively large 
  set.seed(seed = seed)
  Q.beta = quantile(abs(sim.Coef[, 2]), prob = c(0.50, 0.95), na.rm= TRUE)
  id.candi = which(abs(sim.Coef[, 2]) > Q.beta[1] & abs(sim.Coef[, 2]) < Q.beta[2])
  iii = intersect(id.candi, idx.TP)
  
  sim.beta = rep(0, nsites)
  sim.beta[iii] = sim.Coef[iii, 2]
  sim.beta[setdiff(idx.TP, iii)] = runif(length(setdiff(idx.TP, iii)), 1, 2)
  ####
  
  sim.eta = t(model.matrix(model, design)%*% t(cbind(sim.alpha, sim.beta))) 
  sim.mu = exp(sim.eta)/(1 + exp(sim.eta))
  sim.mu[sim.mu < 1e-5] = 0.01; sim.mu[sim.mu > (1- (1e-5))]=0.99
  

  
  Theta0.CtrlInput = Theta0.TrtInput = 
    Theta0.CtrlIP = Theta0.TrtIP = sim.theta
  
  sim.theta = list(Input = cbind( matrix(rep(Theta0.CtrlInput, nreps[1]), 
                                         ncol = nreps[1], byrow = FALSE), 
                                  matrix(rep(Theta0.TrtInput, nreps[2]), 
                                         ncol = nreps[2], byrow = FALSE)),
                   IP = cbind( matrix(rep(Theta0.CtrlIP, nreps[1]), 
                                      ncol = nreps[1], byrow = FALSE), 
                               matrix(rep(Theta0.TrtIP, nreps[2]), 
                                      ncol = nreps[2], byrow = FALSE)))
  colnames(sim.theta$Input) = c(paste0("Ctrl_input_rep", 1:nreps[1]),
                                paste0("Trt_input_rep", 1:nreps[2]))
  colnames(sim.theta$IP) = c(paste0("Ctrl_IP_rep", 1:nreps[1]),
                             paste0("Trt_IP_rep", 1:nreps[2]))
  
  res = list(flag = flag, mu = sim.mu, phi = sim.phi, theta = sim.theta, alpha = sim.alpha, beta = sim.beta)
  
  
  if(adjTheta){
    ############################
    tmp.sim = SimDat.TRESS(nreps = nreps, nsites = nsites, 
                           mu = sim.mu, phi = sim.phi, 
                           theta = sim.theta, 
                           seed = seed)
    adjtheta = adjSimTheta(theta.ini = sim.theta, 
                           Sim.Candidates = tmp.sim, 
                           Real.Candidates = EstiPara$Candidates,
                           Inflate.Input = extraInf.theta.Input,
                           Inflate.IP = extraInf.theta.IP)
    
    res = list(flag = flag, mu = sim.mu, phi = sim.phi, 
               theta = adjtheta, alpha = sim.alpha, 
               beta = sim.beta)
  }
  
  return(res)
}


SimPara.TRESS_uniAlpha <- function(nreps, nsites, p.prop,
                          EstiPara,
                          design, model,
                          PhiStrategy = "NB",
                          adjTheta = TRUE,
                          extraInf.theta.Input = 1,
                          extraInf.theta.IP = 1,
                          seed = 12345){
  
  
  ###
  set.seed(seed)
  nsites = min(nrow(allPara$Coef), nsites)
  idx = sample(1:nrow(allPara$Coef), nsites)
  sim.Coef = allPara$Coef[idx, ]
  
  if(PhiStrategy == "NB"){
    sim.phi = allPara$Dispersion[idx]
    sim.theta = allPara$theta[idx]
    if(mean(log(allPara$Dispersion[idx]), na.rm = TRUE) > -3){
      sim.phi = exp(log(allPara$Dispersion[idx])-2.5)
      sim.theta = allPara$theta[idx]*exp(-2.5)
    }
  }else if(PhiStrategy == "BB"){
    sim.phi = EstiPara$Dispersion.BB[idx]
    sim.theta = EstiPara$theta.BB[idx]
  }
  
  
  ##### generate new mu
  set.seed(seed = seed)
  flag = rep(FALSE, nsites)
  idx.TP = sample(1:nsites, size = round(nsites*p.prop))
  flag[idx.TP] = TRUE
  
  #### alpha and beta
  set.seed(seed = seed)
  sim.alpha_all = sim.Coef[, 1]
  sim.alpha = rep(mean(sim.alpha_all, na.rm = TRUE), length(sim.alpha_all))
  
  ### 1. simulate beta from uniform distribution
  # set.seed(seed = seed)
  # sim.beta = rep(0, nsites)
  # sim.beta[idx.TP[1:(nsites*p.prop/2)]] = runif(nsites*p.prop/2, min = 1, max = 2)
  # sim.beta[idx.TP[(nsites*p.prop/2 + 1):(nsites*p.prop)]] = runif(nsites*p.prop/2,
  #                                                                 min = -2, max = -1)
  
  ### 2. sample from originally estimated beta which are relatively large 
  set.seed(seed = seed)
  Q.beta = quantile(abs(sim.Coef[, 2]), prob = c(0.50, 0.95), na.rm= TRUE)
  id.candi = which(abs(sim.Coef[, 2]) > Q.beta[1] & abs(sim.Coef[, 2]) < Q.beta[2])
  iii = intersect(id.candi, idx.TP)
  
  sim.beta = rep(0, nsites)
  sim.beta[iii] = sim.Coef[iii, 2]
  sim.beta[setdiff(idx.TP, iii)] = runif(length(setdiff(idx.TP, iii)), 1, 2)
  ####
  
  sim.eta = t(model.matrix(model, design)%*% t(cbind(sim.alpha, sim.beta))) 
  sim.mu = exp(sim.eta)/(1 + exp(sim.eta))
  sim.mu[sim.mu < 1e-5] = 0.01; sim.mu[sim.mu > (1- (1e-5))]=0.99
  
  
  
  Theta0.CtrlInput = Theta0.TrtInput = 
    Theta0.CtrlIP = Theta0.TrtIP = sim.theta
  
  sim.theta = list(Input = cbind( matrix(rep(Theta0.CtrlInput, nreps[1]), 
                                         ncol = nreps[1], byrow = FALSE), 
                                  matrix(rep(Theta0.TrtInput, nreps[2]), 
                                         ncol = nreps[2], byrow = FALSE)),
                   IP = cbind( matrix(rep(Theta0.CtrlIP, nreps[1]), 
                                      ncol = nreps[1], byrow = FALSE), 
                               matrix(rep(Theta0.TrtIP, nreps[2]), 
                                      ncol = nreps[2], byrow = FALSE)))
  colnames(sim.theta$Input) = c(paste0("Ctrl_input_rep", 1:nreps[1]),
                                paste0("Trt_input_rep", 1:nreps[2]))
  colnames(sim.theta$IP) = c(paste0("Ctrl_IP_rep", 1:nreps[1]),
                             paste0("Trt_IP_rep", 1:nreps[2]))
  
  res = list(flag = flag, mu = sim.mu, phi = sim.phi, theta = sim.theta, alpha = sim.alpha, beta = sim.beta)
  
  
  if(adjTheta){
    ############################
    tmp.sim = SimDat.TRESS(nreps = nreps, nsites = nsites, 
                           mu = sim.mu, phi = sim.phi, 
                           theta = sim.theta, 
                           seed = seed)
    adjtheta = adjSimTheta(theta.ini = sim.theta, 
                           Sim.Candidates = tmp.sim, 
                           Real.Candidates = EstiPara$Candidates,
                           Inflate.Input = extraInf.theta.Input,
                           Inflate.IP = extraInf.theta.IP)
    
    res = list(flag = flag, mu = sim.mu, phi = sim.phi, 
               theta = adjtheta, alpha = sim.alpha, 
               beta = sim.beta)
  }
  
  return(res)
}




