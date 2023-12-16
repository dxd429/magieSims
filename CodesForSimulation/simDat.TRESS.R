SimDat.TRESS <- function(nreps, 
                         nsites, 
                         mu, 
                         phi,
                         theta,
                         seed = 12345){
  ### Directly sample from NB model.
  
  #### nreps: a vector of length G, containing the number of replicate in all G groups 
  #### mu: nsite*sum(nreps) matrix, denoting the methylation level in all samples, 
  ####     the order of samples are group 1, group2, group3,...
  #### phi: a vector of length nsites, denoting the dispersion of methylation for all sites
  #### theta: a list containing scale parameter of gamma for all sites, cross all conditions
  
  set.seed(seed = seed)
  ######
  sf.x = runif(sum(nreps), 0.5, 2)
  sf.y = runif(sum(nreps), 0.5, 2)
  size.x = as.vector((1-mu)*(phi^{-1} - 1))
  prob.x = as.vector(1/(1+sweep(theta$Input, 2, sf.x , FUN = "*")))
  count.Input = matrix(rnbinom(n = nsites*sum(nreps),size =size.x, prob = prob.x), 
                       nrow = nsites, byrow = FALSE)
  
  set.seed(seed = seed)
  size.y = as.vector(mu*(phi^{-1} - 1))
  prob.y = as.vector(1/(1+sweep(theta$IP, 2, sf.y, FUN = "*")))
  count.IP = matrix(rnbinom(n = nsites*sum(nreps),size =size.y, prob = prob.y), 
                    nrow = nsites, byrow = FALSE)
  ######
  counts = matrix(0, nrow = nsites, ncol = 2*sum(nreps))
  counts[, seq(1, ncol(counts), by = 2)] = count.Input
  counts[, seq(2, ncol(counts), by = 2)] = count.IP
  num.na = sum(is.na(counts))
  if(length(num.na) > 0)
    counts[is.na(counts)] = rpois(num.na, 1e+05)
  
  colnames(counts) = 1:ncol(counts)
  colnames(counts)[seq(1, ncol(counts), 2)] = c(paste0("Ctrl_Input_Rep_",  1:nreps[1]),
                                                paste0("Trt_Input_Rep_",  1:nreps[2]))
  colnames(counts)[seq(2, ncol(counts), 2)] = c(paste0("Ctrl_IP_Rep_",  1:nreps[1]),
                                                paste0("Trt_IP_Rep_",  1:nreps[2]))
  sf = rep(0, 2*sum(nreps))
  sf[seq(1, ncol(counts), 2)] = sf.x
  sf[seq(2, ncol(counts), 2)] = sf.y
  return(list(counts = counts, sf = sf))
}


SimDat.TRESS.Coverage <- function(nreps, nsites, 
                                  mu, 
                                  phi,
                                  theta,
                                  Coverage = 1,
                                  seed = 123){
  ### simulate data based on mu, phi and theta
  #### nreps: a vector of length G, containing the number of replicate in all G groups 
  #### mu: nsite*sum(nreps) matrix, denoting the methylation level in all samples, 
  ####     the order of samples are group 1, group2, group3,...
  #### phi: a vector of length nsites, denoting the dispersion of methylation for all sites
  #### theta: a list containing scale parameter of gamma for all sites, cross all conditions
  
  
  set.seed(seed = seed)
  ##### use mu to generate lambda
  call.rgamma <- function(x, n){
    #rgamma(n, shape = x[1:n], scale = x[length(x)])
    rgamma(n, shape = x[1:n], scale = x[(n+1):(2*n)])
  }
  
  plx = cbind((1- mu)*(1/phi - 1), theta$Input)
  lambda_x = t(apply(plx, 1, call.rgamma, n = sum(nreps)))*Coverage
  
  ply = cbind(mu*(1/phi - 1), theta$IP)
  lambda_y = t(apply(ply, 1, call.rgamma, n = sum(nreps)))*Coverage
  
  ### size factor
  sf.x = runif(sum(nreps), 0.5, 2)
  sf.y = runif(sum(nreps), 0.5, 2)
  ### sample IP and control based on lambda_x, and lambda_y
  para.x = sweep(lambda_x, 2, sf.x, FUN = "*")
  para.y = sweep(lambda_y, 2, sf.y, FUN = "*")
  ### input ~ poison(para.x), 
  ### IP ~ poison(para.y)
  count.Input = matrix(rpois(n = nsites*sum(nreps), as.vector(para.x)), 
                       nrow = nsites, byrow = FALSE)
  count.IP = matrix(rpois(n = nsites*sum(nreps), as.vector(para.y)),
                    nrow = nsites, byrow = FALSE)
  

  ######
  counts = matrix(0, nrow = nsites, ncol = 2*sum(nreps))
  counts[, seq(1, ncol(counts), by = 2)] = count.Input
  counts[, seq(2, ncol(counts), by = 2)] = count.IP
  num.na = sum(is.na(counts))
  if(length(num.na) > 0)
    counts[is.na(counts)] = rpois(num.na, 1e+05)
  
  colnames(counts) = 1:ncol(counts)
  colnames(counts)[seq(1, ncol(counts), 2)] = c(paste0("Ctrl_Input_Rep_",  1:nreps[1]),
                                                paste0("Trt_Input_Rep_",  1:nreps[2]))
  colnames(counts)[seq(2, ncol(counts), 2)] = c(paste0("Ctrl_IP_Rep_",  1:nreps[1]),
                                                paste0("Trt_IP_Rep_",  1:nreps[2]))
  
  sf = rep(0, 2*sum(nreps))
  sf[seq(1, ncol(counts), 2)] = sf.x
  sf[seq(2, ncol(counts), 2)] = sf.y
  return(list(counts = counts, sf = sf, 
              lambda_y = lambda_y, lambda_x = lambda_x,
              para.x = para.x, para.y = para.y))
}


