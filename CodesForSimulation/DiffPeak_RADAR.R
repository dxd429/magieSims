DiffPeak.RADAR <- function(counts, variable){
  
  ### Counts: order of columns are: 
  #### Ctrl.input_1, ..., Ctrl.input_N.reps[1], Trt.input_1, ..., Trt.input_N.reps[2],
  #### Ctrl.ip_1, ..., Ctrl.ip_N.reps[1], Trt.ip_1, ..., Trt.ip_N.reps[2]
  
  ### variable: contain variables in the design, where the first column corresponds to the factor of interest, or the testing factor. 
  #### The rest columns are covariates
  library(RADAR)
  
  RADAR.dat <- MeRIP()
  
  ### added on Feb 2, 2021 in order to avoid error when transforming MeRIP class to MeRIP.RADAR
  ### this txdb will never be used in simulation
  TxDb.Mmusculus.UCSC.mm9.knownGene::TxDb.Mmusculus.UCSC.mm9.knownGene -> txdb
  RADAR.dat@geneModel = exonsBy(txdb,by="gene")
  ###
  
  
  RADAR.dat@reads = as.matrix(counts) 
  rownames(RADAR.dat@reads) = as.character(paste0("Peak_", 1:nrow(counts) ))
  RADAR.dat@samplenames = as.character(paste0("s", 1:nrow(variable)) )
  
  RADAR.dat <- normalizeLibrary( RADAR.dat, boxPlot = FALSE)
  if(length(levels(variable$predictor)) == 0){
    variable$predictor = as.factor(variable$predictor)
  }
  RADAR.dat <- adjustExprLevel( RADAR.dat )
  variable(RADAR.dat) <- variable
  RADAR.dat <- filterBins(RADAR.dat, minCountsCutOff = 15)
  RADAR.dat <- diffIP(RADAR.dat)
  library(readr)
  
  
  res = matrix(NA, nrow = nrow(counts), ncol = ncol(RADAR.dat@test.est))
  rm.idx = parse_number(setdiff(rownames(RADAR.dat@ip_adjExpr), 
                                rownames(RADAR.dat@ip_adjExpr_filtered)))
  if(length(rm.idx) > 0){
    res[-rm.idx, ] = RADAR.dat@test.est
  }else{
    res = RADAR.dat@test.est
  }
  colnames(res) = colnames(RADAR.dat@test.est) 
  rownames(res) = rownames(RADAR.dat@ip_adjExpr)
  res = as.data.frame(res)
  return(res)
  
  # pval.RADAR = rep(NA, nrow(counts))
  # if(length(rm.idx) > 0){
  #   pval.RADAR[-rm.idx] = RADAR.dat@test.est[, "p_value"]
  # }else{
  #   pval.RADAR = RADAR.dat@test.est[, "p_value"]
  # }
  # return(pval.RADAR)
}
