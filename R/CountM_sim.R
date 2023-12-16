CountM_sim <- function(Input.file,
                      IP.file,
                      BamDir,
                      annoDir,
                      variable,
                      bam_factor,
                      nsim = 10,
                      N.reps = c(2, 3, 5),
                      depth_factor = c(1, 2, 5),
                      thres = c(0.01, 0.05, 0.1, 0.2),
                      dmrProp = NULL) {
  ######## 1. Bam to bin-level data
  if (length(IP.file) == 0) {
    stop("Missing IP samples!",
         call. = TRUE, domain = NULL
    )
  }
  
  if (length(Input.file) == 0) {
    stop("Missing Input samples!",
         call. = TRUE, domain = NULL
    )
  }
  
  if (length(Input.file) != length(IP.file)) {
    stop("IP and Input samples are not paired!",
         call. = TRUE, domain = NULL
    )
  }
  
  if (is.na(annoDir)) {
    stop("Please provide the annotation file!",
         call. = TRUE, domain = NULL
    )
  }
  
  if (length(variable) == 0) {
    stop("Please provide variable for model fitting!",
         call. = TRUE, domain = NULL
    )
  }
  
  if (length(bam_factor) == 0) {
    stop("Please provide a bam_factor!",
         call. = TRUE, domain = NULL
    )
  }
  
  
  variable <- data.frame(Trt = variable)
  model <- ~ 1 + Trt
  
  allBins <- TRESS::DivideBins(
    IP.file = IP.file,
    Input.file = Input.file,
    Path_To_AnnoSqlite = annoDir,
    InputDir = BamDir
  )
  
  ######## 2. Expand bin-level data based on bam_factor
  # allBins = list()
  #
  # allBins$binCount = do.call("rbind",
  #                            replicate(round(1/bam_factor),
  #                                      allBins_user$binCount,
  #                                      simplify = FALSE))
  #
  # allBins$bins = do.call("rbind",
  #                        replicate(round(1/bam_factor),
  #                                  allBins_user$bins,
  #                                  simplify = FALSE))
  
  
  ######### 3. Bin-level data to region-level data
  Candidates <- TRESS::CallCandidates(
    Counts = allBins$binCount,
    bins = allBins$bins
  )
  Candidates <- TRESS::filterRegions(Candidates)
  
  # strata_list = GetStrata(do.call("rbind",
  #                                 replicate(round(1/bam_factor),
  #                                           Candidates$Counts,
  #                                           simplify = FALSE)))
  
  ######### 4. Estimation of parameters for simulation
  message("Estimating parameters...")
  ParaEsti <- CallDMRs.paramEsti(
    counts = do.call(
      "rbind",
      replicate(round(1 / bam_factor),
                Candidates$Counts,
                simplify = FALSE
      )
    ),
    sf = Candidates$sf,
    variable = variable,
    model = model,
    shrkPhi = TRUE
  )
  DMR.test <- TRESS::TRESS_DMRtest(DMR = ParaEsti, contrast = c(0, 1))
  message("Parameters estimation has finished...")
  idx.dmr <- which(DMR.test$padj < 0.05)
  #####adjust DMR prop according to user settings
  if (!is.null(dmrProp)){
    currentProp <- length(idx.dmr)/length(DMR.test$padj)
    if (currentProp > dmrProp){
      newLength <- round(dmrProp * length(DMR.test$padj))
      idx.dmr.new <- sample(idx.dmr, newLength)
    } else {
      idx.dmr.new <- idx.dmr
    }
  } else {
    idx.dmr.new <- idx.dmr
  }
  ########## Extra. KL Selection
  KL.list <- KL_cal(
    N_reps = N.reps,
    ParaEsti = ParaEsti,
    idx.dmr = idx.dmr.new,
    Candidates = Candidates,
    sd_multi = 1,
    nsim = 100
  )
  
  message("KL divergence calculation has finished...")
  
  p.value <- wilcox.test(x = KL.list[[1]][, "NB"], y = KL.list[[1]][, "BB"], alternative = "g")$p.value
  
  if (p.value < 0.05) {
    dist_pref <- "BB"
  } else {
    dist_pref <- "NB"
  }
  
  ######### 5. Simulate data and calculate pvalues and results
  SimCounts_all <- lapply(depth_factor, FUN = function(sd_multi) {
    SimCounts(
      N.reps.Ctrl = N.reps,
      N.reps.Trt = N.reps,
      ParaEsti = ParaEsti,
      idx.dmr = idx.dmr.new,
      nsim = nsim,
      Candidates = Candidates,
      sd_multi,
      model_dist = dist_pref
    )
  })
  
  
  names(SimCounts_all) <- paste0(depth_factor, "x")
  SimCounts_all[["flag"]] <- idx.dmr.new
  # save(Power.list, file = paste0(Outdir, "/", "PowerRes_",ExperimentName, ".rdata"))
  ######### 6. Plot results
  
  # PlotRes(Power.list = Power.list,
  #         Outdir = Outdir,
  #         ExperimentName = ExperimentName)
  # print("The .rdata file and plots are saved.")
  return(SimCounts_all)
}
