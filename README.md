
# Workflow to obtain simulated data and apply DMR detection methods for power evaluation

This workflow acts as a complement to our R package `magpie`, providing an alternative way for users who wish to conduct power analysis using their DMR detection methods of choice. `magpie` currently includes the DMR testing methods `TRESS` and `exomePeak2`, and is available at https://bioconductor.org/packages/magpie/.

Below, we present a detailed breakdown of data simulation and DMR testing.

## Getting Started

### Load Necessary Libraries and Source Local Functions

First, it's essential to load the necessary libraries and source the local functions for the workflow.

### Installing Example Data

To access the example data, the data package `magpieData` needs to be installed from GitHub. Use the `devtools` package to do this:

```R
library(devtools)
install_github("dxd429/magpieData")
```

### Library Loading

Proceed by loading the following libraries:

```R
library(miceadds)
source.all("./R")
library(magpieData)
library(magpie)
library(BiocParallel)
library(Matrix)
library(matrixStats)
library(reshape2)
library(RColorBrewer)
```

### Retrieve Example Data

Retrieve the example data from the `magpieData` package using:

```R
BAM_path <- getBAMpath()
```

## Simulation Parameters

Set the parameters required for simulation. Feel free to adjust these parameters according to your needs:

```R
bam_factor <- 0.03
nsim <- 1
N.reps <- c(2, 3, 5, 7)
depth_factor <- c(1, 2)
thres <- c(0.01, 0.05, 0.1)
```

## Simulate Count Matrices for Testing

To simulate count matrices, use the following code:

```R
Dat_sims <- CountM_sim(
   Input.file = c("Ctrl1.chr15.input.bam", "Ctrl2.chr15.input.bam", "Case1.chr15.input.bam", "Case2.chr15.input.bam"),
    IP.file = c("Ctrl1.chr15.ip.bam", "Ctrl2.chr15.ip.bam", "Case1.chr15.ip.bam", "Case2.chr15.ip.bam"),
    BamDir = BAM_path,
    annoDir = paste0(BAM_path, "/hg18_chr15.sqlite"),
    variable = rep(c("Ctrl", "Trt"), each = 2),
    bam_factor = bam_factor,
    nsim = nsim,
    N.reps = N.reps,
    depth_factor = depth_factor,
    thres = thres
 )
```
Here, `Dat_sims` is a list containing two main elements, each representing simulated count matrices at different sequencing depths: `1x` and `2x`. Within each of these elements, there are four sub-lists corresponding to varying sample sizes per group—specifically 2, 3, 5, and 7, parameters defined in `Simulation Parameters`. For each combination of sequencing depth and sample size, currently, there's only one matrix since we set `nsim <- 1`, indicating a single iteration for each scenario. Had `nsim` been assigned a different value, the number of matrices in each sub-list would align with the value of `nsim`.

## Apply Testing Method

For illustration purposes, we apply TRESS in this section. It's crucial to adjust your DMR detection function accordingly so that the input is a count matrix with rows representing regions and parameters related to the experimental design. The output should ideally be p-values or adjusted p-values, as used here:

```R
PVALS <- vector("list", length = length(N.reps))
names(PVALS) <- paste0("n = ", N.reps)
for (i in 1:length(Dat_sims[[1]])){
  p_df <- NULL
  for (j in length(Dat_sims[[1]][[i]])){
    res.sim <- Dat_sims[[1]][[i]][[j]]
    design <- data.frame(Trt = rep(c("Ctrl", "Trt"), each = N.reps[i]))
    TRESS.DMR <- CallDMRs.paramEsti(
      counts = res.sim$counts,
      sf = res.sim$sf,
      variable = design,
      model = ~ 1 + Trt
    )
    Contrast <- c(c(0, 1), rep(0, ncol(design) - 1))
    res.test <- TRESS::TRESS_DMRtest(DMR = TRESS.DMR, contrast = Contrast)
    p_df <- data.frame(cbind(p_df, res.test$padj))
  }
  PVALS[[i]] <- p_df
}
```

## Investigation and Plotting

After applying the testing method, you can either perform your own analysis on the testing results or utilize the plotting functions available in the `magpie` package to create insightful line plots. The following code can be used for calculating power and generating a line plot for `FDR`:

```R
Power.list <- Power.cal(
  PVALS = PVALS,
  idx.dmr = Dat_sims$flag,
  N.reps = N.reps ,
  thre = thres
)
Power.list <- list("1x" = Power.list)
plotRes(Power.list, value_option = "FDR")
```

![FDR](images/fdr.pdf)
