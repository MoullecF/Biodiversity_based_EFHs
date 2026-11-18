# ------------------------------------------------------------------------------
# Fit abundance model using HPC cluster with GPUs and high RAM
# ------------------------------------------------------------------------------

# Load libraries
source("Scripts_R/0-Load libraries.R")

# Settings
nSamples <- 250
thin <- 4000
nChains <- 4
verbose <- 10
transient <- ceiling(0.5*nSamples*thin)

modelnames <- "normal"
filename_abu <- paste("../Models/init_file_abu_thin_", as.character(thin), 
                     "_samples_", as.character(nSamples),
                     "_chains_", as.character(nChains),
                     ".rds", sep = "")

load(file = "../Models/unfitted_models")

# ------------------------------------------------------------------------------
# Create initialization object for HPC cluster
# ------------------------------------------------------------------------------
# This step prepares the model for fitting on a High-Performance Computing (HPC)
# cluster with GPU acceleration and high RAM. The init_obj_abu contains all
# necessary information to run MCMC sampling remotely.

### Abundance model
init_obj_abu <- sampleMcmc(models[[2]], samples = nSamples, thin = thin,
                          adaptNf = rep(transient, models[[2]]$nr),
                          transient = transient, 
                          nChains = nChains,
                          nParallel = nChains,
                          useSocket = FALSE,
                          verbose = verbose,
                          initPar = "fixed effects",
                          engine = "HPC")

# Save initialization file for HPC submission
init_file_path <- file.path(filename_abu)
saveRDS(to_json(init_obj_abu), file = init_file_path)

# ------------------------------------------------------------------------------
# HPC cluster job submission
# ------------------------------------------------------------------------------
# Once the init_file is created, submit the job to the HPC cluster.
# The cluster will use GPUs and high RAM to perform MCMC sampling efficiently.
# Use the SLURM script "run_hmsc_abu.slurm" located in the "models" folder
# to submit the job on the Jean Zay HPC cluster.
# After completion, retrieve the posterior chain files from the "models" folder.

# ------------------------------------------------------------------------------
# Import calculated posterior samples from HPC
# ------------------------------------------------------------------------------

chainList <- vector("list", nChains)
for(cInd in 1:nChains){
  chain_file_path <- file.path(paste("../Models/abu", "_thin_", as.character(thin), "_samples_", as.character(nSamples), "_chains_", as.character(nChains), "_post_chain0", cInd-1, "_file.rds", sep = ""))
  chainList[[cInd]] <- from_json(readRDS(file = chain_file_path)[[1]])[[1]]
}

# Combine chains into fitted model object
fitSepTF_abu <- importPosteriorFromHPC(models[[2]], chainList, nSamples, thin, transient)

# ------------------------------------------------------------------------------
# Save posterior samples
# ------------------------------------------------------------------------------

filename <- paste("../Models/model_abu_thin_", as.character(thin),
                 "_samples_", as.character(nSamples),
                 "_chains_", as.character(nChains),
                 ".RData", sep = "")

save(fitSepTF_abu, file = filename)
