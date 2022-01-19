
library(devtools)
#install_github("LucianoRogerio/AlphaSimHlpR", force = T, ref = "MarninPkg")
#install_github("wolfemd/genomicMateSelectR", ref = "HEAD", force = T)

suppressMessages(library(AlphaSimHlpR))
suppressMessages(library(tidyverse))
suppressMessages(library(genomicMateSelectR))
library(dplyr); library(furrr); library(here)
select <- dplyr::select
mutate <- dplyr::mutate
schemeDF <- read.csv(here::here("data","baselineScheme - IITA.csv"),
                     header = T, stringsAsFactors = F, sep = ";")

bsp1 <- specifyBSP(schemeDF = schemeDF,
                   nTrainPopCycles = 4, nYrsAsCandidates = 3, maxTrainingPopSize = 1500,
                   nChr = 18, effPopSize = 200, quickHaplo = F,
                   segSites = 1650, nQTL = 1000, nSNP = 500, genVar = 1500,
                   gxeVar = NULL, gxyVar = 1500, gxlVar = 750, gxyxlVar = 300,
                   meanDD = 0.23, varDD = 0.05, relAA = 0.5,
                   stageToGenotype = "CET",
                   nParents = 100, nCrosses = 400, nProgeny = 30, nClonesToNCRP = 3,
                   phenoF1toStage1 = F, errVarPreStage1 = 17500,
                   useCurrentPhenoTrain = F,
                   nCyclesToKeepRecords = 6,
                   selCritPipeAdv = selCritIID,
                   selCritPopImprov =  selCritIID)
source(here::here("code","runBurnInSchemes.R"))

start <- proc.time()[3]
burnIn_IITA_PS <- runBurnInSchemes(bsp = bsp1,
                                   nBurnInCycles=26,
                                   selCritPop="parentSelCritBLUP",
                                   selCritPipe="productSelCritBLUP",
                                   iniFunc="initializeScheme",
                                   productFunc="productPipeline",
                                   popImprovFunc="popImprovByParentSel",
                                   nReplications = 30, nSimCores = 10,
                                   nBLASthreads = 1,nThreadsMacs2 = 1)
end <- proc.time()[3]; print(paste0(floor((end - start)/60), "mins ", floor(((end - start)%%60)), "s elapsed - bsp1"))
saveRDS(burnIn_IITA_PS,file = here::here("output","BurnIn_IITA_PS.rds"))
rm(start); rm(end)

source(here::here("code","runSchemesPostBurnIn.R"))

## This is finished, and it is okay
burnIn_IITA_PS <- readRDS(here::here("output", "BurnIn_IITA_PS.rds"))
start <- proc.time()[3]
postBurnIn_IITA_PS <- runSchemesPostBurnIn(simulations = burnIn_IITA_PS,
                                           nPostBurnInCycles = 12,
                                           selCritPop = "parentSelCritBLUP",
                                           selCritPipe = "productSelCritBLUP",
                                           productFunc = "productPipelinePostBurnIn",
                                           popImprovFunc = "popImprovByParentSel",
                                           nSimCores = 1,
                                           nBLASthreads = 1)
end <- proc.time()[3]; print(paste0(floor((end - start)/60), "mins ", floor(((end - start)%%60)), "s elapsed - bspP"))
saveRDS(postBurnIn_IITA_PS,file = here::here("output","postBurnIn_IITA_PS.rds"))
rm(postBurnIn_IITA_PS)

## This is running, it is okay
start <- proc.time()[3]
postBurnIn_IITA_GS <- runSchemesPostBurnIn(simulations = burnIn_IITA_PS,
                                           nPostBurnInCycles = 12,
                                           selCritPop="parentSelCritGEBV",
                                           selCritPipe="productSelCritBLUP",
                                           productFunc="productPipelinePostBurnIn",
                                           popImprovFunc="popImprovByParentSel",
                                           nSimCores = 1,
                                           nBLASthreads = 1)
end <- proc.time()[3]; print(paste0(floor((end - start)/60), "mins ", floor(((end - start)%%60)), "s elapsed - Newbsp1"))
saveRDS(postBurnIn_IITA_GS,file = here::here("output","postBurnIn_IITA_GS.rds"))
rm(postBurnIn_IITA_GS)


NschemeDF1 <- read.csv(here::here("data","baselineScheme - IITA - NVDP1.csv"),
                       header = T, stringsAsFactors = F, sep = ";")
bsp2 <- specifyBSP(schemeDF = NschemeDF1,
                   nTrainPopCycles = 4, nYrsAsCandidates = 3, maxTrainingPopSize = 1500,
                   nChr = 18, effPopSize = 200, quickHaplo = F,
                   segSites = 1650, nQTL = 1000, nSNP = 500, genVar = 1500,
                   gxeVar = NULL, gxyVar = 1500, gxlVar = 750, gxyxlVar = 300,
                   meanDD = 0.23, varDD = 0.05, relAA = 0.5,
                   stageToGenotype = "CET",
                   nParents = 100, nCrosses = 400, nProgeny = 30, nClonesToNCRP = 3,
                   phenoF1toStage1 = F, errVarPreStage1 = 17500,
                   useCurrentPhenoTrain = F,
                   nCyclesToKeepRecords = 6,
                   selCritPipeAdv = selCritIID,
                   selCritPopImprov =  selCritIID)

## Mitigation used, I just remove the last stage UYT2

#bsp2 <- specifyBSP(schemeDF = schemeDF %>% filter(stageNames!="UYT2"),
#                   nTrainPopCycles = 4, nYrsAsCandidates = 3, maxTrainingPopSize = 1500,
#                   nChr = 18, effPopSize = 200, quickHaplo = F,
#                   segSites = 1650, nQTL = 1000, nSNP = 500, genVar = 1500,
#                   gxeVar = NULL, gxyVar = 1500, gxlVar = 750, gxyxlVar = 300,
#                   meanDD = 0.23, varDD = 0.05, relAA = 0.5,
#                   stageToGenotype = "CET",
#                   nParents = 100, nCrosses = 400, nProgeny = 30, nClonesToNCRP = 3,
#                   phenoF1toStage1 = F, errVarPreStage1 = 17500,
#                   useCurrentPhenoTrain = F,
#                   nCyclesToKeepRecords = 6,
#                   selCritPipeAdv = selCritIID,
#                   selCritPopImprov =  selCritIID)

## This is giving the message "Trying to select invalid individuals",
## when I remove the last stage (UYT2) and increase the number of entries of
## the remaining stages
start <- proc.time()[3]
postBurnIn_NVDP1_IITA_PS <- runSchemesPostBurnIn(simulations = burnIn_IITA_PS,
                                                 newBSP = bsp2, # so you can change the scheme after burn-in
                                                 nPostBurnInCycles = 12,
                                                 selCritPop = "parentSelCritBLUP",
                                                 selCritPipe = "productSelCritBLUP",
                                                 productFunc = "productPipelinePostBurnIn",
                                                 popImprovFunc = "popImprovByParentSel",
                                                 nSimCores = 1,
                                                 nBLASthreads = 1)
end <- proc.time()[3]; print(paste0(floor((end - start)/60), "mins ", floor(((end - start)%%60)), "s elapsed - bspP"))
saveRDS(postBurnIn_NVDP1_IITA_PS,file = here::here("output","postBurnIn_NVDP1_IITA_PS.rds"))
rm(postBurnIn_NVDP1_IITA_PS)
