#########################
rm(list = ls())
resPath <- "~/scratch-data/neta-metabric/disc-cv-vote/her2/04/"
load(file.path(resPath, "her2-04.rda"))
source("~/repos/spacemap/R/cvVoteRDS.R")
minIndex <- cvsmap$minIndex
cvsmap$bestFoldFiles <- list.files(path = resPath, 
                                   pattern =  paste0("tuneid_", sprintf("%03d", minIndex), "_fold"),
                                   full.names = TRUE)
cvsmap$cvVote <- cvVoteRDS(files = cvsmap$bestFoldFiles, tol = 0.0, thresh = 0.5, method = "spacemap")
save.image(file = file.path(resPath, "her2-04.rda"))


cvsmap$cvVote$nfits
nonZeroUpper(cvsmap$cvVote$yy, 0.0)
#number of edges x-y
nonZeroWhole(cvsmap$cvVote$xy, 0.0)
#X-hub degree distribution
hist(rowSums(cvsmap$cvVote$xy)[rowSums(cvsmap$cvVote$xy) != 0])
cvsmap$minTune
summary(tmap)


######################
rm(list = ls())
resPath <- "~/scratch-data/neta-metabric/disc-cv-vote/basal/03/"
load(file.path(resPath, "basal-03.rda"))
source("~/repos/spacemap/R/cvVoteRDS.R")
minIndex <- cvsmap$minIndex
cvsmap$bestFoldFiles <- list.files(path = resPath, 
                                   pattern =  paste0("tuneid_", sprintf("%03d", minIndex), "_fold"),
                                   full.names = TRUE)
cvsmap$cvVote <- cvVoteRDS(files = cvsmap$bestFoldFiles, tol = 0.0, thresh = 0.5, method = "spacemap")
save.image(file = file.path(resPath, "basal-03.rda"))



nonZeroUpper(cvsmap$cvVote$yy, 0.0)
#number of edges x-y
nonZeroWhole(cvsmap$cvVote$xy, 0.0)
#X-hub degree distribution
hist(rowSums(cvsmap$cvVote$xy)[rowSums(cvsmap$cvVote$xy) != 0])
cvsmap$minTune
