outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/her2/04"
load(file = file.path(outdir, "her2-04.rda"))


nsplits <- sapply(testSets, length)
cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam1, tuneParam1Name = "lam1")
library(gridExtra)
grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)

library(spacemap)
nonZeroUpper(cvsmap$cvVote$yy,0.0)
nonZeroWhole(cvsmap$cvVote$xy,0.0)




outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/her2/03"
load(file = file.path(outdir, "her2-03.rda"))


nsplits <- sapply(testSets, length)
cvsmap$metricScores$dfParCor <- cvsmap$metricScores$dfPqarCor
cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam1, tuneParam1Name = "lam1")
cvVis2 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam2, tuneParam1Name = "lam2")
cvVis3 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam3, tuneParam1Name = "lam3")

