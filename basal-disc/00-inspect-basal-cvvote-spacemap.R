

############
#   TRY 01 #
############

outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/basal/03"
load(file = file.path(outdir, "basal-03.rda"))

#selected lambda's
cvsmap$minTune

#score
cvsmap$logcvScore

library(spacemap)
#number of edges y-y
nonZeroUpper(cvsmap$cvVote$yy, 0.0)
#number of edges x-y
nonZeroWhole(cvsmap$cvVote$xy, 0.0)

#X-hub degree distribution
hist(rowSums(cvsmap$cvVote$xy)[rowSums(cvsmap$cvVote$xy) != 0])

## Inspect model fits under CV-selected tuning parameters
tmp <- sub(pattern = "/home/cconley/", replacement = "~/", cvsmap$bestFoldFiles)
init_sigs <- sapply(tmp, function(f) { 
  out <-readRDS(f)
  out$sig
})
sighat <- apply(X = init_sigs, MARGIN = 1, mean)
init <- readRDS(tmp[[1]])
rhohat <- as.matrix(init$ThetaYY)

#How much of a drop from individual fit to CV.Vote?
xy <- sapply(tmp, function(f) { 
  out <-readRDS(f)
  nonZeroWhole(as.matrix(out$ThetaXY), 1e-6)
})
yy <- sapply(tmp, function(f) { 
  out <-readRDS(f)
  nonZeroUpper(as.matrix(out$ThetaYY), 1e-6)
})
plot(xy, yy, pch = 19)


## Visualize
nsplits <- sapply(testSets, length)
cvsmap$metricScores$dfParCor <- cvsmap$metricScores$dfPqarCor
cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam1, tuneParam1Name = "lam1")
cvVis2 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam2, tuneParam1Name = "lam2")
cvVis3 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam3, tuneParam1Name = "lam3")

library(gridExtra)
grid.arrange(cvVis1[[1]], cvVis2[[1]],cvVis3[[1]], ncol = 2)
grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)
grid.arrange(cvVis2[[1]], cvVis2[[2]],cvVis2[[3]], cvVis2[[4]], ncol = 2)
grid.arrange(cvVis3[[1]], cvVis3[[2]],cvVis3[[3]], cvVis3[[4]], ncol = 2)

cvVis2[[4]] + coord_cartesian(ylim = c(0, 7000))

## Next Grid
tmap2 <- rbind(expand.grid(lam1 = c(67, 68, 69, 71, 72),
                           lam2 = seq(25.5, 27.5, length = 5), 
                           lam3 = seq(0, 20, length= 5)),
               expand.grid(lam1 =  c(67, 68, 69, 71, 72),
                           lam2 = seq(20, 23, length = 3), 
                           lam3 = seq(35, 80, length= 5)))



plot(tmap$lam1, tmap$lam2, xlim = c(67, 85), ylim = c(19, 33))
points(tmap2$lam1, tmap2$lam2, pch = 19)
plot(tmap$lam2, tmap$lam3, xlim = c(19, 33), ylim = c(0, 80))
points(tmap2$lam2, tmap2$lam3, pch = 19)



############
#   TRY 02 #
############

outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/basal/02"
load(file = file.path(outdir, "basal-02.rda"))

#selected lambda's
cvsmap$minTune

#score
cvsmap$logcvScore

library(spacemap)
#number of edges y-y
nonZeroUpper(cvsmap$cvVote$yy, 0.0)
#number of edges x-y
nonZeroWhole(cvsmap$cvVote$xy, 0.0)

#X-hub degree distribution
hist(rowSums(cvsmap$cvVote$xy)[rowSums(cvsmap$cvVote$xy) != 0])

## Visualize
nsplits <- sapply(testSets, length)
cvsmap$metricScores$dfParCor <- cvsmap$metricScores$dfPqarCor
cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam1, tuneParam1Name = "lam1")
cvVis2 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam2, tuneParam1Name = "lam2")
cvVis3 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam3, tuneParam1Name = "lam3")

library(gridExtra)
grid.arrange(cvVis1[[1]], cvVis2[[1]],cvVis3[[1]], ncol = 2)
grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)
grid.arrange(cvVis2[[1]], cvVis2[[2]],cvVis2[[3]], cvVis2[[4]], ncol = 2)
grid.arrange(cvVis3[[1]], cvVis3[[2]],cvVis3[[3]], cvVis3[[4]], ncol = 2)



library(ggplot2)
cvVis1[[1]] + coord_cartesian(ylim= c(7.18,7.20))
cvVis2[[1]] + coord_cartesian(ylim= c(7.18,7.20))
cvVis3[[1]] + coord_cartesian(ylim= c(7.15,7.20))

spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                  tuneParam1 = tmap$lam2, tuneParam1Name = "lam2",
                  tuneParam2 = tmap$lam3)

## Next Grid (Try 03)
idx <- which(rowMeans(cvsmap$metricScores$dfParCor) > 4500 &  rowMeans(cvsmap$metricScores$dfGamma) > 2000)

idx <- which(rowMeans(cvsmap$metricScores$dfParCor) > 4500 &  rowMeans(cvsmap$metricScores$dfGamma) > 6000)

tmap[idx,]
cvsmap$metricScores$dfParCor[idx,]
cvsmap$metricScores$dfGamma[idx,]
cvsmap$metricScores$dfGamma[cvsmap$minIndex,]

tmap2 <- expand.grid(lam1 = seq(63, 78, length = 10),
                     lam2 = seq(15, 28, length = 16), 
                     lam3 = seq(0, 80, length= 16))
tmap2 <- tmap2[!(tmap2$lam2 < 23 & tmap2$lam3 < 30),]
tmap2 <- tmap2[!(tmap2$lam2 > 22 & tmap2$lam3 > 30),]
plot(tmap$lam1, tmap$lam2, xlim = c(50, 85), ylim = c(10, 33))
points(tmap2$lam1, tmap2$lam2, pch = 19)
plot(tmap$lam2, tmap$lam3, xlim = c(10, 33), ylim = c(0, 80))
points(tmap2$lam2, tmap2$lam3, pch = 19)


######################################################
#basal 04



tmap4 <- expand.grid(lam1 = exp(seq(log(60), log(70), length = 10)),
                     lam2 = seq(20^4, 29^4, length = 10)^(1/4),
                     lam3 = exp(seq(log(3), log(25), length = 10)))


##################
# lam2 vs lam3   #
##################

library(RColorBrewer)
brewer.pal(name = "Dark2", n = 8)
outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/basal/01"
load(file = file.path(outdir, "basal-01.rda"))
plot(tmap$lam2, tmap$lam3, xlim = c(15, 35), ylim = c(0, 80), col = "#1B9E77", pch = 20)
cvsmap$logcvScore
points(cvsmap$minTune$lam2, cvsmap$minTune$lam3, pch = "*", col = "#1B9E77",cex = 3)

outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/basal/02"
load(file = file.path(outdir, "basal-02.rda"))
points(tmap$lam2, tmap$lam3, col = "#D95F02", pch = 20)
points(cvsmap$minTune$lam2, cvsmap$minTune$lam3, pch = "*", col = "#D95F02", cex = 3)
cvsmap$logcvScore

outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/basal/03"
load(file = file.path(outdir, "basal-03.rda"))
points(tmap$lam2, tmap$lam3, col = "#7570B3", pch = 20)
points(cvsmap$minTune$lam2, cvsmap$minTune$lam3, pch ="*", col = "#7570B3", cex = 3)
cvsmap$logcvScore

points(tmap4$lam2, tmap4$lam3, col = "#A6761D", pch = 19)


##################
# lam1 vs lam2   #
##################
library(RColorBrewer)
brewer.pal(name = "Dark2", n = 4)
outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/basal/01"
load(file = file.path(outdir, "basal-01.rda"))
plot(tmap$lam1, tmap$lam2, xlim = c(55, 85), ylim = c(15, 35), col = "#1B9E77", pch = 20)
cvsmap$logcvScore
points(cvsmap$minTune$lam1, cvsmap$minTune$lam2, pch = "*", col = "#1B9E77", cex = 3)

outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/basal/02"
load(file = file.path(outdir, "basal-02.rda"))
points(tmap$lam1, tmap$lam2, col = "#D95F02", pch = 20)
points(cvsmap$minTune$lam1, cvsmap$minTune$lam2, pch = "*", col = "#D95F02", cex = 3)
cvsmap$logcvScore
outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/basal/03"
load(file = file.path(outdir, "basal-03.rda"))
points(tmap$lam1, tmap$lam2, col = "#7570B3", pch = 20)
points(cvsmap$minTune$lam1, cvsmap$minTune$lam2, pch ="*", col = "#7570B3", cex = 3)
cvsmap$logcvScore

points(tmap4$lam1, tmap4$lam2, col = "#A6761D", pch = 19)

##################
# lam1 vs lam3   #
##################
library(RColorBrewer)
brewer.pal(name = "Dark2", n = 4)
outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/basal/01"
load(file = file.path(outdir, "basal-01.rda"))
plot(tmap$lam1, tmap$lam3, xlim = c(55, 85), ylim = c(0, 80), col = "#1B9E77", pch = 20)
cvsmap$logcvScore
points(cvsmap$minTune$lam1, cvsmap$minTune$lam3, pch = "*", col = "#1B9E77", cex = 3)

outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/basal/02"
load(file = file.path(outdir, "basal-02.rda"))
points(tmap$lam1, tmap$lam3, col = "#D95F02", pch = 20)
points(cvsmap$minTune$lam1, cvsmap$minTune$lam3, pch = "*", col = "#D95F02", cex = 3)
cvsmap$logcvScore
outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/basal/03"
load(file = file.path(outdir, "basal-03.rda"))
points(tmap$lam1, tmap$lam3, col = "#7570B3", pch = 20)
points(cvsmap$minTune$lam1, cvsmap$minTune$lam3, pch ="*", col = "#7570B3", cex = 3)
cvsmap$logcvScore


points(tmap4$lam1, tmap4$lam3, col = "#A6761D", pch = 19)

#next grid
idx <- which(rowMeans(cvsmap$metricScores$dfGamma) < 10000 & rowMeans(cvsmap$metricScores$dfParCor) < 10000)
tmlim <- tmap[idx,]




plot(tmap$lam2, tmap$lam3)
points(tmlim$lam2, tmlim$lam3, pch = 19)

plot(tmap$lam1, tmap$lam3)
points(tmlim$lam1, tmlim$lam3, pch = 19)


tmap <- expand.grid(lam1 = seq(55, 75, length = 10),
                    lam2 = seq(15, 35, length = 15),
                    lam3 = seq(0, 100, length = 15))


####################################

######################################################
#her2 03
#next grid


idx <- which(rowMeans(cvsmap$metricScores$dfParCor) > 4500 &  rowMeans(cvsmap$metricScores$dfGamma) > 2000)

tmap[idx,]
cvsmap$metricScores$dfParCor[idx,]
cvsmap$metricScores$dfGamma[idx,]
cvsmap$metricScores$dfGamma[cvsmap$minIndex,]

tmap2 <- expand.grid(lam1 = seq(50, 74, length = 10),
                     lam2 = seq(10, 26, length = 15), 
                     lam3 = seq(0, 65, length= 15))
tmap2 <- tmap2[!(tmap2$lam2 < 20 & tmap2$lam3 < 20),]

plot(tmap$lam1, tmap$lam2, xlim = c(50, 85), ylim = c(10, 33))
points(tmap2$lam1, tmap2$lam2, pch = 19)
plot(tmap$lam2, tmap$lam3, xlim = c(10, 33), ylim = c(0, 80))
points(tmap2$lam2, tmap2$lam3, pch = 19)


