outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/her2/05"
load(file = file.path(outdir, "her2-05.rda"))


nsplits <- sapply(testSets, length)
cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                            tuneParam1 = tmap$lam1, tuneParam1Name = "lam1")
library(gridExtra)
grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)

library(spacemap)
nonZeroUpper(cvsmap$cvVote$yy,0.0)
nonZeroWhole(cvsmap$cvVote$xy,0.0)

bv <- readRDS(file.path("~/scratch-data/neta-metabric/disc-cv-vote/her2/06/","her2-06-boot-vote.rds"))
nonZeroUpper(bv$bv$yy,0.0)
nonZeroWhole(bv$bv$xy,0.0)


hist(rowSums(bv$bv$xy))
str(bv$bdeg$xy)

str(bv$bc$xy)



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



#################################################
#from cv 
###################################################

head(avg)

min(log(avg[,"rss"]))


#################################################
#from tuneVis
###################################################

testSetLen <- sapply(testSets, length)
cvOut <- cvsmap
total <- sum(testSetLen)
foldconvmat <- !is.na(cvOut$metricScores[["rss"]])
nfoldconv <- rowSums(foldconvmat)
ntestconv <- apply(X = foldconvmat, MARGIN  = 1, FUN = function(not_na) sum(testSetLen[not_na]))
requireNamespace("foreach")
#for R CMD check NOTE passing
i <- NULL;
avgScores <- foreach(i = seq_along(cvOut$metricScores), .combine = 'cbind') %do% {
  score <- cvOut$metricScores[[i]]
  stype <- names(cvOut$metricScores)[[i]]
  msum <- apply(X = score , MARGIN = 1, FUN = sum, na.rm = TRUE) 
  
  ###Selection of tuning grid restricted to case when majority of test hold outs
  ###were evaluated (chunked by folds) 
  
  # divide rss by number of test sets with corresponding trained convergence
  #divide other metrics by the number of converged folds
  
  if (stype == "rss") { 
    ifelse(ntestconv > 0.5*total, msum / ntestconv, Inf)
  } else {
    ifelse(ntestconv > 0.5*total, msum / nfoldconv, Inf)
  }
}
colnames(avgScores) <- names(cvOut$metricScores)

head(avgScores)
head(avg)


grid.arrange(cvVis1[[1]] + ,
             , cvVis2[[1]],cvVis3[[1]], ncol = 2)
grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)
grid.arrange(cvVis2[[1]], cvVis2[[2]],cvVis2[[3]], cvVis2[[4]], ncol = 2)
grid.arrange(cvVis3[[1]], cvVis3[[2]],cvVis3[[3]], cvVis3[[4]], ncol = 2)