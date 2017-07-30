averageScores <- function(metricScores, testIds) { 
  ##Average cross-validated scores across folds
  #while accounting for non-convergent folds
  testSetLen <- sapply(testIds, length)
  total <- sum(testSetLen)
  foldconvmat <- !is.na(metricScores$rss)
  nfoldconv <- rowSums(foldconvmat)
  ntestconv <- apply(X = foldconvmat, MARGIN  = 1, FUN = function(not_na) sum(testSetLen[not_na]))
  i <- 0L #for R CMD check 
  library(foreach)
  metricScoresAvg <- foreach(i = seq_along(metricScores), .combine = 'cbind') %do% {
    score <- metricScores[[i]]
    stype <- names(metricScores)[[i]]
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
  colnames(metricScoresAvg) <- names(metricScores)
  metricScoresAvg
}


dxsmap <- function(file, lims = T, .probs = c(0, 0.25)) { 
  load(file)
  library(spacemap)
  
  #number of edges y-y
  nyy <- nonZeroUpper(cvsmap$cvVote$yy, 0.0)
  #number of edges x-y
  nxy <- nonZeroWhole(cvsmap$cvVote$xy, 0.0)
  
  #############
  library(ggplot2)
  #X-hub degree distribution
  xhubdeg <- rowSums(cvsmap$cvVote$xy)[rowSums(cvsmap$cvVote$xy) != 0]
  gxhub <- qplot(xhubdeg, xlab = "X-hub Degree") + 
    geom_histogram(binwidth = 1) + theme_bw()
  
  
  ## Inspect model fits under CV-selected tuning parameters
  tmp <- sub(pattern = "/home/cconley/", replacement = "~/", cvsmap$bestFoldFiles)
  #How much of a drop from individual fit to CV.Vote?
  xy <- sapply(tmp, function(f) { 
    out <-readRDS(f)
    nonZeroWhole(as.matrix(out$ThetaXY), 1e-6)
  })
  yy <- sapply(tmp, function(f) { 
    out <-readRDS(f)
    nonZeroUpper(as.matrix(out$ThetaYY), 1e-6)
  })
  cvdrop <- data.frame(xy = xy, yy = yy, type = "fold")
  cvdrop <- rbind(data.frame(xy = nxy, yy = nyy, type = "cv.vote"),
                  cvdrop)
  gcvdrop <- ggplot(data = cvdrop, mapping = aes(x = xy, y = yy, shape = type)) + 
    geom_point() + theme_bw()
  

  ## CV curve and model size
  #fix bug
  if(!is.null(cvsmap$metricScores$dfPqarCor)) { 
    cvsmap$metricScores$dfParCor <- cvsmap$metricScores$dfPqarCor
    cvsmap$metricScores$dfPqarCor <- NULL
  }
  nsplits <- sapply(testSets, length)
  cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                              tuneParam1 = tmap$lam1, tuneParam1Name = "lam1")
  cvVis2 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                              tuneParam1 = tmap$lam2, tuneParam1Name = "lam2")
  cvVis3 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                              tuneParam1 = tmap$lam3, tuneParam1Name = "lam3")
  
  avg <- averageScores(cvsmap$metricScores, testSets)
  library(gridExtra)
  if (lims) {
    lims1 <- coord_cartesian(ylim = quantile(log(avg[,"rss"]), probs = .probs))
    grid.arrange(cvVis1[[1]] + lims1, 
                 cvVis2[[1]] + lims1,
                 cvVis3[[1]] + lims1, ncol = 2)
  } else { 
    grid.arrange(cvVis1[[1]], cvVis2[[1]], cvVis3[[1]], ncol = 2)
  }

  
  #grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)
  #grid.arrange(cvVis2[[1]], cvVis2[[2]],cvVis2[[3]], cvVis2[[4]], ncol = 2)
  #grid.arrange(cvVis3[[1]], cvVis3[[2]],cvVis3[[3]], cvVis3[[4]], ncol = 2)
  list(nxy = nxy, nyy = nyy, gxhub = gxhub, gcvdrop = gcvdrop,
       v1 = cvVis1, v2 = cvVis2, v3 = cvVis3, avg = as.data.frame(avg))
} 


dxspace <- function(file, lims = T, .probs = c(0, 0.25)) { 
  load(file)
  library(spacemap)

  #number of edges y-y
  nyy <- nonZeroUpper(cvsmap$cvVote$yy, 0.0)

  #############
  library(ggplot2)

  ## Inspect model fits under CV-selected tuning parameters
  tmp <- sub(pattern = "/home/cconley/", replacement = "~/", cvsmap$bestFoldFiles)
  yy <- sapply(tmp, function(f) { 
    out <-readRDS(f)
    nonZeroUpper(as.matrix(out$ThetaYY), 1e-6)
  })
  cvdrop <- data.frame(yy = yy, type = "fold")
  gcvdrop <- ggplot(data = cvdrop, mapping = aes(x = yy)) + 
    geom_histogram() + geom_vline(xintercept = nyy) + theme_bw()
  
  ## CV curve and model size
  nsplits <- sapply(testSets, length)
  #fix bug
  if(!is.null(cvsmap$metricScores$dfPqarCor)) { 
    cvsmap$metricScores$dfParCor <- cvsmap$metricScores$dfPqarCor
    cvsmap$metricScores$dfPqarCor <- NULL
  }
  cvVis1 <- spacemap::tuneVis(cvOut = cvsmap, testSetLen = nsplits, 
                              tuneParam1 = tmap$lam1, tuneParam1Name = "lam1")
  
  avg <- averageScores(cvsmap$metricScores, testSets)
  library(gridExtra)
  if (lims) {
    lims1 <- coord_cartesian(ylim = quantile(log(avg[,"rss"]), probs = .probs))
    grid.arrange(cvVis1[[1]] + lims1, 
                 cvVis1[[2]], cvVis1[[3]],
                 cvVis1[[4]], ncol = 2)
  } else { 
    grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)
  }
  
  #grid.arrange(cvVis1[[1]], cvVis1[[2]],cvVis1[[3]], cvVis1[[4]], ncol = 2)
  #grid.arrange(cvVis2[[1]], cvVis2[[2]],cvVis2[[3]], cvVis2[[4]], ncol = 2)
  #grid.arrange(cvVis3[[1]], cvVis3[[2]],cvVis3[[3]], cvVis3[[4]], ncol = 2)
  list(nyy = nyy, gcvdrop = gcvdrop,
       v1 = cvVis1, avg = as.data.frame(avg))
} 






