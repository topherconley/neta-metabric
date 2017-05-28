
####################
# PARALLEL BACKEND #
####################


#Take in arguments from sbatch. 
args <- commandArgs(TRUE)

if (length(args) == 0) { 
  ncores <- 1
  comp <- "singlenode"
} else {
  ncores <- as.integer(args[1])
  comp <- "gauss" 
}

if (comp == "gauss") {
  
  #====Setup for running on Gauss... ======#
  
  getnodelist <- function(maxpernode=30, f = "nodelist.txt") {
    nodelist <- Sys.getenv('SLURM_NODELIST')
    system(paste0("scontrol show hostname ", nodelist, " > ", f))
    x <- if (nzchar(nodelist)) readLines(f) else rep('localhost', 3)
    d <- as.data.frame(table(x), stringsAsFactors=FALSE)
    rep(d$x, pmax(d$Freq, maxpernode))
  }
  
  library(doParallel)
  nodelist <- getnodelist(maxpernode=ncores, f = "nodelist_01.txt")
  print(nodelist)
  cl <- makePSOCKcluster(nodelist, outfile='')
  registerDoParallel(cl)
} else if (comp == "singlenode") {
  #===== Setup for running on a multicore machine =====#
  suppressPackageStartupMessages(library(doParallel))
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
} else {
  stop("Argument: comp is not specified correctly.")
}


##################
#     Data       #
##################
filterdir <- "/home/cconley/scratch-data/neta-metabric/disc-filtered/"
fy <- "disc-mrna-eset-nostd-union-dropout-std-Basal.rds"
fx <- "disc-cna-eset-nostd-union-dropout-multi-std-Basal.rds"
library(Biobase)
xset <- readRDS(file.path(filterdir, fx))
X <- t(exprs(xset))
yset <- readRDS(file.path(filterdir, fy))
Y <- t(exprs(yset))
#assert the sample names match
stopifnot(all(sampleNames(xset) == sampleNames(yset)))

##################
#   INIT FIT     #
##################

#DEFINE lam1 from previous efforts 
#load("~/scratch-data/old-neta-metabric/disc-cv-vote/space_yy/Basal/04/workspace.RData")
#cvSpaceYY$tuneGrid
# 73.22034
#define grid
# 
# tmap1 <- expand.grid(lam1 = 80, 
#                      lam2 = c(20, 30, 35), 
#                      lam3 = 20)
# tmap1 <- expand.grid(lam1 = 80, 
#                      lam2 = c(25,27, 29), 
#                      lam3 = 20)
# library(spacemap)
# ntmap1 <- initFit(Y = Y, X = X, method = "spacemap", tuneGrid = tmap1,
#                   iscale = FALSE, tol = 1e-4, cdmax = 6e7)
# 
# 
# library(ggplot2)
# qplot(x = tmap1$lam2, y = ntmap1$nxy,
#       ylab = "No. of X->Y edges", xlab = "lam2") +
#   theme_bw()


##################
#   TUNING       #
##################

#training/test setse
trainSets <- readRDS(file.path(filterdir, "discovery_train_sets_dropoutBasal.rds"))
testSets <- readRDS(file.path(filterdir, "discovery_test_sets_dropoutBasal.rds"))

##################
#   TRY O1       #
##################

#GRID
tmap <- expand.grid(lam1 = seq(70, 85, length = 5),
                    lam2 = seq(25, 33, length = 5),
                    lam3 = seq(10, 70, length = 5))
#result directory
respath <- "/home/cconley/scratch-data/neta-metabric/disc-cv-vote/basal/01/"
library(spacemap)
tictoc <- system.time({cvsmap <- cvVote(Y = Y, X = X,
                                        trainIds = trainSets, testIds = testSets,
                                        method = "spacemap", tuneGrid = tmap,
                                        resPath = respath,
                                       tol = 1e-4, cdmax = 6e7)})
save.image(file = file.path(respath, "basal-01.rda"))

stopCluster(cl)
