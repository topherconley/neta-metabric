
####################
# PARALLEL BACKEND #
####################


#Take in arguments from sbatch. 
args <- commandArgs(TRUE)

if (length(args) == 0) { 
  ncores <- 3
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
fy <- "disc-mrna-eset-nostd-union-dropout-std-LumB.rds"
fx <- "disc-cna-eset-nostd-union-dropout-multi-std-LumB.rds"
library(Biobase)
xset <- readRDS(file.path(filterdir, fx))
X <- t(exprs(xset))
yset <- readRDS(file.path(filterdir, fy))
Y <- t(exprs(yset))
#assert the sample names match
stopifnot(all(sampleNames(xset) == sampleNames(yset)))


#training/test setse
trainSets <- readRDS(file.path(filterdir, "discovery_train_sets_dropoutLumB.rds"))
testSets <- readRDS(file.path(filterdir, "discovery_test_sets_dropoutLumB.rds"))

Ytrain <- Y[trainSets$Fold01,]
Xtrain <- X[trainSets$Fold01,]
##################
#   INIT FIT     #
##################

#DEFINE lam1 from previous efforts
# load("~/scratch-data/old-neta-metabric/disc-cv-vote/space_yy/LumB/03/workspace.RData")
# cvSpaceYY$tuneGrid
# 142.8571

tmap <- expand.grid(lam1 = seq(125, 160, length = 5),
                    lam2 = seq(75, 100, length = 5), 
                    lam3 = seq(0, 100, length = 5))

library(spacemap)
tictoc <- system.time({ntmap1 <- initFit(Y = Ytrain, X = Xtrain, method = "spacemap", tuneGrid = tmap,
                                         iscale = FALSE, tol = 1e-3, cdmax = 10e7)})
respath <- "/home/cconley/scratch-data/neta-metabric/disc-cv-vote/lumb/01/"
save.image(file = file.path(respath, "lumb-01.rda"))
stopCluster(cl)



