
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


##################
#   TUNING       #
##################

#training/test setse
trainSets <- readRDS(file.path(filterdir, "discovery_train_sets_dropoutLumB.rds"))
testSets <- readRDS(file.path(filterdir, "discovery_test_sets_dropoutLumB.rds"))

##################
#   TRY 04       #
##################

#GRID
tmap <- expand.grid(lam1 = seq(115^2, 134^2, length = 10)^(1/2),
                    lam2 = seq(35^2, 48^2, length = 10)^(1/2),
                    lam3 = seq(3, 30, length = 5))

# plot(tmap$lam1, tmap$lam2, xlim  = c(110, 180), ylim  = c(20, 80))
# points(tmap2$lam1, tmap2$lam2, pch = 19)
# 
# plot(tmap$lam2, tmap$lam3, ylim  = c(0, 100), xlim  = c(20, 80))
# points(tmap2$lam2, tmap2$lam3, pch = 19)

#result directory
respath <- "/home/cconley/scratch-data/neta-metabric/disc-cv-vote/lumb/04/"
library(spacemap)
tictoc <- system.time({cvsmap <- cvVote(Y = Y, X = X,
                                        trainIds = trainSets, testIds = testSets,
                                        method = "spacemap", tuneGrid = tmap,
                                        resPath = respath,
                                        tol = 1e-4, cdmax = 40e7)})
save.image(file = file.path(respath, "lumb-04.rda"))
stopCluster(cl)
