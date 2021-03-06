
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
#   TUNING       #
##################

#training/test setse
trainSets <- readRDS(file.path(filterdir, "discovery_train_sets_dropoutBasal.rds"))
testSets <- readRDS(file.path(filterdir, "discovery_test_sets_dropoutBasal.rds"))


##################
#   TRY 02       #
##################

#GRID
tmap <- rbind(expand.grid(lam1 = c(67, 68, 69, 71, 72),
                          lam2 = seq(25.5, 27.5, length = 5), 
                          lam3 = seq(0, 20, length= 5)),
              expand.grid(lam1 =  c(67, 68, 69, 71, 72),
                          lam2 = seq(20, 23, length = 3), 
                          lam3 = seq(35, 80, length= 5)))
#do most intensive computations first
tmap <- tmap[with(data = tmap, order(lam1, lam2, lam3)),]
#result directory
respath <- "/home/cconley/scratch-data/neta-metabric/disc-cv-vote/basal/02/"

library(spacemap)
tictoc <- system.time({cvsmap <- cvVote(Y = Y, X = X, 
                                        trainIds = trainSets, testIds = testSets, 
                                        method = "spacemap", tuneGrid = tmap, 
                                        resPath = respath,
                                        tol = 1e-4, cdmax = 6e7)})

save.image(file = file.path(respath, "basal-02.rda"))
stopCluster(cl)
