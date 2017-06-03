
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
filterdir <- "~/scratch-data/neta-metabric/disc-filtered/"
fy <- "disc-mrna-eset-nostd-union-dropout-std-Her2.rds"
fx <- "disc-cna-eset-nostd-union-dropout-multi-std-Her2.rds"
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
trainSets <- readRDS(file.path(filterdir, "discovery_train_sets_dropoutHer2.rds"))
testSets <- readRDS(file.path(filterdir, "discovery_test_sets_dropoutHer2.rds"))

##################
#   TRY 05       #
##################

tmap <- expand.grid(lam1 = exp(seq(log(50), log(70), length = 15)),
                    lam2 = seq(10, 20, length = 10), 
                    lam3 = c(sqrt(seq(2^2, 10^2, length= 5)),
                             exp(seq(log(11), log(20), length= 5))))

#result directory
respath <- "/home/cconley/scratch-data/neta-metabric/disc-cv-vote/her2/05/"
library(spacemap)
tictoc <- system.time({cvsmap <- spacemap::cvVote(Y = Y, X = X,
                                                  trainIds = trainSets, testIds = testSets, 
                                                  method = "spacemap", tuneGrid = tmap, 
                                                  resPath = respath,
                                                  tol = 1e-4, cdmax = 30e7)})
save.image(file = file.path(respath, "her2-05.rda"))
stopCluster(cl)
