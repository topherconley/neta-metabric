
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
#   TRY O4       #
##################


tmap <- expand.grid(lam1 = seq(50, 74, length = 10),
                    lam2 = seq(10, 26, length = 15), 
                    lam3 = seq(0, 65, length= 15))
tmap <- tmap[!(tmap$lam2 < 20 & tmap$lam3 < 20),]

#result directory
respath <- "/home/cconley/scratch-data/neta-metabric/disc-cv-vote/her2/04/"
library(spacemap)
tictoc <- system.time({cvsmap <- spacemap::cvVote(Y = Y, X = X,
                                                  trainIds = trainSets, testIds = testSets, 
                                                  method = "spacemap", tuneGrid = tmap, 
                                                  resPath = respath,
                                                  tol = 1e-4, cdmax = 30e7)})
save.image(file = file.path(respath, "her2-04.rda"))
stopCluster(cl)
