
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
#PREVIOUS TUNING #
##################

prev <- new.env()
tunedir <- "~/scratch-data/neta-metabric/disc-cv-vote/lumb/04"
load(file = file.path(tunedir, "lumb-04-dropNA.rda"), envir = prev)
tune <- prev$cvsmap$minTune
tune

##################
#    ENSEMBLE    #
##################

#result directory
respath <- "/home/cconley/scratch-data/neta-metabric/disc-boot-vote/lumb/01/"

library(spacemap)
seed <- 971500
tictoc <- system.time({ens <- spacemap::bootEnsemble(Y = Y, X = X, tune = tune,
                                                     method = "spacemap", B = 250,
						     resPath = respath,
                                                     seed = seed, p0 = 0.90,
                                                     tol = 1e-4, cdmax = 90e7)})
save.image(file = file.path(respath, "lumb-boot-vote-p90-01.rda"))
#stop cluster to not run out of memory
stopCluster(cl)
#re-register the sequential backend
registerDoSEQ()
object.size(ens)
bv <- bootVote(ens)

##################
#  SAVE RESULTS  #
##################
saveRDS(object = bv, file = file.path(respath, "lumb-boot-vote-p90-01.rds"))
save.image(file = file.path(respath, "lumb-boot-vote-p90-01.rda"))
