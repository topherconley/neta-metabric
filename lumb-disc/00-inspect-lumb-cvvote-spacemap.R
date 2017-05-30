
outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/lumb/01/"
load(file.path(outdir, "lumb-01.rda"))


outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/lumb/01/"
load(file.path(outdir, "lumb-01a.rda"))


plot(tmap$lam2, ntmap1$nxy)
plot(tmap$lam3, ntmap1$nxy)
plot(tmap$lam1, ntmap1$nyy)
