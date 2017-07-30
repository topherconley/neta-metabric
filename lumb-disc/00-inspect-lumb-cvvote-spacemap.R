
outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/lumb/01/"
load(file.path(outdir, "lumb-01.rda"))


outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/lumb/01/"
load(file.path(outdir, "lumb-01a.rda"))


plot(tmap$lam2, ntmap1$nxy)
plot(tmap$lam3, ntmap1$nxy)
plot(tmap$lam1, ntmap1$nyy)


outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/lumb/02/"
load(file.path(outdir, "lumb-02-dropNA.rda"))

outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/lumb/03/"
load(file.path(outdir, "lumb-03-dropNA.rda"))


outdir <- "~/scratch-data/neta-metabric/disc-cv-vote/lumb/04/"
load(file.path(outdir, "lumb-04-dropNA.rda"))


cvsmap$logcvScore

metricScores <- cvsmap$metricScores
testSetLen <- nsplits

rss <- log(ifelse(ntestconv > 0.5*total, msum / ntestconv, Inf))

sorted <- sort(rss, index.return  = T, decreasing = F)
tmap[head(sorted$ix),]
plot(rss)


i02$cvsmap$minTune
i03$cvsmap$minTune
i04$cvsmap$minTune


i02$cvsmap$logcvScore
i03$cvsmap$logcvScore
i04$cvsmap$logcvScore


