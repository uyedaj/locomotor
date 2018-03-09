## Running analysis of locomotor data
library(bayou)
library(treeplyr)
setwd("~/repos/locomotor/R/")

tds <- readRDS("../output/tds.rds")

## Analysis of Maximum Speed
trait <- "Maximum Speed"
td <- tds[[trait]]
td <- treeply(td, multi2di, random=FALSE)
td <- filter(td, !is.na(LnMass))
tree <- td$phy
tree$edge.length[tree$edge.length==0] <- .Machine$double.eps
dat <- log(td[['TempTraitValue']])
pred <- select(td$dat, LnMass, Environment, TrophicGroup, Thermy)


chainNN <- list()
chainNN[[1]] <- readRDS('../output/LnMass/chain.modelNN_lnorm_r1099.rds')
chainNN[[2]] <- readRDS('../output/LnMass/chain.modelNN_lnorm_r4363.rds')
chainNN[[3]] <- readRDS('../output/LnMass/chain.modelNN_lnorm_r4499.rds')
chainNN[[4]] <- readRDS('../output/LnMass/chain.modelNN_lnorm_r5043.rds')
chainNN <- combine.chains(chainNN, burnin.prop=0.3)
plot(chainNN)

pdf("../output/branchheatmapsNN.pdf", width=10, height=20)
par(mfrow=c(1,2), mar=c(0,0,0,0))
plotBranchHeatMap(td$phy, chainNN, variable="theta", burnin=0, pal=viridis::viridis,edge.width=2, cex=0.5, x.lim=c(-20, 600))
par(mar=c(0,0,0,0))
plotBranchHeatMap(td$phy, chainNN, variable="beta_LnMass", burnin=0, pal=viridis::viridis,edge.width=2, cex=0.5, x.lim=c(-20, 600))
dev.off()



chainN1 <- list()
chainN1[[1]] <- readRDS('../output/LnMass/chain.modelN1_lnorm_r0672.rds')
chainN1[[2]] <- readRDS('../output/LnMass/chain.modelN1_lnorm_r0815.rds')
chainN1[[3]] <- readRDS('../output/LnMass/chain.modelN1_lnorm_r3174.rds')
chainN1[[4]] <- readRDS('../output/LnMass/chain.modelN1_lnorm_r5988.rds')
chainN1 <- combine.chains(chainN1, burnin.prop=0.3)
plot(chainN1)

pdf("../output/branchheatmapsN1.pdf", width=10, height=20)
par(mfrow=c(1,1), mar=c(0,0,0,0))
plotBranchHeatMap(td$phy, chainN1, variable="theta", burnin=0, pal=viridis::viridis,edge.width=2, cex=0.5, x.lim=c(-20, 600))
#par(mar=c(0,0,0,0))
#plotBranchHeatMap(td$phy, chainNN, variable="beta_LnMass", burnin=0, pal=viridis::viridis,edge.width=2, cex=0.5, x.lim=c(-20, 600))
dev.off()

mcmc.N1 <- readRDS("../output/LnMass/mcmc.modelN1_lnorm_r0672.rds")
mcmc.NN <- readRDS("../output/LnMass/mcmc.modelNN_lnorm_r1099.rds")

shiftsum.N1 <- shiftSummaries(chainN1, mcmc.N1, pp.cutoff=0.5)
shiftsum.NN <- shiftSummaries(chainNN, mcmc.NN, pp.cutoff=0.5)
pdf("../output/shiftsumN1.pdf")
plotShiftSummaries(shiftsum.N1, pal=viridis::viridis)
dev.off()

pdf("../output/shiftsumNN.pdf")
plotShiftSummaries(shiftsum.NN, pal=viridis::viridis)
dev.off()



## wH4

mcmcH4 <- list()
mcmcH4[[1]] <-   mcmc.wH4 <- bayou.makeMCMC(td$phy, dat, pred=.pred, SE=0.05, prior=prior.wH4, model=model.weightedH4, startpar=startpar,
                                            new.dir="../output/wH4/", outname="modelwH4_r88202", plot.freq=NULL, samp=100, ticker.freq=1000) # Set up the MCMC
mcmcH4[[2]] <-   mcmc.wH4 <- bayou.makeMCMC(td$phy, dat, pred=.pred, SE=0.05, prior=prior.wH4, model=model.weightedH4, startpar=startpar,
                                            new.dir="../output/wH4/", outname="modelwH4_r73077", plot.freq=NULL, samp=100, ticker.freq=1000) # Set up the MCMC
mcmcH4[[3]] <-   mcmc.wH4 <- bayou.makeMCMC(td$phy, dat, pred=.pred, SE=0.05, prior=prior.wH4, model=model.weightedH4, startpar=startpar,
                                            new.dir="../output/wH4/", outname="modelwH4_r63235", plot.freq=NULL, samp=100, ticker.freq=1000) # Set up the MCMC
mcmcH4[[4]] <-   mcmc.wH4 <- bayou.makeMCMC(td$phy, dat, pred=.pred, SE=0.05, prior=prior.wH4, model=model.weightedH4, startpar=startpar,
                                            new.dir="../output/wH4/", outname="modelwH4_r94537", plot.freq=NULL, samp=100, ticker.freq=1000) # Set up the MCMC
mcmcH4[[5]] <-   mcmc.wH4 <- bayou.makeMCMC(td$phy, dat, pred=.pred, SE=0.05, prior=prior.wH4, model=model.weightedH4, startpar=startpar,
                                            new.dir="../output/wH4/", outname="modelwH4_r101194", plot.freq=NULL, samp=100, ticker.freq=1000) # Set up the MCMC
mcmcH4[[6]] <-   mcmc.wH4 <- bayou.makeMCMC(td$phy, dat, pred=.pred, SE=0.05, prior=prior.wH4, model=model.weightedH4, startpar=startpar,
                                            new.dir="../output/wH4/", outname="modelwH4_r116622", plot.freq=NULL, samp=100, ticker.freq=1000) # Set up the MCMC
mcmcH4[[7]] <-   mcmc.wH4 <- bayou.makeMCMC(td$phy, dat, pred=.pred, SE=0.05, prior=prior.wH4, model=model.weightedH4, startpar=startpar,
                                            new.dir="../output/wH4/", outname="modelwH4_r122159", plot.freq=NULL, samp=100, ticker.freq=1000) # Set up the MCMC
mcmcH4[[8]] <-   mcmc.wH4 <- bayou.makeMCMC(td$phy, dat, pred=.pred, SE=0.05, prior=prior.wH4, model=model.weightedH4, startpar=startpar,
                                            new.dir="../output/wH4/", outname="modelwH4_r51493", plot.freq=NULL, samp=100, ticker.freq=1000) # Set up the MCMC


chainH4 <- lapply(mcmcH4, function(x) x$load())
chain.wH4 <- combine.chains(chainH4, burnin.prop = 0.3, thin=10)
plot(chain.wH4)



