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

prior.11 <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dlnorm", dsig2="dhalfcauchy", dbeta_LnMass="dnorm",
                                  dsb="fixed", dk="fixed", dtheta="dnorm"), 
                       param=list(dalpha=list(meanlog=-3.5, sdlog=2.5), dsig2=list(scale=0.1),
                                  dbeta_LnMass=list(mean=0.1, sd=0.25),
                                  dtheta=list(mean=0, sd=2.5)),
                       fixed=list(k=0, sb=numeric(0))
)

prior.N1 <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dlnorm", dsig2="dhalfcauchy", dbeta_LnMass="dnorm",
                                  dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                       param=list(dalpha=list(meanlog=-3.5, sdlog=2.5), dsig2=list(scale=0.1),
                                  dbeta_LnMass=list(mean=0.1, sd=0.25),
                                  dk=list(lambda=17.35, kmax=69.4),
                                  dtheta=list(mean=0, sd=2.5))
)


prior.NN <- make.prior(tree, plot.prior = FALSE, 
                       dists=list(dalpha="dlnorm", dsig2="dhalfcauchy", dbeta_LnMass="dnorm",
                                  dsb="dsb", dk="cdpois", dtheta="dnorm"), 
                       param=list(dalpha=list(meanlog=-3.5, sdlog=2.5), dsig2=list(scale=0.1),
                                  dbeta_LnMass=list(mean=0.1, sd=0.25),
                                  dk=list(lambda=17.35, kmax=69.4), 
                                  dtheta=list(mean=0, sd=2.5))
)

D11 = list(alpha=1, sig2=1, beta_LnMass=0.1, k=1, theta=0.5, slide=1)
DN1 = list(alpha=1, sig2=1, beta_LnMass=0.1, k=1, theta=2, slide=1)
DNN = list(alpha=1, sig2=1, beta_LnMass=0.6, k=c(1,1), theta=3, slide=1)

MEvar <- 0.05
gens <- 1000000
## Make MCMC objects:

runmodel <- function(i, gens=1000000){
  if(i==1){
  model.11 <- makeBayouModel(dat ~ LnMass, rjpars = c(), 
                               tree=tree, dat=dat, pred=pred, SE=MEvar, prior=prior.11, D=D11, slopechange="alphaweighted")
  name.11 <- paste0("model11_lnorm_r", paste0(sample(0:9,4,replace=TRUE), collapse=""))
  mcmc.11 <- bayou.makeMCMC(tree, dat, pred=pred, SE=MEvar, model=model.11$model, prior=prior.11, samp = 100, startpar=model.11$startpar, new.dir="../output/LnMass/", outname=name.11, plot.freq=NULL)
  mcmc.11$run(gens)
  chain.11 <- set.burnin(mcmc.11$load(), 0.3)
  saveRDS(chain.11, file=paste("../output/LnMass/chain.", name.11, ".rds", sep=""))
  saveRDS(mcmc.11, file=paste("../output/LnMass/mcmc.", name.11, ".rds", sep=""))
  return(list(mcmc.11, chain.11, name.11))
  }

  #endpar.11 <<- pull.pars(length(chain.11$gen), chain.11, model=model.11$model)
  #startpar.N1 <<- model.N1$startpar
  #startpar.N1$theta[1] <- endpar.11$theta
  #startpar.N1$beta_LnMass <- endpar.11$beta_LnMass
  if(i==2){
  model.N1 <- makeBayouModel(dat ~ LnMass, rjpars = c("theta"),  
                               tree=tree, dat=dat, pred=pred, SE=MEvar, prior=prior.N1, D=DN1, slopechange="alphaweighted")
  name.N1 <- paste0("modelN1_lnorm_r", paste0(sample(0:9,4,replace=TRUE), collapse=""))
  mcmc.N1 <- bayou.makeMCMC(tree, dat, pred=pred, SE=MEvar, model=model.N1$model, prior=prior.N1, samp = 100, startpar=model.N1$startpar, new.dir="../output/LnMass/", outname=name.N1, plot.freq=NULL)
  #mcmc.N1 <- readRDS("../output/LnMass/mcmc.N1.r001.rds")
  mcmc.N1$run(gens)
  chain.N1 <- set.burnin(mcmc.N1$load(), 0.3)
  saveRDS(chain.N1, file=paste("../output/LnMass/chain.", name.N1, ".rds", sep=""))
  saveRDS(mcmc.N1, file=paste("../output/LnMass/mcmc.", name.N1, ".rds", sep=""))
  return(list(mcmc.N1, chain.N1, name.N1))
  }

#endpar.N1 <- pull.pars(length(chain.N1$gen), chain.N1, model=model.N1$model)
#startpar.NN <- endpar.N1
#startpar.NN$beta_LnMass <- rep(startpar.NN$beta_LnMass, startpar.NN$ntheta)
  if(i==3){
  model.NN <- makeBayouModel(dat ~ LnMass, rjpars = c("theta", "LnMass"),  
                               tree=tree, dat=dat, pred=pred, SE=MEvar, prior=prior.NN, D=DNN, slopechange="alphaweighted")
  name.NN <- paste0("modelNN_lnorm_r", paste0(sample(0:9,4,replace=TRUE), collapse=""))
  mcmc.NN <- bayou.makeMCMC(tree, dat, pred=pred, SE=MEvar, model=model.NN$model, prior=prior.NN, samp = 100, startpar=model.NN$startpar, new.dir="../output/LnMass/", outname=name.NN, plot.freq=NULL)
  mcmc.NN$run(gens)
  chain.NN <- set.burnin(mcmc.NN$load(), 0.3)
  saveRDS(chain.NN, file=paste("../output/LnMass/chain.", name.NN, ".rds", sep=""))
  saveRDS(mcmc.NN, file=paste("../output/LnMass/mcmc.", name.NN, ".rds", sep=""))
  return(list(mcmc.NN, chain.NN, name.NN))
  }
}

library(foreach)
library(doParallel)
ncores <- 50
nsteps <- 50
registerDoParallel(cores=ncores)
res <- foreach(i=c(1,2,2,2,2,3,3,3,3)) %dopar% runmodel(i, gens=1000000)

mcmc.11 <- res[[1]][[1]]
mcmc.N1 <- res[[2]][[1]]
mcmc.NN <- res[[6]][[1]]

chain.11 <- res[[1]][[2]]
chain.N1 <- combine.chains(lapply(2:5, function(x) res[[x]][[2]]), burnin.prop=0.3)
chain.NN <- combine.chains(lapply(6:9, function(x) res[[x]][[2]]), burnin.prop=0.3)

ngens <- 500000
Bk <- qbeta(seq(0,1, length.out=nsteps), 0.3,1)

ss.11 <- mcmc.11$steppingstone(ngens, chain.11, Bk, burnin=0.3, plot=FALSE)
saveRDS(ss.11, file="../output/LnMass/ss.11.r005.rds")

ss.N1 <- mcmc.N1$steppingstone(ngens, chain.N1, Bk, burnin=0, plot=FALSE)
saveRDS(ss.N1, file="../output/LnMass/ss.N1.r005.rds")

ss.NN <- mcmc.NN$steppingstone(ngens, chain.NN, Bk, burnin=0, plot=FALSE)
saveRDS(ss.NN, file="../output/LnMass/ss.NN.r005.rds")

ss.11
ss.N1
ss.NN


shiftsum.N1 <- shiftSummaries(chain.N1, mcmc.N1, pp.cutoff=0.5)
shiftsum.NN <- shiftSummaries(chain.NN, mcmc.NN, pp.cutoff=0.5)
pdf("../output/shiftsumN1.pdf")
plotShiftSummaries(shiftsum.N1, pal=viridis::viridis)
dev.off()

pdf("../output/shiftsumNN.pdf")
plotShiftSummaries(shiftsum.NN, pal=viridis::viridis, single.plot=TRUE)
dev.off()
