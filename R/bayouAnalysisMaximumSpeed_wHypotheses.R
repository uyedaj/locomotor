## Running analysis of locomotor data
library(bayou)
library(treeplyr)
library(diversitree)
library(DirichletReg)
.dirichletMove <- function(cache=NULL, pars, d=NULL, move, ct=NULL, prior=NULL){
  prop <- MCMCpack::rdirichlet(1, pars[[move]]*d)
  if(any(prop==0)){
    lnHastingsRatio <- -Inf
    prop[prop==0] <- .Machine$double.eps
    prop <- prop/sum(prop)
  } else {
    lnHastingsRatio <- DirichletReg::ddirichlet(matrix(pars[[move]], nrow=1), prop, log=TRUE)  - DirichletReg::ddirichlet(matrix(prop, nrow=1), pars[[move]], log=TRUE)
  }
  pars.new <- pars
  pars.new[[move]] <- prop
  return(list(pars=pars.new, hr=lnHastingsRatio))
}

ddirichlet <- function(x, alpha, log=TRUE, sum.up=FALSE){
  DirichletReg::ddirichlet(matrix(x, nrow=1), alpha, log=log, sum.up=sum.up)
}
prop.W <- NULL
W <- NULL

setwd("~/repos/locomotor/R/")

tds <- readRDS("../output/tds.rds")

## Analysis of Maximum Speed
trait <- "Maximum Speed"
td <- tds[[trait]]
td <- treeply(td, multi2di, random=FALSE)
td <- filter(td, !is.na(LnMass))
td$phy$edge.length[td$phy$edge.length==0] <- .Machine$double.eps
nH <- nodeHeights(td$phy)
externalEdge <- which(td$phy$edge[,2] <= length(td$phy$tip.label))
td$phy$edge.length[externalEdge] <- td$phy$edge.length[externalEdge] + (max(nH) - nH[externalEdge,2])
td <- reorder(td, "postorder")
is.ultrametric(td$phy)

dat <- log(td[['TempTraitValue']])
pred <- select(td$dat, LnMass, Environment, TrophicGroup, Thermy)
pred$TrophicGroup <- recode_factor(pred$TrophicGroup, "Carnivore"="carnivore", "Omnivore"="omnivore", "Producer"="producer")
pred$TrophicGroup <- factor(pred$TrophicGroup, levels=c("herbivore", "omnivore","carnivore"))
pred <- mutate(pred, "fly"=as.numeric(Environment=="Air"), "swim"=as.numeric(Environment=="Water"))
plot(pred$LnMass, dat, pch=21, bg=pred$TrophicGroup)


## Predictors Environment, Trophic Group & Thermy; Generate evolutionary history

## Generate hypothesized history of environment on phylogeny
plot(td$phy, type="fan", show.tip.label=FALSE)
tiplabels(pch=21, bg=pred$Environment)
tipEnv <- as.numeric(factor(pred$Environment, levels=c("Water", "Land", "Air")))
likEnvironment <- make.mkn(td$phy, setNames(as.numeric(factor(pred$Environment, levels=c("Water", "Land", "Air"))), td$phy$tip.label),3)
likEnvironment <- constrain(likEnvironment, q13~0, q31~0)
mleEnvironment <- find.mle(likEnvironment, c(0.01, 0.01, 0.01, 0.01), method="optim")
envJoint <- round(apply(asr.joint(likEnvironment, mleEnvironment$par, n=100),2,mean),0)
plot(td$phy, show.tip.label=FALSE, type="fan")
nodelabels(pch=21, bg=envJoint)
tiplabels(pch=21, bg=as.numeric(factor(pred$Environment, levels=c("Water", "Land", "Air"))), col=as.numeric(factor(pred$Environment, levels=c("Water", "Land", "Air"))))
statesEnv <- c(tipEnv, envJoint)
H_env <- list(sb=which(statesEnv[td$phy$edge[,1]] != statesEnv[td$phy$edge[,2]]))
H_env$k <- length(H_env$sb)
H_env$ntheta <- 3
H_env$t2 <- statesEnv[td$phy$edge[H_env$sb,2]]
H_env$loc <- td$phy$edge.length[H_env$sb]/2
pdf("../output/envRecon.pdf", height=20, width=10); plotBayoupars(H_env, td$phy, cex=0.5); dev.off()
plotBayoupars(H_env, td$phy, cex=0.5)

## Generate hypothesized history of thermy on phylogeny
plot(td$phy, type="fan", show.tip.label=FALSE)
tiplabels(pch=21, bg=pred$Thermy)
tipThermy <- as.numeric(factor(pred$Thermy, levels=c("Ectotherm", "Endotherm")))
likThermy <- make.mkn(td$phy, setNames(tipThermy, td$phy$tip.label),2)
mleThermy <- find.mle(likThermy, c(0.01, 0.01), method="optim")
thermJoint <- round(apply(asr.joint(likThermy, mleThermy$par, n=100),2,mean),0)
plot(td$phy, show.tip.label=FALSE, type="fan")
nodelabels(pch=21, bg=thermJoint)
tiplabels(pch=21, bg=tipThermy, col=tipThermy)
statesTherm <- c(tipThermy, thermJoint)
H_therm <- list(sb=which(statesTherm[td$phy$edge[,1]] != statesTherm[td$phy$edge[,2]]))
H_therm$k <- length(H_therm$sb)
H_therm$ntheta <- 2
H_therm$t2 <- statesTherm[td$phy$edge[H_therm$sb,2]]
H_therm$loc <- td$phy$edge.length[H_therm$sb]/2
pdf("../output/thermRecon.pdf", height=20, width=10); plotBayoupars(H_therm, td$phy, cex=0.5); dev.off()
plotBayoupars(H_therm, td$phy, cex=0.5)

## Generate hypothesized history of trophic group on phylogeny
plot(td$phy, type="fan", show.tip.label=FALSE)
tiplabels(pch=21, bg=pred$TrophicGroup)
tipTrophicGroup <- as.numeric(factor(pred$TrophicGroup, levels=c("herbivore", "omnivore", "carnivore")))
likTrophicGroup <- make.mkn(td$phy, setNames(tipTrophicGroup, td$phy$tip.label),3)
#likTrophicGroup <- constrain(likTrophicGroup, q13~0, q31~0)
mleTrophicGroup <- find.mle(likTrophicGroup, c(0.01, 0.01, 0.01, 0.01,0.01, 0.01), method="optim")
trophJoint <- round(apply(asr.joint(likTrophicGroup, mleTrophicGroup$par, n=100),2,mean),0)
plot(td$phy, show.tip.label=FALSE, type="fan")
nodelabels(pch=21, bg=trophJoint)
tiplabels(pch=21, bg=tipTrophicGroup, col=tipTrophicGroup)
statestroph <- c(tipTrophicGroup, trophJoint)
H_troph <- list(sb=which(statestroph[td$phy$edge[,1]] != statestroph[td$phy$edge[,2]]))
H_troph$k <- length(H_troph$sb)
H_troph$ntheta <- 3
H_troph$t2 <- statestroph[td$phy$edge[H_troph$sb,2]]
H_troph$loc <- td$phy$edge.length[H_troph$sb]/2
pdf("../output/trophRecon.pdf", height=20, width=10); plotBayoupars(H_troph, td$phy, cex=0.5); dev.off()
plotBayoupars(H_troph, td$phy, cex=0.5)

pars <- list()
pars$w <- matrix(rep(1,4)/4, nrow=1)

H4.lik <- function(pars, cache, X, model="Custom"){
  n <- cache$n
  X <- cache$dat
  pred <- cache$pred
  EX.map1 <- bayou:::C_weightmatrix(cache, pars)$E
  
  pars_env <- list(alpha=pars$alpha, sig2=pars$sig2, k=H_env$k, ntheta=H_env$ntheta, theta=c(pars$thEnv1, pars$thEnv2, pars$thEnv3), sb=H_env$sb, loc=H_env$loc, t2=H_env$t2)
  EX.map2 <- bayou:::C_weightmatrix(cache, pars_env)$E  
  
  
  pars_therm <- list(alpha=pars$alpha, sig2=pars$sig2, k=H_therm$k, ntheta=H_therm$ntheta, theta=c(pars$thTherm1, pars$thTherm2), sb=H_therm$sb, loc=H_therm$loc, t2=H_therm$t2)
  EX.map3 <- bayou:::C_weightmatrix(cache, pars_therm)$E 
  
  pars_troph <- list(alpha=pars$alpha, sig2=pars$sig2, k=H_troph$k, ntheta=H_troph$ntheta, theta=c(pars$thTroph1, pars$thTroph2, pars$thTroph3), sb=H_troph$sb, loc=H_troph$loc, t2=H_troph$t2)
  EX.map4 <- bayou:::C_weightmatrix(cache, pars_troph)$E 
  
  X.c <- X - pars$beta_LnMass*pred$LnMass - pars$w[1]*(EX.map1) - pars$w[2]*(EX.map2) - pars$w[3]*(EX.map3) - pars$w[4]*(EX.map4)
  ## This part adds the endothermy parameter to the theta for Mammal and Bird branches
  #dpars <- pars1
  #dpars$theta[dpars$t2[which(dpars$sb %in% c(857, 1719))]] <- dpars$theta[dpars$t2[which(dpars$sb %in% c(857, 1719))]]+dpars$endo
  ### The part below mostly does not change
  transf.phy <- bayou:::C_transf_branch_lengths(cache, 1, X.c, pars$alpha)
  transf.phy$edge.length[cache$externalEdge] <- transf.phy$edge[cache$externalEdge] + cache$SE[cache$phy$edge[cache$externalEdge, 2]]^2*(2*pars$alpha)/pars$sig2
  comp <- bayou:::C_threepoint(list(n=n, N=cache$N, anc=cache$phy$edge[, 1], des=cache$phy$edge[, 2], diagMatrix=transf.phy$diagMatrix, P=X.c, root=transf.phy$root.edge, len=transf.phy$edge.length))
  if(pars$alpha==0){
    inv.yVy <- comp$PP
    detV <- comp$logd
  } else {
    inv.yVy <- comp$PP*(2*pars$alpha)/(pars$sig2)
    detV <- comp$logd+n*log(pars$sig2/(2*pars$alpha))
  }
  llh <- -0.5*(n*log(2*pi)+detV+inv.yVy)
  #llh <- llh + gs.lik(c(pars$pred.sig2, pars$pred.root), root=ROOT.GIVEN) #$impute
  return(list(loglik=llh, theta=pars$theta,resid=X.c, comp=comp, transf.phy=transf.phy))
}

H4.monitor = function(i, lik, pr, pars, accept, accept.type, j){
  names <- c("gen", "lnL", "prior", "alpha","sig2", "rtheta", "rEnv", "rTherm", "rTG","wEnv", "wTherm", "wTG", "k")
  string <- "%-8i%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8.2f%-8i"
  acceptratios <- tapply(accept, accept.type, mean)
  names <- c(names, names(acceptratios))
  if(j==0){
    cat(sprintf("%-7.7s", names), "\n", sep=" ")                           
  }
  cat(sprintf(string, i, lik, pr, pars$alpha, pars$sig2, pars$theta[1], pars$thEnv1, pars$thTherm1, pars$thTroph1, pars$w[2], pars$w[3], pars$w[4], pars$k), sprintf("%-8.2f", acceptratios),"\n", sep="")
}

model.weightedH4 <- list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", beta_LnMass=".slidingWindowProposal",
                                      k=".splitmergePrior", theta=".adjustTheta", thEnv1=".slidingWindowProposal", 
                                      thEnv2=".slidingWindowProposal", thEnv3=".slidingWindowProposal",
                                      thTherm1=".slidingWindowProposal", thTherm2=".slidingWindowProposal",
                                      thTroph1=".slidingWindowProposal", thTroph2=".slidingWindowProposal",
                                      thTroph3=".slidingWindowProposal",
                                      w= ".dirichletMove", slide=".slide"),
control.weights = list(alpha=5, sig2=3, beta_LnMass=3, theta=8, thEnv1=3,
                       thEnv2=3, thEnv3=3, thTherm1=3, thTherm2=3, thTroph1=3, thTroph2=3, thTroph3=3, k=10, w=3, slide=1),
D = list(alpha=0.5, sig2= 0.5, k=1, beta_LnMass=0.02, theta=0.5, thEnv1=0.5,thEnv2=0.5, thEnv3=0.5, thTherm1=0.5, thTherm2=0.5, thTroph1=0.5, thTroph2=0.5, thTroph3=0.5, w =100, slide=1),
parorder = c("alpha", "sig2", "w", "beta_LnMass","thEnv1", "thEnv2","thEnv3", "thTherm1", "thTherm2","thTroph1", "thTroph2", "thTroph3", "k", "ntheta",  "theta", "sb", "loc", "t2"),
rjpars = c("theta"),
shiftpars = c("sb", "loc", "t2"),
monitor.fn = H4.monitor,
lik.fn = H4.lik
) 

prior.wH4 <- make.prior(td$phy, plot.prior = FALSE, 
                        dists=list(dalpha="dhalfcauchy", dsig2="dhalfcauchy",dbeta_LnMass="dnorm",
                                   dk="cdpois", dtheta="dnorm", dthEnv1="dnorm",
                                   dthEnv2="dnorm", dthEnv3="dnorm", dthTherm1="dnorm", dthTherm2="dnorm",
                                   dthTroph1="dnorm", dthTroph2="dnorm", dthTroph3="dnorm",
                                   dw="ddirichlet", dsb="dsb", dloc="dloc"
                        ), 
                        param=list(dalpha=list(scale=0.1), dsig2=list(scale=0.1), 
                                   dbeta_LnMass=list(mean=0.1, sd=0.25),
                                   dk=list(lambda=0.5, kmax=10), dtheta=list(mean=0, sd=2.5),
                                   dthEnv1=list(mean=0, sd=2.5), 
                                   dthEnv2=list(mean=0, sd=2.5), 
                                   dthEnv3=list(mean=0, sd=2.5), 
                                   dthTherm1=list(mean=0, sd=2.5),
                                   dthTherm2=list(mean=0, sd=2.5),
                                   dthTroph1=list(mean=0, sd=2.5),
                                   dthTroph2=list(mean=0, sd=2.5),
                                   dthTroph3=list(mean=0, sd=2.5),
                                   dw=list(alpha=c(1,1,1,1)/4), dsb=list(bmax=1, prob=1), dloc=list(min=0, max=1)
                        )
                        #fixed=list(sb=startpar$sb, k=startpar$k)
)


startpar <- list(alpha=1, sig2=1, beta_LnMass=0.14, k=0, ntheta=1, theta=-1, 
                 thEnv1=-0.5, thEnv2=0, thEnv3=1.5, 
                 thTherm1=-1, thTherm2=0, 
                 thTroph1=-1, thTroph2=0, thTroph3=-0.5, 
                 w=matrix(c(0.25, 0.25, 0.25, 0.25), nrow=1), 
                 sb=numeric(0),
                 t2=numeric(0),
                 loc=numeric(0))

.pred <- select(pred, LnMass)
cache <- bayou:::.prepare.ou.univariate(td$phy, dat, pred=.pred)

model.weightedH4$lik.fn(startpar, cache, dat)$loglik
prior.wH4(startpar, cache)

mcmc.wH4 <- bayou.makeMCMC(td$phy, dat, pred=.pred, SE=0, prior=prior.wH4, model=model.weightedH4, startpar=startpar,
                               new.dir="../output/wH4/", outname="wH4_r001", plot.freq=NULL, samp=100, ticker.freq=1000) # Set up the MCMC


mcmc.wH4$run(10000)












###########
      
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
