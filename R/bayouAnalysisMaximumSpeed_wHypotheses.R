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
td <- filter(td, TempKelvin < 280)
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
#pdf("../output/envRecon.pdf", height=20, width=10); plotBayoupars(H_env, td$phy, cex=0.5); dev.off()
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
#pdf("../output/thermRecon.pdf", height=20, width=10); plotBayoupars(H_therm, td$phy, cex=0.5); dev.off()
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
#pdf("../output/trophRecon.pdf", height=20, width=10); plotBayoupars(H_troph, td$phy, cex=0.5); dev.off()
plotBayoupars(H_troph, td$phy, cex=0.5)

pars <- list()
pars$w <- matrix(rep(1,4)/4, nrow=1)
cache <- bayou:::.prepare.ou.univariate(td$phy, dat, SE=0.05, pred=select(pred, LnMass))

H4.lik <- function(pars, cache, X, model="Custom"){
  n <- cache$n
  X <- cache$dat
  pred <- cache$pred
  W0 <- bayou:::C_weightmatrix(cache, pars)$W
  EX.map0 <- W0 %*% pars$beta_LnMass * pred[,1] + W0 %*% pars$theta
  
  pars_env <- list(alpha=pars$alpha, sig2=pars$sig2, k=H_env$k, ntheta=H_env$ntheta, beta_LnMass = pars$beta_Env, theta=c(pars$thEnv1, pars$thEnv2, pars$thEnv3), sb=H_env$sb, loc=H_env$loc, t2=H_env$t2)
  W1 <- bayou:::C_weightmatrix(cache, pars_env)$W
  EX.map1 <- W1 %*% pars_env$beta_LnMass * pred[,1] + W1 %*% pars_env$theta
  
  
  pars_therm <- list(alpha=pars$alpha, sig2=pars$sig2, k=H_therm$k, ntheta=H_therm$ntheta, beta_LnMass = pars$beta_Therm, theta=c(pars$thTherm1, pars$thTherm2), sb=H_therm$sb, loc=H_therm$loc, t2=H_therm$t2)
  W2 <- bayou:::C_weightmatrix(cache, pars_therm)$W
  EX.map2 <- W2 %*% pars_therm$beta_LnMass * pred[,1] + W2 %*% pars_therm$theta
  
  pars_troph <- list(alpha=pars$alpha, sig2=pars$sig2, k=H_troph$k, ntheta=H_troph$ntheta, beta_LnMass = pars$beta_Troph, theta=c(pars$thTroph1, pars$thTroph2, pars$thTroph3), sb=H_troph$sb, loc=H_troph$loc, t2=H_troph$t2)
  W3 <- bayou:::C_weightmatrix(cache, pars_troph)$W
  EX.map3 <- W3 %*% pars_troph$beta_LnMass * pred[,1] + W3 %*% pars_troph$theta
  
  X.c <- X - pars$w[1]*EX.map0 - pars$w[2]*EX.map1 - pars$w[3]*EX.map2 - pars$w[4]*EX.map3
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

model.weightedH4 <- list(moves = list(alpha=".multiplierProposal", sig2=".multiplierProposal", beta_LnMass=".vectorSlidingWindow",
                                      beta_Env=".vectorSlidingWindow", beta_Therm=".vectorSlidingWindow", beta_Troph=".vectorSlidingWindow",
                                      k=".splitmergePrior", theta=".adjustTheta", thEnv1=".slidingWindowProposal", 
                                      thEnv2=".slidingWindowProposal", thEnv3=".slidingWindowProposal",
                                      thTherm1=".slidingWindowProposal", thTherm2=".slidingWindowProposal",
                                      thTroph1=".slidingWindowProposal", thTroph2=".slidingWindowProposal",
                                      thTroph3=".slidingWindowProposal",
                                      w= ".dirichletMove", slide=".slide"),
control.weights = list(alpha=5, sig2=3, beta_LnMass=5, beta_Env=5, beta_Therm=5, beta_Troph=5, theta=8, thEnv1=3,
                       thEnv2=3, thEnv3=3, thTherm1=3, thTherm2=3, thTroph1=3, thTroph2=3, thTroph3=3, k=10, w=3, slide=1),
D = list(alpha=0.7, sig2= 0.7, k=c(1,1), beta_LnMass=0.1, beta_Env=0.1, beta_Therm=0.1, beta_Troph=0.1, theta=0.5, thEnv1=0.75,thEnv2=0.75, thEnv3=0.75, thTherm1=0.5, thTherm2=0.5, thTroph1=0.7, thTroph2=0.7, thTroph3=0.7, w =100, slide=1),
parorder = c("alpha", "sig2", "w", "beta_Env", "beta_Therm", "beta_Troph","thEnv1", "thEnv2","thEnv3", "thTherm1", "thTherm2","thTroph1", "thTroph2", "thTroph3", "k", "ntheta", "beta_LnMass",  "theta", "sb", "loc", "t2"),
rjpars = c("theta", "beta_LnMass"),
shiftpars = c("sb", "loc", "t2"),
monitor.fn = H4.monitor,
call = dat ~ LnMass,
lik.fn = H4.lik
) 

prior.wH4 <- make.prior(td$phy, plot.prior = FALSE, 
                        dists=list(dalpha="dlnorm", dsig2="dhalfcauchy",dbeta_LnMass="dnorm", 
                                   dbeta_Env="dnorm", dbeta_Therm="dnorm", dbeta_Troph="dnorm",
                                   dk="cdpois", dtheta="dnorm", dthEnv1="dnorm",
                                   dthEnv2="dnorm", dthEnv3="dnorm", dthTherm1="dnorm", dthTherm2="dnorm",
                                   dthTroph1="dnorm", dthTroph2="dnorm", dthTroph3="dnorm",
                                   dw="ddirichlet", dsb="dsb", dloc="dloc"
                        ), 
                        param=list(dalpha=list(meanlog=-3.5, sdlog=2.5), dsig2=list(scale=0.1), 
                                   dbeta_LnMass=list(mean=0.1, sd=0.4),
                                   dbeta_Env=list(mean=0.1, sd=0.4),
                                   dbeta_Therm=list(mean=0.1, sd=0.4),
                                   dbeta_Troph=list(mean=0.1, sd=0.4),
                                   dk=list(lambda=2, kmax=50), 
                                   dtheta=list(mean=-1, sd=3),
                                   dthEnv1=list(mean=-1, sd=3), 
                                   dthEnv2=list(mean=-1, sd=3), 
                                   dthEnv3=list(mean=-1, sd=3), 
                                   dthTherm1=list(mean=-1, sd=3),
                                   dthTherm2=list(mean=-1, sd=3),
                                   dthTroph1=list(mean=-1, sd=3),
                                   dthTroph2=list(mean=-1, sd=3),
                                   dthTroph3=list(mean=-1, sd=3),
                                   dw=list(alpha=c(1,1,1,1)), dsb=list(bmax=1, prob=1), dloc=list(min=0, max=1)
                        )
                        #fixed=list(sb=startpar$sb, k=startpar$k)
)




#.pred <- select(pred, LnMass)
#cache <- bayou:::.prepare.ou.univariate(td$phy, dat, pred=.pred)

#model.weightedH4$lik.fn(startpar, cache, dat)$loglik
#prior.wH4(startpar, cache)
library(foreach)
library(doParallel)
registerDoParallel(cores=50)
.pred <- select(pred, LnMass)
dum <- foreach(qqq=3:12) %dopar% {
  startpar <- priorSim(prior.wH4, tree=td$phy)$pars[[1]]
  startpar$beta_LnMass <- rnorm(startpar$ntheta, 0.1, 0.4)
  startpar$beta_Env <- rnorm(3, 0.1, 0.4)
  startpar$beta_Therm <- rnorm(2, 0.1, 0.4)
  startpar$beta_Troph <- rnorm(3, 0.1, 0.4)
  
  name.x <- paste0("modelwH4_r", qqq, paste0(sample(0:9,4,replace=TRUE), collapse=""))
  mcmc.wH4 <- bayou.makeMCMC(td$phy, dat, pred=.pred, SE=0.05, prior=prior.wH4, model=model.weightedH4, startpar=startpar,
                           new.dir="../output/wH4/", outname=name.x, plot.freq=NULL, samp=100, ticker.freq=1000) # Set up the MCMC
  mcmc.wH4$run(100000)
  saveRDS(mcmc.wH4, file=paste("../output/wH4/mcmc.", name.x, ".rds", sep=""))
  chain.wH4 <- set.burnin(mcmc.wH4$load(), 0.3)
  saveRDS(chain.wH4, file=paste("../output/wH4/chain.", name.x, ".rds", sep=""))

  pdf(paste("../output/figures.", name.x, ".pdf", sep=""), height=20, width=15)
  plot(chain.wH4)
  plotSimmap.mcmc(chain.wH4, pp.cutoff=0.1, cex=0.5)
  plotBranchHeatMap(td$phy, chain.wH4, variable="theta", nn=1000, pal=viridis::viridis)
  plotBranchHeatMap(td$phy, chain.wH4, variable="beta_LnMass", nn=1000, pal=viridis::viridis)
  try({shiftsum.wH4 <- bayou::shiftSummaries(chain.wH4, mcmc.wH4, pp.cutoff=0.2);
  plotShiftSummaries(shiftsum.wH4, pal=viridis::viridis);
  plotShiftSummaries(shiftsum.wH4, pal=viridis::viridis, single.plot=TRUE, label.pts=FALSE)
  })
  par(mfrow=c(2,2))

  wchain <- do.call(rbind, chain.wH4$w)
  densW <- apply(wchain, 2, density)
  plot(densW[[1]], xlim=c(0,1), type="n", ylim=c(0,8), main="Hypothesis weights")
  polycol <- sapply(viridis::viridis(4), function(x) bayou::makeTransparent(x, 100))
  polygon(densW[[1]], col = polycol[1])
  polygon(densW[[2]], col = polycol[2])
  polygon(densW[[3]], col = polycol[3])
  polygon(densW[[4]], col = polycol[4])
  polygon(density(rdirichlet(nrow(wchain), alpha=c(1,1,1,1))[,1]), lty=2)
  legend(0.6,8, legend=c("prior", "Trophic Group", "Thermy", "Environment", "Clade-Specific RJ"),pt.cex=2, pt.bg=c(0,rev(polycol)), pch=c(NA, 22,22,22,22), lty=c(2,0,0,0,0))
  
  apply(wchain, 2, quantile)
  w1 <- wchain[,2]/length(wchain[,2]); w2 <-  wchain[,3]/length(wchain[,3]); w3 <-  wchain[,4]/length(wchain[,4]); w0 <-  wchain[,1]/length(wchain[,1])
  w1 <- w1/sum(w1); w2 <- w2/sum(w2); w3 <- w3/sum(w3); w0 <- w0/sum(w0)
  m0 <- sum(w0*sapply(chain.wH4$beta_LnMass, function(x) x[1]))
  b0 <- sum(w0*sapply(chain.wH4$theta, function(x) x[1]))
  postburn <- floor(0.3*length(chain.wH4$gen)):length(chain.wH4$gen)
  samp <- sort(sample(floor(0.3*length(chain.wH4$gen)):length(chain.wH4$gen), 100, replace=FALSE))
  
  plot(pred$LnMass, dat, pch=21, bg=viridis::viridis(3)[tipEnv], col=viridis::viridis(3)[tipEnv], xlim=c(-10, 20), ylim=c(-6, 4))
  dum <- lapply(samp, function(i) curve(chain.wH4$beta_Env[[i]][1]*x + chain.wH4$thEnv1[[i]], add=TRUE, col=bayou::makeTransparent(viridis::viridis(3)[1], alpha = 100*chain.wH4$w[[i]][2])))
  dum <- lapply(samp, function(i) curve(chain.wH4$beta_Env[[i]][2]*x + chain.wH4$thEnv2[[i]], add=TRUE, col=bayou::makeTransparent(viridis::viridis(3)[2], alpha = 100*chain.wH4$w[[i]][2])))
  dum <- lapply(samp, function(i) curve(chain.wH4$beta_Env[[i]][3]*x + chain.wH4$thEnv3[[i]], add=TRUE, col=bayou::makeTransparent(viridis::viridis(3)[3], alpha = 100*chain.wH4$w[[i]][2])))
  m <- apply(do.call(rbind, chain.wH4$beta_Env[postburn]), 2, function(x) sum(w1[postburn]*x))
  b <- apply(cbind(chain.wH4$thEnv1[postburn], chain.wH4$thEnv2[postburn], chain.wH4$thEnv3[postburn]), 2, function(x) sum(w1[postburn]*x))
  lapply(1:3, function(x) lines(seq(-10,20,length.out=1000), m[x]*seq(-10,20,length.out=1000)+b[x], col=viridis::viridis(3)[x], lwd=2))


  plot(pred$LnMass, dat, pch=21, bg=viridis::viridis(3)[tipThermy])
  dum <- lapply(samp, function(i) curve(chain.wH4$beta_Therm[[i]][1]*x + chain.wH4$thTherm1[[i]], add=TRUE, col=bayou::makeTransparent(viridis::viridis(3)[1], alpha = 100*chain.wH4$w[[i]][3])))
  dum <- lapply(samp, function(i) curve(chain.wH4$beta_Therm[[i]][2]*x + chain.wH4$thTherm2[[i]], add=TRUE, col=bayou::makeTransparent(viridis::viridis(3)[2], alpha = 100*chain.wH4$w[[i]][3])))
  m <- apply(do.call(rbind, chain.wH4$beta_Therm[postburn]), 2, function(x) sum(w2[postburn]*x))
  b <- apply(cbind(chain.wH4$thTherm1[postburn], chain.wH4$thTherm2[postburn]), 2, function(x) sum(w2[postburn]*x))
  lapply(1:2, function(x) lines(seq(-10,20,length.out=1000), m[x]*seq(-10,20,length.out=1000)+b[x], col=viridis::viridis(3)[x], lwd=2))

  plot(pred$LnMass, dat, pch=21, bg=viridis::viridis(3)[tipTrophicGroup])
  dum <- lapply(samp, function(i) curve(chain.wH4$beta_Troph[[i]][1]*x + chain.wH4$thTroph1[[i]], add=TRUE, col=bayou::makeTransparent(viridis::viridis(3)[1], alpha = 100*chain.wH4$w[[i]][3])))
  dum <- lapply(samp, function(i) curve(chain.wH4$beta_Troph[[i]][2]*x + chain.wH4$thTroph2[[i]], add=TRUE, col=bayou::makeTransparent(viridis::viridis(3)[2], alpha = 100*chain.wH4$w[[i]][3])))
  dum <- lapply(samp, function(i) curve(chain.wH4$beta_Troph[[i]][3]*x + chain.wH4$thTroph3[[i]], add=TRUE, col=bayou::makeTransparent(viridis::viridis(3)[3], alpha = 100*chain.wH4$w[[i]][3])))
  m <- apply(do.call(rbind, chain.wH4$beta_Troph[postburn]), 2, function(x) sum(w3[postburn]*x))
  b <- apply(cbind(chain.wH4$thTroph1[postburn], chain.wH4$thTroph2[postburn], chain.wH4$thTroph3[postburn]), 2, function(x) sum(w3[postburn]*x))
  lapply(1:3, function(x) lines(seq(-10,20,length.out=1000), m[x]*seq(-10,20,length.out=1000)+b[x], col=viridis::viridis(3)[x], lwd=2))

  dev.off()
}
