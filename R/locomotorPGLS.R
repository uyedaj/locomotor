---
title: "PGLS locomtor performance"
output: html_notebook
---

Phylogenetic analysis of locomotor performance across all life.  

```{r}
library(phylolm)
library(treeplyr)
library(foreach)
library(doParallel)
registerDoParallel(cores=50)
tds <- readRDS("../output/fullTreeTreeDataObjects.rds")
tds[[1]]
```
  
*Single variable fits w/Brownian Motion*
  
```{r}
par(mfrow=c(3,2))
singlefits <- foreach(i=1:length(tds)) %do% {
  td <- filter(tds[[i]], !is.na(TempTraitValue), !is.na(LnMass))
  td$dat <- mutate(td$dat, "LnTempTraitValue"=log(TempTraitValue))
  LnTempTraitValue <- td$dat[['LnTempTraitValue']]
  LnMass <- td$dat[['LnMass']]
  plot(LnMass, LnTempTraitValue)
  phylolm(LnTempTraitValue~LnMass, phy=td$phy, model="BM")
  
}

```


```{r}
names(singlefits) <- sapply(tds, function(x) x$dat$Locomotor.Trait[1])
lapply(singlefits, summary)

```
*Quadratic fits*

```{r}

par(mfrow=c(3,2))
quadfits <- foreach(i=1:length(tds)) %do% {
  td <- filter(tds[[i]], !is.na(TempTraitValue), !is.na(LnMass))
  td$dat <- mutate(td$dat, "LnTempTraitValue"=log(TempTraitValue))
  LnTempTraitValue <- td$dat[['LnTempTraitValue']]
  LnMass <- td$dat[['LnMass']]
  LnMass2 <- LnMass^2
  phylolm(LnTempTraitValue~LnMass + LnMass2, phy=td$phy, model="BM")
  
}


```

```{r}
names(quadfits) <- sapply(tds, function(x) x$dat$Locomotor.Trait[1])
lapply(quadfits, summary)
```

*Environmental effects*

```{r}
par(mfrow=c(3,2))
varenv <- which(sapply(tds, function(x) length(unique(x[['Environment']]))>1))

envfits <- foreach(i=1:length(tds)) %do% {
  if(i %in% varenv){
  td <- filter(tds[[i]], !is.na(TempTraitValue), !is.na(LnMass))
  td$dat <- mutate(td$dat, "LnTempTraitValue"=log(TempTraitValue))
  LnTempTraitValue <- td$dat[['LnTempTraitValue']]
  LnMass <- td$dat[['LnMass']]
  Environment <- td$dat[['Environment']]
  phylolm(LnTempTraitValue~LnMass+Environment, phy=td$phy, model="BM")
  } else {
    list()
  }
  
}
names(envfits) <-  sapply(tds, function(x) x$dat$Locomotor.Trait[1])
```

```{r}
lapply(envfits, summary)
```

*Interaction with environment*

```{r}
par(mfrow=c(3,2))
varenv <- which(sapply(tds, function(x) length(unique(x[['Environment']]))>1))

envIntfits <- foreach(i=1:length(tds)) %do% {
  if(i %in% varenv){
  td <- filter(tds[[i]], !is.na(TempTraitValue), !is.na(LnMass))
  td$dat <- mutate(td$dat, "LnTempTraitValue"=log(TempTraitValue))
  LnTempTraitValue <- td$dat[['LnTempTraitValue']]
  LnMass <- td$dat[['LnMass']]
  Environment <- td$dat[['Environment']]
  phylolm(LnTempTraitValue~LnMass*Environment, phy=td$phy, model="BM")
  } else {
    list()
  }
  
}
names(envIntfits) <-  sapply(tds, function(x) x$dat$Locomotor.Trait[1])
```

```{r}
lapply(envIntfits, summary)
```


```{r}
par(mfrow=c(3,2))
varenv <- which(sapply(tds, function(x) length(unique(x[['Environment']]))>1))

envQuadfits <- foreach(i=1:length(tds)) %do% {
  if(i %in% varenv){
  td <- filter(tds[[i]], !is.na(TempTraitValue), !is.na(LnMass))
  td$dat <- mutate(td$dat, "LnTempTraitValue"=log(TempTraitValue))
  LnTempTraitValue <- td$dat[['LnTempTraitValue']]
  LnMass <- td$dat[['LnMass']]
  LnMass2 <- LnMass^2
  Environment <- td$dat[['Environment']]
  phylolm(LnTempTraitValue~LnMass+LnMass2+Environment, phy=td$phy, model="BM")
  } else {
    list()
  }
  
}
names(envQuadfits) <-  sapply(tds, function(x) x$dat$Locomotor.Trait[1])
```


```{r}
lapply(envQuadfits, summary)
```


*Determine which is the best*

```{r}
aic1 <- sapply(singlefits, function(x) x$aic)
aic2 <- sapply(quadfits, function(x) x$aic)
aic3 <- sapply(envfits, function(x) x$aic)
aic4 <- sapply(envIntfits, function(x) x$aic)
aic5 <- sapply(envQuadfits, function(x) x$aic)
aicTable <- do.call(rbind, list(aic1, aic2, aic3, aic4, aic5))
rownames(aicTable) <- c("singlefits", "quadfits", "envfits", "envIntfits", "envQuadfits")

aicTable
```

For the following models, it looks like these are the best fits:
Angular speed: Quadratic w/ lnMass
Maximum acceleration: Quadratic w/ lnMass
Maximum Speed: Separate environment intercepts
Minimum turn radius: Linear with lnMass
Routine Acceleration: Linear with lnMass
Routine Speed: Separate environment intercepts


*Compare different residual structures*

Make it one big function
```{r}
fitModels <- function(model){
  singlefits <- foreach(i=1:length(tds)) %do% {
    td <- filter(tds[[i]], !is.na(TempTraitValue), !is.na(LnMass))
    td$dat <- mutate(td$dat, "LnTempTraitValue"=log(TempTraitValue))
    LnTempTraitValue <- td$dat[['LnTempTraitValue']]
    LnMass <- td$dat[['LnMass']]
    phylolm(LnTempTraitValue~LnMass, phy=td$phy, model=model)
  }
  quadfits <- foreach(i=1:length(tds)) %do% {
    td <- filter(tds[[i]], !is.na(TempTraitValue), !is.na(LnMass))
    td$dat <- mutate(td$dat, "LnTempTraitValue"=log(TempTraitValue))
    LnTempTraitValue <- td$dat[['LnTempTraitValue']]
    LnMass <- td$dat[['LnMass']]
    LnMass2 <- LnMass^2
    phylolm(LnTempTraitValue~LnMass + LnMass2, phy=td$phy, model=model)
  }
  names(singlefits) <- names(quadfits) <-sapply(tds, function(x) x$dat$Locomotor.Trait[1])
  
varenv <- which(sapply(tds, function(x) length(unique(x[['Environment']]))>1))
  
envfits <- foreach(i=1:length(tds)) %do% {
  if(i %in% varenv){
  td <- filter(tds[[i]], !is.na(TempTraitValue), !is.na(LnMass))
  td$dat <- mutate(td$dat, "LnTempTraitValue"=log(TempTraitValue))
  LnTempTraitValue <- td$dat[['LnTempTraitValue']]
  LnMass <- td$dat[['LnMass']]
  Environment <- td$dat[['Environment']]
  phylolm(LnTempTraitValue~LnMass+Environment, phy=td$phy, model=model)
  } else {
    list()
  }
  
}
names(envfits) <-  sapply(tds, function(x) x$dat$Locomotor.Trait[1])

envIntfits <- foreach(i=1:length(tds)) %do% {
  if(i %in% varenv){
  td <- filter(tds[[i]], !is.na(TempTraitValue), !is.na(LnMass))
  td$dat <- mutate(td$dat, "LnTempTraitValue"=log(TempTraitValue))
  LnTempTraitValue <- td$dat[['LnTempTraitValue']]
  LnMass <- td$dat[['LnMass']]
  Environment <- td$dat[['Environment']]
  phylolm(LnTempTraitValue~LnMass*Environment, phy=td$phy, model=model)
  } else {
    list()
  }
  
}
names(envIntfits) <-  sapply(tds, function(x) x$dat$Locomotor.Trait[1])

envQuadfits <- foreach(i=1:length(tds)) %do% {
  if(i %in% varenv){
  td <- filter(tds[[i]], !is.na(TempTraitValue), !is.na(LnMass))
  td$dat <- mutate(td$dat, "LnTempTraitValue"=log(TempTraitValue))
  LnTempTraitValue <- td$dat[['LnTempTraitValue']]
  LnMass <- td$dat[['LnMass']]
  LnMass2 <- LnMass^2
  Environment <- td$dat[['Environment']]
  phylolm(LnTempTraitValue~LnMass+LnMass2+Environment, phy=td$phy, model=model)
  } else {
    list()
  }
  
}
names(envQuadfits) <-  sapply(tds, function(x) x$dat$Locomotor.Trait[1])

return(list(singlefits=singlefits, quadfits=quadfits, envfits=envfits, envIntfits=envIntfits, envQuadfits=envQuadfits))
}







```

Run it over different error structures

```{r}

fitsOUr <- fitModels("OUrandomRoot")
fitsOUf <- fitModels("OUfixedRoot")
fitsLambda <- fitModels("lambda")

```


```{r}
print("BM")
aicTable

print("OUrandomRoot")
do.call(rbind, lapply(fitsOUr, function(x) sapply(x, function(y) y$aic)))

print("OUfixedRoot")
do.call(rbind, lapply(fitsOUf, function(x) sapply(x, function(y) y$aic)))

print("Lambda")
do.call(rbind, lapply(fitsLambda, function(x) sapply(x, function(y) y$aic)))

```

For the following models, it looks like these are the best fits:
Angular speed: Quadratic w/OU residuals
Maximum acceleration: Quadratic & environment w/ OU residuals
Maximum Speed: Quadratic & environment w/Lambda residuals
Minimum turn radius: Environment with either OU or Lambda residuals
Routine Acceleration: Linear with BM residuals
Routine Speed: Environment and Quadratic with OU residuals

*Plotting fitted models*

```{r}
## Angular Speed
fit <- fitsOUr$quadfits$`Angular Speed`
plot(fit$X[,2], fit$y, main="Angular Speed", pch=21, bg=tds[[1]]$dat$Environment)
o <- order(fit$X[,2])
lines(fit$X[o,2], fit$fitted.values[o], lty=2)
coef <- round(fit$coefficients,2)
text(min(fit$X[,2]), min(fit$y), labels=paste(coef[1], " + ",coef[2], "*x + ", coef[3],"*x^2", sep=""))

## Maximum Acceleration
fit <- fitsOUr$envQuadfits$`Maximum Acceleration`
plot(fit$X[,2], fit$y, main="Maximum Acceleration", pch=21, bg=tds[[2]]$dat$Environment)
o <- order(fit$X[,2])
coef <- fit$coefficients
curve(coef[1]+coef[2]*x+coef[3]*x^2, add=TRUE, lty=2, col=1)
curve(coef[1]+coef[2]*x+coef[3]*x^2+coef[4], add=TRUE, lty=2, col=2)
curve(coef[1]+coef[2]*x+coef[3]*x^2+coef[5], add=TRUE, lty=2, col=3)

## Maximum Speed
fit <- fitsLambda$envQuadfits$`Maximum Speed`
plot(fit$X[,2], fit$y, main="Maximum Speed", pch=21, bg=tds[[3]]$dat$Environment)
o <- order(fit$X[,2])
coef <- fit$coefficients
curve(coef[1]+coef[2]*x+coef[3]*x^2, add=TRUE, lty=2, col=1)
curve(coef[1]+coef[2]*x+coef[3]*x^2+coef[4], add=TRUE, lty=2, col=2)
curve(coef[1]+coef[2]*x+coef[3]*x^2+coef[5], add=TRUE, lty=2, col=3)

## Minimum turn radius
fit <- fitsLambda$envfits$`Minimum Turn Radius`
plot(fit$X[,2], fit$y, main="Minimum Turn Radius", pch=21, bg=tds[[4]]$dat$Environment)
o <- order(fit$X[,2])
coef <- fit$coefficients
curve(coef[1]+coef[2]*x, add=TRUE, lty=2, col=1)
curve(coef[1]+coef[2]*x+coef[3], add=TRUE, lty=2, col=2)
curve(coef[1]+coef[2]*x+coef[4], add=TRUE, lty=2, col=3)

## Routine acceleration
fit <- singlefits$`Routine Accleration`
plot(fit$X[,2], fit$y, main="Routine Acceleration", pch=21, bg=tds[[5]]$dat$Environment)
o <- order(fit$X[,2])
coef <- fit$coefficients
curve(coef[1]+coef[2]*x, add=TRUE, lty=2, col=1)

## Routine Speed
fit <- fitsOUr$envQuadfits$`Routine Speed`
plot(fit$X[,2], fit$y, main="Routine Speed", pch=21, bg=tds[[6]]$dat$Environment)
o <- order(fit$X[,2])
coef <- fit$coefficients
curve(coef[1]+coef[2]*x+coef[3]*x^2, add=TRUE, lty=2, col=1)
curve(coef[1]+coef[2]*x+coef[3]*x^2+coef[4], add=TRUE, lty=2, col=2)
curve(coef[1]+coef[2]*x+coef[3]*x^2+coef[5], add=TRUE, lty=2, col=3)



```




