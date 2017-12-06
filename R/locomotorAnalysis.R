library(treeplyr)
devtools::install_github("uyedaj/bayou", ref="dev")
library(bayou)
library(rotl)
library(phytools)
library(Matrix)
library(phylolm)
setwd("~/repos/locomotor/R/")

tree <- read.tree(file="../data/tetrapods.tre")

locomotor <- read.csv("../data/Locomotor_Performance_Database_002.csv")
locomotor <- mutate(locomotor, genspec = gsub(" " ,"_", TaxaID))
locomotor$LnMass[locomotor$LnMass=="Err:502"] <- NA
locomotor$Thermy[locomotor$Thermy=="ectotherm"] <- "Ectotherm"
locomotor$Thermy[locomotor$Thermy=="endotherm"] <- "Endotherm"
td <- make.treedata(tree, locomotor)

traits <- levels(locomotor$Locomotor.Trait)

tds <- list()
for(i in 1:length(traits)){
  .data <- filter(locomotor, Locomotor.Trait==traits[i])
  tds[[i]] <- try(make.treedata(tree, .data))
  if(class(tds[[i]])=="try-error"){
   tds[[i]] <- .data
  }
}

names(tds) <- traits
plot(tds$`Maximum Speed`[['LnMass']], tds$`Maximum Speed`[['TempTraitValue']], pch=21, bg=as.numeric(tds$`Maximum Speed`[['Order']])) 


do.call(rbind, lapply(tds, function(x) length(x$phy$tip.label)))

saveRDS(tds, "../output/tds.rds")

