library(treeplyr)
#devtools::install_github("uyedaj/bayou", ref="dev")
library(bayou)
library(rotl)
library(phytools)
library(Matrix)
library(phylolm)
library(treetimer)
library(phyndr)

data("ttolData")

setwd("~/repos/locomotor/R/")

tree <- read.tree(file="../data/tetrapods.tre")

locomotor <- read.csv("../data/Locomotor_Performance_Database_002.csv")
locomotor <- mutate(locomotor, genspec = gsub(" " ,"_", TaxaID))
locomotor$LnMass[locomotor$LnMass=="Err:502"] <- NA
locomotor$Thermy[locomotor$Thermy=="ectotherm"] <- "Ectotherm"
locomotor$Thermy[locomotor$Thermy=="endotherm"] <- "Endotherm"
td <- make.treedata(ttolData$phy, locomotor)
taxalist <- unique(locomotor$genspec)
ttolData$phy$edge.length[ttolData$phy$edge.length==0] <- .Machine$double.eps
nH <- nodeHeights(ttolData$phy)
TL <- max(nH)
externalEdge <- which(ttolData$phy$edge[,2] <= length(ttolData$phy$tip.label))
diff <- TL - nH[externalEdge, 2]
ttolData$phy$edge.length[externalEdge] <- ttolData$phy$edge.length[externalEdge] + diff
is.ultrametric(ttolData$phy)
phyndtree <-phyndr::phyndr_genus(ttolData$phy, taxalist[!is.na(taxalist)])

td <- make.treedata(ttolData$phy, locomotor)
ptree <- phyndr::phyndr_sample(phyndtree)
tmp <- make.treedata(ptree, locomotor)
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

saveRDS(tds, "../output/tds_ttol.rds")

