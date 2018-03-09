library(treeplyr)
#devtools::install_github("uyedaj/bayou", ref="dev")
library(bayou)
library(rotl)
library(phytools)
library(Matrix)
library(phylolm)
library(treetimer)
library(phyndr)

getOttIds <- function (taxalist, ncores = 1, context = NULL) {
  scipen <- options()$scipen
  digits <- options()$digits
  options(scipen = 100, digits = 4)
  .taxalist <- gsub("_", " ", taxalist)
  #.taxalist <- gsub(" sp$", "", .taxalist)
  tax <- parallel::mclapply(1:length(taxalist), function(i) try(rotl::tnrs_match_names(.taxalist[i], 
                                                                                       do_approximate_matching = FALSE, context_name = context)), 
                            mc.cores = ncores)
  failed <- which(sapply(tax, function(x) class(x)[1] == "try-error"))
  if (length(failed) > 0) {
    tax[failed] <- parallel::mclapply(failed, function(i) try(rotl::tnrs_match_names(.taxalist[i], 
                                                                                     do_approximate_matching = TRUE, context_name = context)), 
                                      mc.cores = ncores)
  }
  stillfailed <- which(sapply(tax, function(x) if (class(x)[1] == 
                                                   "try-error") {
    TRUE
  }
  else {
    is.na(x$ott_id)
  }))
  if (length(stillfailed > 0)) {
    tax[stillfailed] <- lapply(stillfailed, function(x) data.frame(search_string = .taxalist[x], 
                                                                   unique_name = .taxalist[x], approximate_match = NA, 
                                                                   ott_id = NA, is_synonym = NA, flags = NA, number_matches = 0))
  }
  tax <- do.call(rbind, tax)
  genspec <- unname(sapply(tax[, 2], function(x) paste(strsplit(x, 
                                                                split = " ")[[1]][1:2], collapse = " ")))
  genspec <- gsub(" (genus", " sp.", genspec, fixed = TRUE)
  genspec <- gsub(" NA", " sp.", genspec, fixed = TRUE)
  if (sum(duplicated(genspec)) > 0) {
    cat("Dropping duplicated taxa: ", paste(taxalist[duplicated(genspec)], 
                                            collapse = ", "), "\n")
  }
  if (sum(is.na(tax$ott_id)) > 0) {
    cat("No ott ids found for taxa: ", paste(taxalist[is.na(tax$ott_id)], 
                                             collapse = ", "), "\n")
  }
  tax_unique <- tax[!(duplicated(genspec) | is.na(tax$ott_id)), 
                    ]
  tax_unique$ottids <- as.character(tax_unique$ott_id)
  options(scipen = scipen, digits = digits)
  tax_unique[, 1] <- gsub(" ", "_", tax_unique[, 1])
  tax_unique[, 1] <- sapply(tax_unique[, 1], function(x) simpleCap(x))
  return(tax_unique)
}

data("ttolData")

setwd("~/repos/locomotor/R/")

tree <- read.tree(file="../data/tetrapods.tre")

locomotor <- read.csv("../data/Locomotor_Performance_Database_002.csv")
locomotor <- mutate(locomotor, genspec = gsub(" " ,"_", TaxaID))
locomotor$LnMass[locomotor$LnMass=="Err:502"] <- NA
locomotor$LnMass <- as.numeric(as.character(locomotor$LnMass))
locomotor$Thermy[locomotor$Thermy=="ectotherm"] <- "Ectotherm"
locomotor$Thermy[locomotor$Thermy=="endotherm"] <- "Endotherm"
locomotor <- filter(locomotor, !is.na(genspec))
#td <- make.treedata(ttolData$phy, locomotor)


locomotor$genspec <- gsub(".", "", locomotor$genspec, fixed=TRUE)
locomotor$genspec <- gsub("_Spp$", "", locomotor$genspec)
locomotor$genspec <- gsub("_Sp$", "", locomotor$genspec)
locomotor$genspec <- gsub("_sp$", "", locomotor$genspec)
locomotor$genspec <- gsub("_spp$", "", locomotor$genspec)
locomotor$genspec <- gsub("\\n.*$", "", locomotor$genspec)
locomotor$genspec <- gsub("Mesentotoma", "Haloentomobrya", locomotor$genspec)
locomotor$genspec <- gsub("Meroles", "Ichnotropis", locomotor$genspec)
locomotor$genspec <- gsub("Meroles", "Ichnotropis", locomotor$genspec)
locomotor$genspec <- gsub("_p_", "_pygmaeus_", locomotor$genspec)
locomotor$genspec <- gsub("(.*?_.*?)_.*$", "\\1", locomotor$genspec)

taxalist <- sort(unique(locomotor$genspec))

otts <- getOttIds(taxalist, ncores=80)
ind <- 1:length(taxalist)
notfound <- taxalist[ind[!ind %in% attributes(otts)$row.names]]
notfound

library(foreach)
library(doParallel)
ncores =80
registerDoParallel(cores=ncores)

tree <- read.tree("http://www.timetree.org/public/data/TimetreeOfLife2015.nwk")
refTTOL <- treetimer::makeReferenceTree(tree, cores=ncores)
tiplabels <- gsub(" ", "_", refTTOL$dat$unique_name)
tiplabels <- gsub("(.*?_.*?)_.*$", "\\1", tiplabels)
tiplabels <- gsub("__", "_", tiplabels)
refTTOL$phy$tip.label <- tiplabels
saveRDS(refTTOL, "../output/refTTOL2018.rds")


dataOttTable <- otts
rownames(dataOttTable) <- 1:nrow(dataOttTable)
ttolObject <- refTTOL
tree <- ttolObject$phy
dat <- ttolObject$dat
lineages <- ttolObject$lineages
names(lineages) <- tree$tip.label
#tax <- sliceTaxonomyTable(timeslices, tree, lookupLICAs = FALSE)
rm(ttolObject)
#data_taxa <- gsub(" ", "_", dataOttTable$unique_name)
missing <- setdiff(1:length(taxalist), as.numeric(rownames(dataOttTable)))
nas <- which(is.na(dataOttTable$ott_id))

unique_names <- gsub(" ", "_", dataOttTable$unique_name)
unique_names <- gsub("(.*?_.*?)_.*$", "\\1", unique_names)
old_names <- dataOttTable$unique_name
dataOttTable$unique_name <- unique_names

exact_matches <- which(unique_names %in% tree$tip.label)
exact_tax_matches <- which(tree$tip.label %in% unique_names)
exact_taxa <- unique_names[exact_matches]

extra_tree_taxa <- tree$tip.label[!tree$tip.label %in% exact_taxa]

extra_data_taxa <- which(!unique_names %in% exact_taxa)

extractLineages <- function (otts, ncores, ottnames = NULL) {
  .extractFromRaw <- function(raw) {
    data.frame(do.call(rbind, lapply(raw[[1]]$lineage, function(y) data.frame(otName = y$name, 
                                                                              rank = y$rank, ottid = y$ott_id, unique_name = y$unique_name))))
  }
  lineages_raw <- foreach(j=1:length(otts)) %dopar% {rotl::taxonomy_taxon_info(otts[j], include_lineage = TRUE)}
  lineageTable <- lapply(lineages_raw, .extractFromRaw)
  if (!is.null(ottnames)) {
    names(lineageTable) <- ottnames
  }
  return(lineageTable)
}


dataLineages <- extractLineages(dataOttTable$ott_id[extra_data_taxa], ncores = ncores)
names(dataLineages) <- unique_names[extra_data_taxa]
cbind(1:length(dataLineages), names(dataLineages), sapply(dataLineages, function(x) as.character(x[1,1])))
#

exactLineages <- lineages[exact_taxa]
exactLineageMatches <- foreach(i=1:length(lineages)) %dopar% {
  sapply(exactLineages, function(x) sum(lineages[[i]]$ottid %in% x$ottid))
}
names(exactLineageMatches) <- tree$tip.label

query_newLineages <- foreach(i=1:length(dataLineages)) %dopar% {
  focalTaxa <- unique_names[extra_data_taxa][i]
  dL <- dataLineages[[i]]
  lineageMatches <- sapply(lineages, function(x) sum(dL$ottid %in% x$ottid))
  match_taxa <- tree$tip.label[which(lineageMatches==max(lineageMatches))]
  match_lineages <- lineages[match_taxa]
  exact_comparisons <- lapply(exactLineageMatches[match_taxa], function(x) x[exact_taxa])
  exact_max <- sapply(exact_comparisons, function(x) max(x))
  if(max(lineageMatches) > max(sapply(exact_comparisons, function(x) max(x)))){
    add=TRUE
  } else {
    add=FALSE
  }
  return(list(ott= dataOttTable$ott_id[extra_data_taxa][i],focal=focalTaxa, lineage=dL, taxa=match_taxa, max=max(lineageMatches), exact=exact_max, add=add))
}
names(query_newLineages) <- unique_names[extra_data_taxa]

add_taxa <- which(sapply(query_newLineages, function(x) x$add))

add_taxa_names <- lapply(add_taxa, function(x) query_newLineages[[x]]$taxa)
names(add_taxa_names) <- unique_names[extra_data_taxa][add_taxa]

add_overlap <- foreach(i=1:length(add_taxa_names)) %dopar% {
  sapply(add_taxa_names, function(x) sum(add_taxa_names[[i]] %in% x))
}

count_overlap <- lapply(1:length(add_overlap), function(x) add_overlap[[x]][add_overlap[[x]]>0 & 1:length(add_overlap)!=x])
names(count_overlap) <- names(add_taxa_names)
safe_add <- which(sapply(count_overlap, length)==0)
check_add <- which(sapply(count_overlap, length)>0)

genus_compare <- function(x,y){
  xg <- strsplit(x, "_")[[1]][1]
  yg <- unname(sapply(y, function(i) strsplit(i, "_")[[1]][1]))
  sdist <- stringdist::stringdist(xg, yg)
  y[which(sdist==min(sdist))[1]]
}
repOtt <- unname(sapply(add_taxa, function(x) query_newLineages[[x]]$ott))
replacements <- do.call(rbind, lapply(1:length(add_taxa_names), function(x) data.frame("ottid"=repOtt[x], "taxa"=names(add_taxa_names)[x], "replacement"=genus_compare(names(add_taxa_names)[x], add_taxa_names[[x]]))))
replaceTable <- replacements[!duplicated(replacements$replacement),]

tmp <- dataOttTable[match(replaceTable$ottid, dataOttTable$ott_id),]

alltaxa <- exact_taxa

ptree <- tree
ptree$tip.label[match(replaceTable$replacement, ptree$tip.label)] <- dataOttTable$search_string[match(replaceTable$ottid, dataOttTable$ott_id)]

traits <- levels(locomotor$Locomotor.Trait)

tds <- list()

for(i in 1:length(traits)){
  .data <- filter(locomotor, Locomotor.Trait==traits[i])
  .data <- group_by(.data, genspec) %>% summarize(., SI_Mass=mean(SI_Mass, na.rm=TRUE), LnMass=mean(LnMass, na.rm=TRUE), Locomotor.Trait=unique(Locomotor.Trait)[1],
                                         TraitValue=mean(TraitValue, na.rm=TRUE), TempTraitValue=mean(TempTraitValue, na.rm=TRUE), Temperature=mean(Temperature, na.rm=TRUE),
                                         TempType=unique(TempType)[1], Environment=unique(Environment)[1], Thermy=unique(Thermy), TrophicGroup=unique(TrophicGroup)[1])
  tds[[i]] <- try(make.treedata(ptree, .data))
  if(class(tds[[i]])=="try-error"){
    tds[[i]] <- .data
  }
}

saveRDS(tds, file="../output/fullTreeTreeDataObjects.rds")




############
times <- c(max(branching.times(ttolData$phy))-1, 2000, 1000, 750, 600, 500, 400, 350, 300,275, 250,225, 200,175, 150, 125, 100, 75, 50, 40, 30, 20, 10, 5, 2, 1)
phyndrTree <- phyndrTTOL(ttolData, otts$ott_id, times, ncores=50, ottids=TRUE, prune=FALSE)

phyndr_sample(phyndrTree$phyndr)

td <- make.treedata(ttolData$phy, locomotor)

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

