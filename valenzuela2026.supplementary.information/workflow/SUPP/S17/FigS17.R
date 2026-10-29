######################################
##COMPONENTS STUDY FOR TRAITS GROUPS##
###, BASED ON MERGED MULTIPLE GROUPS##
###IN THE ANALYSIS OF CATEGORIES #####
###AND VARIABLES #####################
######################################

library(ggplot2)
library(gridExtra)
library(tibble)
library(data.table)
library(ggrepel)
library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)
library(dendextend)
library(ComplexHeatmap)
library("factoextra")
library("circlize")
library("RColorBrewer")
library(ggnewscale)
library(scales)
library(cowplot)
library(ggdendro)
library(viridis)
library(cowplot)
library(FactoMineR)
library(missMDA)
library(psych)
library(network)
library(ape)
library(ggtree)
library(randomcoloR)
library(devtools)
library(phytools)
library(fastDummies)
options(scipen=0)



#A) READ LAST TRAITS FILE and metadata

all_traits <- read.csv("../../Figures/Fig1/PRIMATES_PHENOMIC_DATASET.tsv", header=T, sep = "\t")

#Now let's bind the Group info



#B) SELECT DESIRED SPECIES (sequenced_ones)


primates_palette <- readRDS("../../Figures/Fig1/primate_families_palette.rds")
seq_species <- read.csv("../../Figures/Fig1/names_FINAL.txt",
                        sep = "\t", header = TRUE)
colnames(seq_species)<- c("Tree_Species", "Genome_species", "Trait_species")
covered_species_all <- as.character(seq_species$Trait_species[!is.na(seq_species$Tree_Species)])
covered_species <- as.character(covered_species_all[!is.na(covered_species_all)])
sequenced_traits <- all_traits[which(all_traits$SpeciesBROAD %in% covered_species),]

groups_list <- strsplit(as.character(sequenced_traits$GroupName), "_")
groups <- c()
for(i in 1:length(groups_list)){
  group <- paste(groups_list[[i]][1], "_", groups_list[[i]][2], sep="")
  groups[i] <- group}
groups
sequenced_traits <- cbind(sequenced_traits, groups)
sequenced_traits <- sequenced_traits %>%
  mutate_at(vars(`IUCN_ER`), funs(case_when(
    . %in% c("5", "4", "3") ~ "1", 
    . %in% c("2", "1") ~ "0")))


#read redundant measurements

redundant <- read.csv("Redundant_measurements.tsv", header=FALSE, sep = "\t")

#D) BINARIZE THE MATRIX

binarized_traits <- sequenced_traits[,-which(colnames(sequenced_traits) %in% c("SpeciesBROAD",  "GroupName",
                                                                               "groups") |
                                               colnames(sequenced_traits) %in% redundant$V1)] %>%
  mutate_all(funs(ifelse(is.na(.), 0, 1)))


#E) STUDY THE CLUSTERS OF HOMOGENEITY AND CONSISTENCY

#binarized_traits <- cbind(binarized_traits, groups)

rownames(binarized_traits) <- sequenced_traits[,1]


##############################################
#A) ANALYSIS OF CORRELATIONS WITH CLUSTERS###
##############################################


coverage_traits_clustered <- read.csv("../../Figures/Fig1/Traits_clusters_FINAL_JAN23.tsv", sep = "\t")

traits_clust <- coverage_traits_clustered[which(coverage_traits_clustered$Cluster ==1),]

data_clust <- binarized_traits[,which(colnames(binarized_traits) %in% traits_clust$Trait)]

num_data_clust <- sapply(data_clust, as.numeric)
rownames(num_data_clust) <- rownames(data_clust)
##Subset as usual
best_species_clust <- rownames(num_data_clust[which(rowSums(num_data_clust) >= 0.90*ncol(num_data_clust)),])
best_species_clust
all_traits_clust <- binarized_traits[which(rownames(binarized_traits) %in% best_species_clust),]


#totals <- colSums(sapply(all_traits_clust, as.numeric))
best_traits_clust <- all_traits_clust[, which(colSums(all_traits_clust) >= 0.90*nrow(all_traits_clust))]
best_traits_clust <- colnames(best_traits_clust)

#MULTIVARIATE ANALYSIS
trues <- sequenced_traits$SpeciesBROAD %in% best_species_clust
#length(trues[trues==TRUE])

selected_subset <- sequenced_traits[which(sequenced_traits$SpeciesBROAD %in% best_species_clust),
                                    which(colnames(sequenced_traits) %in% best_traits_clust)] 

rownames(selected_subset) <- sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD %in% best_species_clust]
selected_subset <- selected_subset[,which(!colnames(selected_subset) %in% as.character(redundant$V1))]



library(missMDA)
no_missing <- selected_subset %>% drop_na()

res.famd <- FAMD(no_missing) #cluster6 | cluster7 | cluster 8

nrow(selected_subset)
#sequenced_traits[,1]
fviz_famd_var(res.famd, "quali.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)
options(ggrepel.max.overlaps = Inf)

p1 <- fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)



#2) make correlations for quantitative sets of data

var_cor <- as.data.frame(res.famd$quanti.var)

primates_tree <- read.newick("../../Figures/Fig1/FINAL_TREE_BL_WITHOUTNODELAB.newick")

#CLUSTS
quanti_clust <- selected_subset[,which(colnames(selected_subset) %in% rownames(var_cor)
                                         & colnames(selected_subset) != "IUCN_ER"
                                       & colnames(selected_subset) != "R_Min_LitterSize")] %>% drop_na()
matched_rows <- match(rownames(quanti_clust),
                      primates_tree$tip.label)
matched_rows <- matched_rows[!is.na(matched_rows)]
clust.tree<-drop.tip(primates_tree,primates_tree$tip.label[-matched_rows])
quanti_clust_phylo <- quanti_clust[which(rownames(quanti_clust) %in% primates_tree$tip.label[matched_rows] ),] 

#imputed_DAT <- imputeFAMD(quanti_clust_phylo)
correl_struct <- phyl.vcv(as.matrix(quanti_clust_phylo),vcv(clust.tree),1)
correl <- cov2cor(correl_struct$R)
#correl_struct <- phyl.vcv(as.matrix(imputed_DAT$completeObs),vcv(clust.tree),1)

correl <- correl[rowSums(is.na(correl)) != ncol(correl)-1,]
correl <- correl[,colSums(is.na(correl)) != nrow(correl)]



t.xy<-correl*sqrt((Ntip(clust.tree)-2)/(1-correl^2))
#View(t.xy)
P.xy<-2*pt(abs(t.xy),df=Ntip(clust.tree)-2,lower.tail=F)
P.xy

MELTED <- melt(P.xy)
MELTED$fdr <- p.adjust(MELTED$value)
FDR <- MELTED[,c(1,2,4)]
FDR_MATRIX <- acast(FDR, Var1~Var2, value.var='fdr', fill=0)


library(GGally)
#Correlation matrices

#View(correl)
p1 <- ggcorr(data = NULL, 
       cor_matrix = correl,
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 2,
       label_round = 2,
       color = "black",
       low = "darkred",
       mid = "#EEEEEE", high = "darkgreen",
       layout.exp = 0.4,
       size=3) 

saveRDS(p1, "FigureS17")
