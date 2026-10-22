######################################
##COMPONENTS STUDY FOR TRAITS GROUPS##
###BASED ON PRESENCE/ABSCENSE, TO ####
###DETERMINE CONSISTENCY FROM DIFFER##
## # SOURCE DATASETS AND KNOW THE TRAITS#
##WE CAN COMPARE BETWEEN EACH OTHER DUE TO###
#THE SPECIES THEY CONBTAIN ###########
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
library(dendextend)
library(ape)
library(randomcoloR)
library(tidyverse)
library(magrittr)
library(cluster)
library(cluster.datasets)
library(cowplot)
library(NbClust)
library(clValid)
library(ggfortify)
library(clustree)
library(dendextend)
library(factoextra)
library(FactoMineR)
library(corrplot)


#A) READ LAST TRAITS FILE and metadata

all_traits <- read.csv("../../Figures/Fig1//PRIMATES_PHENOMIC_DATASET.tsv", header=T, sep = "\t")





#B) SELECT DESIRED SPECIES (sequenced_ones)
#I think it's better to interpret the trait from the whole set of trait data , then
#we can select in the tree those species which have been sequenced


#B) SELECT DESIRED SPECIES (sequenced_ones)

primates_palette <- readRDS("../../Figures/Fig1/primate_families_palette.rds")
seq_species <- read.csv("../../Figures/Fig1/names_FINAL.txt",
                        sep = "\t", header = TRUE)
colnames(seq_species)<- c("Tree_Species", "Genome_species", "Trait_species")
covered_species_all <- as.character(seq_species$Trait_species[!is.na(seq_species$Tree_Species)])
covered_species <- as.character(covered_species_all[!is.na(covered_species_all)])
sequenced_traits <- all_traits[which(all_traits$SpeciesBROAD %in% covered_species),]
levels(sequenced_traits$SpeciesBROAD) <- c(levels(sequenced_traits$SpeciesBROAD), "Tarsius_bancanus", "Leontocebus_fuscicolis",
                                           "Plecturocebus_hoffmansi", "Lagothrix_lagothricha") 
sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Cephalopachus_bancanus"] <- "Tarsius_bancanus"
sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Leontocebus_fuscicollis"] <- "Leontocebus_fuscicolis"
sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Plecturocebus_hoffmannsi"] <- "Plecturocebus_hoffmansi"
sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Lagothrix_lagotricha"] <- "Lagothrix_lagothricha"

groups_list <- strsplit(as.character(sequenced_traits$GroupName), "_")
groups <- c()
for(i in 1:length(groups_list)){
  group <- paste(groups_list[[i]][1], "_", groups_list[[i]][2], sep="")
  groups[i] <- group}
groups
sequenced_traits <- cbind(sequenced_traits, groups)

nrow(sequenced_traits)
ncol(sequenced_traits)
#write.table(as.data.frame(cbind(as.character(sequenced_traits$SpeciesBROAD), 
 #                               as.character(sequenced_traits$groups))),
  #          "conversion_table.txt", quote=FALSE,  col.names = TRUE, row.names=FALSE, sep="\t")

#C) SELECT TRAITS OF INTEREST

#binary_traits <- read.csv("~/Documents/LInux_DIC_2020/Documents/PhD_EvoGenomics/3rd_year/GENOME-PHENOME-PAPER_NOV2021/results/Metadata/Metadata_binary_traits_FINAL.tsv", sep = "\t")
#binary_traits$nature <- "binary"

#categorical_traits <- read.csv("~/Documents/LInux_DIC_2020/Documents/PhD_EvoGenomics/3rd_year/GENOME-PHENOME-PAPER_NOV2021/results/Metadata/Metadata_categorical_traits_FINAL.tsv", sep = "\t")
#categorical_traits$nature <- "categorical"

#numeric_traits <- read.csv("~/Documents/LInux_DIC_2020/Documents/PhD_EvoGenomics/3rd_year/GENOME-PHENOME-PAPER_NOV2021/results/Metadata/Metadata_continous_traits_FINAL.tsv", sep = "\t")
#numeric_traits$nature <- "numeric"

#integer_traits <- read.csv("~/Documents/LInux_DIC_2020/Documents/PhD_EvoGenomics/3rd_year/GENOME-PHENOME-PAPER_NOV2021/results/Metadata/Metadata_discrete_traits_FINAL.tsv", sep = "\t")
#integer_traits$nature <- "integer"

#build global dataset (filter non-informative categories)

domains <- read.csv("../../Figures/Fig1/Traits_domains_sorted_FIG1.txt",
                    sep="\t")
#selected_traits <- domains$Measurement[!domains$Domain_Secondary_Domain %in% c("Parasite-Host_Health", "Morphology_Geometric Morphometrics")]


#Trait_domains <- c(as.character(binary_traits$Trait_domain),
 #                  as.character(categorical_traits$Trait_domain),
  #                 as.character(numeric_traits$Trait_domain))
#Measurements <- c(as.character(binary_traits$Measurement),
 #                 as.character(categorical_traits$Measurement),
  #                as.character(numeric_traits$Measurement))
#Natures <- c(as.character(binary_traits$nature),
 #            as.character(categorical_traits$nature),
  #           as.character(numeric_traits$nature))
#all_metadata <- as.data.frame(cbind(Trait_domains,
 #                             Measurements,
  #                            Natures))

#D) BINARIZE THE MATRIX

#before calling MCA algorithm for factor analysis of traitsxspecies,
#let's transform the data into presence/absence of traits

binarized_traits <- sequenced_traits[,which(!colnames(sequenced_traits) %in% c("SpeciesBROAD",  "GroupName"))] %>%
  mutate_all(funs(ifelse(is.na(.), 0, 1)))

#binarized_traits <- binarized_traits[,which(colnames(binarized_traits) %in% all_metadata$Measurements)]

#E) STUDY THE CLUSTERS OF HOMOGENEITY AND CONSISTENCY

binarized_traits <- binarized_traits %>%
  mutate_all(as.character)
rownames(binarized_traits) <- sequenced_traits[,1]
binarized_res.mca <- MCA(binarized_traits[,-ncol(binarized_traits)], 
                         graph=FALSE, ncp = 10)

num_matrix <-  mutate_all(as.data.frame(t(binarized_traits[,-ncol(binarized_traits)])), 
                          function(x) as.numeric(as.character(x)))


trans_binarized_res.mca <- MCA(t(binarized_traits[,-ncol(binarized_traits)]), 
                               graph=FALSE, ncp = 13)
trans_res.HCPC <- HCPC(trans_binarized_res.mca, nb.clust = 13, max = NULL, graph = TRUE)
fviz_cluster(trans_res.HCPC,
             show.clust.cent = TRUE, # Show cluster centers
             palette = distinctColorPalette(20),         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             geom = "point",
             main = "Factor map"
)



#H) FIND THE SPECIES CLUSTERS
# DECIDE THE OPTIMAL CLUSTERING METHOD NUMBER (ELBOW PROCESS) or (GAP STATISTIC) for species

res.HCPC <- HCPC(binarized_res.mca, nb.clust = -1, min = 3, max = NULL, graph = FALSE)
p2 <- fviz_cluster(res.HCPC,
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco", # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map",
             geom = c("point")
             
)

saveRDS(p2, "FigureS3")





