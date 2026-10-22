########################################################
##DIFFERENT METHODS TO DETECT PRIOR GROUPS FOR CAAS#####
########################################################

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
library(network)
library(dendextend)
library(ape)
library(randomcoloR)
library(phytools)
library(scales)
library(stringr)

#READ TREE
primates_tree <- read.newick("../../Figures/Fig1/FINAL_TREE_BL_WITHOUTNODELAB.newick")
primates_tree <- root(primates_tree, outgroup = "Mus_musculus")




#A) READ LAST TRAITS FILE and metadata

all_traits <- read.csv("../../Figures/Fig1/PRIMATES_PHENOMIC_DATASET.tsv", header=T, sep = "\t")

#Now let's bind the Group info


#B) READ CONFIG_FILES

#READ ALL SPECIES FILE
#Here your path for files
path <- "Configs_raw"

#get species names
list.configs <- list.files(path)
fields <- strsplit(list.configs, "/", fixed = TRUE)
genes <- c()

ldf <- list()
count <- 1
for (k in 1:(length(list.configs))){
  print(k)
  for (j in k:(length(list.configs))){
    print(j)
  # read pair of trait files
  trait1 <- list.configs[k]
  trait2 <- list.configs[j]
  trait_name1 <- strsplit(trait1, "\\.")[[1]][1]
  trait_name2 <- strsplit(trait2, "\\.")[[1]][1]
 filepath1 <- paste(path,  trait1, sep="/")
 filepath2 <- paste(path, trait2, sep="/")
  config_df1  <- read.csv(filepath1, sep ="\t", na.strings = c('NaN', ''), header = FALSE)
  config_df2  <- read.csv(filepath2, sep ="\t", na.strings = c('NaN', ''), header = FALSE)
 #compare tops
  top1 <- config_df1$V1[config_df1$V2 == "1"] 
  top2 <- config_df2$V1[config_df2$V2 == "1"] 
  setA <- setdiff(top1, top2)
  setB <- setdiff(top2, top1)
  nonoverlap <- unique(c(setA, setB))
  overlap <- intersect(top1, top2)
  total <- c(nonoverlap, overlap)
  distance_top <- 1- (length(overlap)/(length(overlap) + length(nonoverlap)))
  #compare bottoms
 bottom1 <- config_df1$V1[config_df1$V2 == "0"] 
  bottom2 <- config_df2$V1[config_df2$V2 == "0"] 
  setC <- setdiff(bottom1, bottom2)
  setD <- setdiff(bottom2, bottom1)
  nonoverlap_bot <- unique(c(setC, setD))
  overlap_bot <- intersect(bottom1, bottom2)
  total_bot <- c(nonoverlap_bot, overlap_bot)
  distance_bot <- 1- (length(overlap_bot)/(length(overlap_bot) + length(nonoverlap_bot)))
  #merge in dataframe
  ldf[[count]] <- as.data.frame(cbind( trait_name1,
                                       trait_name2,
                                      distance_bot,
                                      distance_top))

   #trait_names <- rep(trait_name, nrow(config_df))
  
 # top <- config_df$V1[config_df$V2 == "1"] 
  #    ldf[[count]] <- as.data.frame(cbind(trait_names,
   #                                       config_df))
      count <- count + 1
  }}



distances_configs <- do.call(rbind, ldf)
distances_configs <- distances_configs[which(distances_configs$trait_name1 != distances_configs$trait_name2),]


#PLOT COMPARISONS

subset_top <- distances_configs[which(as.numeric(as.character( distances_configs$distance_bot)) < 0.3 |
                                        as.numeric(as.character(distances_configs$distance_top)) < 0.3), ]
Fig2 <- ggplot(data=distances_configs,aes(x=as.numeric(as.character(distance_bot)), 
                                  y=as.numeric(as.character(distance_top)))) + 
  geom_point(col="gold") +
  xlab("distance_bottom") + ylab("distance_top") +
  theme_minimal() 

saveRDS(Fig2, "FigureS27")
