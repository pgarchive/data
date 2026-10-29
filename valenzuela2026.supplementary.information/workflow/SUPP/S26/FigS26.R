########################################################
##DIFFERENT METHODS TO DETECT PRIOR GROUPS FOR CAAS#####
########################################################

library(rvcheck)
library(ggtree)
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
library(fastDummies)

all_traits <- read.csv("../../Figures/Fig1/PRIMATES_PHENOMIC_DATASET.tsv", header=T, sep = "\t")

all_categorical_traits_extremes_filt <- read.csv("Qualitative_extremes", sep="\t", header=F)
#ggtree plot spreadness

primates_tree <- read.newick("../../Figures/Fig1/FINAL_TREE_BL_WITHOUTNODELAB.newick")

primates_tree$tip.label
tree_y <-  function(ggtree, data){
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  left_join(dplyr::select(data, label), dplyr::select(ggtree$data, label, y)) %>%
    pull(y)
}

#library(ggtree)
p <- ggtree(primates_tree)

families <- c()
for(i in 1:length(p$data$label)){
  if (as.character(p$data$label[i]) %in% as.character(sequenced_traits$SpeciesBROAD)){
    family <- as.character(sequenced_traits$groups[sequenced_traits$SpeciesBROAD == as.character(p$data$label[i])])
    families <- c(families, family)}
  else{
    family <- NA
    families <- c(families, family)
  }}

p$data$Family <- families
#View(p$data)
library(randomcoloR)

#PLOT SPREADNESS SPECIES

primate_families_palette <- readRDS("../../Figures/Fig1/primate_families_palette.rds")
p1 <-  p + geom_tippoint(aes(color=Family), size=1.5,
                         show.legend = FALSE) +
  scale_color_manual(values = primate_families_palette) 


all_categorical_traits_extremes_filt <- readRDS("all_categorical_traits_extremes_2BG_FINAL")
all_categorical_traits_extremes_filt$label[all_categorical_traits_extremes_filt$label == "Lagothrix_lagothricha"] <- "Lagothrix_lagotricha"

#View(all_categorical_traits_extremes_filt)
p2 <- ggplot(data = all_categorical_traits_extremes_filt) + geom_point(shape = 21, size = 5,
                                                   mapping = aes(x = as.numeric(as.character(value)),
                                                                 y= tree_y(p1, all_categorical_traits_extremes_filt),
                                                                 color = factor(Family.x),
                                                                 fill = dplyr::case_when(extremes == "high_extreme" ~ "darkgreen", 
                                                                                         extremes == "low_extreme" ~ "darkred",
                                                                                         TRUE ~ "black")), show.legend = FALSE,
                                                   position=position_jitter(width = 0.1, height = 0.1)) +
  scale_fill_manual(values = c("white", "darkgreen", "darkred", "darkgreen", "darkred")) +
  scale_color_manual(values = primate_families_palette) +
  facet_wrap(~Trait, scales = "free") +
  theme(axis.title.y =element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()) 
saveRDS(p2, "FigureS26")  

 # geom_text(data = count_subset, 
 #           aes(x = 1.5,
  #              y= 230,
   #             fill = extremes,
    #            label = paste(n, "spc", "")),
     #       size = 3) +
  #scale_fill_manual(values = c("white", "darkgreen", "darkred", "darkgreen", "darkred")) +
  #scale_color_manual(values = primate_families_palette) +
  #facet_wrap(~Trait, scales = "free") +
  #xlim(0, 1.75)
  
