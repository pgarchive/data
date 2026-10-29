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
added_BG <- read.csv("ALLOMETRIC_EXTREMES", header = T, sep="\t")


#ggtree plot spreadness

tree_y <-  function(ggtree, data){
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  left_join(dplyr::select(data, label), dplyr::select(ggtree$data, label, y)) %>%
    pull(y)
}

library(ggtree)
primates_tree <- read.newick("../../Figures/Fig1/FINAL_TREE_BL_WITHOUTNODELAB.newick")

#setdiff(added_BG_filt$label, primates_tree$tip.label)
  #  View(unique(primates_tree$tip.label))    
p <- ggtree(primates_tree)

families <- c()
for(i in 1:length(p$data$label)){
  if (as.character(p$data$label[i]) %in% as.character(sequenced_traits$label)){
    family <- as.character(sequenced_traits$Family[sequenced_traits$label == as.character(p$data$label[i])])
    families <- c(families, family)}
  else{
    family <- NA
    families <- c(families, family)
  }}

p$data$Family <- families
library(randomcoloR)
library(ggtree)

#PLOT SPREADNESS SPECIES

primate_families_palette<- readRDS("../../Figures/Fig1/primate_families_palette.rds")
p1 <-  p + geom_tippoint(aes(color=Family), size=1.5,
                         show.legend = FALSE) +
  scale_color_manual(values = primate_families_palette) 

#View(Traits_top[which(Traits_top$variable == "LitterSz"),])
#View(final_ALL[which(final_ALL$variable == "LitterSz"),])

#unique(Traits_top$variable)
added_BG_filt <- added_BG[which(!added_BG$variable %in% c("BrainMass_g",
                                                          "meanMaleBodyMass.in.g")),]
p2 <- ggplot(data = added_BG_filt) + geom_point(shape = 22, size = 4,
                                                mapping = aes(x = as.numeric(as.character(value)),
                                                              y= tree_y(p1, added_BG_filt),
                                                              color = factor(Family.x),
                                                              fill = dplyr::case_when(!is.na(condition_FG) ~ "darkgreen", 
                                                                                      !is.na(condition_BG) ~ "darkred",
                                                                                      TRUE ~ "black")), show.legend = F) +
  #geom_text( mapping = aes(x = as.numeric(as.character(value)),
  #                        y= tree_y(p1, added_BG_filt),
  #                       label = label),
  #        size=3) +
  
  scale_fill_manual(values = c("white", "darkgreen", "darkred")) +
  scale_color_manual(values = primate_families_palette) +
  facet_wrap(~variable, scales = "free") +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())

saveRDS(p2, "FigureS25")
