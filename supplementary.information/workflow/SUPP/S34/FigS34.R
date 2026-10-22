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


summary_stats_df <- read.csv("ALL_MAMMALIAN_INFO_MAY23_FILTTOT.txt", sep="\t", header = FALSE)
colnames(summary_stats_df) <- c("Gene", "Trait", "Pos", "AA", "Chrom",
                                "Coord", "Type", "Freq_top", "Freq_taxon", "Freq_distinct_aa")

summary_stats_df_contr <- summary_stats_df[which(!summary_stats_df$Trait %in% traits_with_no_contrast),]


summary_stats_trait_independent <- summary_stats_df_contr %>%
  distinct(Gene, Pos, AA, Chrom, Coord, Type, Freq_top, Freq_taxon, Freq_distinct_aa)

#Freq of covered codon positions

table_covered <- summary_stats_trait_independent %>% drop_na() 

p2 <- ggplot(data = table_covered) + geom_density(aes(x=Freq_top),
                                                                 fill="darkgreen") + theme_minimal() +
  xlab("Frequency AA TOP in mammalian gene alignments (Zoonomia)") + ylab("Total density")

saveRDS(p2, "FigureS34")
