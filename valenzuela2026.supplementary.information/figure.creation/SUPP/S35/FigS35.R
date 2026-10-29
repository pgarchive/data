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

#ggplot stuff


#Freq of covered codon positions

table_covered <- summary_stats_trait_independent %>% drop_na() 

p2 <- ggplot(data = table_covered) + geom_density(aes(x=Freq_top),
                                                                 fill="darkgreen") + theme_minimal() +
  xlab("Frequency AA TOP in mammalian gene alignments (Zoonomia)") + ylab("Total density")

#Freq top aa

table_covered <- summary_stats_trait_independent %>% drop_na() 

#freq orders


#X cuts

b <- seq(-0.9, 1, 0.1)
table_covered$Freq_top_cuts <- cut(as.numeric(as.character(table_covered$Freq_top)), 
                                         breaks = b)


#CROSSCHECK WITH PERM VALUES
pvalues_perm <- read.csv("../../Figures/Fig2/CAAS_FINAL_MAR23.tab", sep ="\t", header=FALSE)
colnames(pvalues_perm) <- c("Gene", "Trait", "Pos", "AA", "pvalue", "scenario", "Top_n",
                            "Bottom_n")

summary_stats_df_filt <- summary_stats_df_contr %>%
  drop_na()
total_df <- left_join(summary_stats_df_filt, pvalues_perm)

test <- total_df[which(total_df$pvalue<0.05),]
nrow(test)


b <- seq(-0.9, 1, 0.1)
total_df$Freq_top_cuts <- cut(as.numeric(as.character(total_df$Freq_top)), 
                                   breaks = b)



##################################
#COMBINATION WITH SCENARIOS
#PERMUL_RESULTS <- read.csv("ALL_RESULTS_PERMUL.tab", sep="\t", header = FALSE)
#colnames(PERMUL_RESULTS) <- c("Gene", "Trait", "Pos", "AA", "pvalue_perm", "scenario")

summary_stats_df_filt <- summary_stats_df_contr %>%
  drop_na()
total_df <- left_join(summary_stats_df_filt, pvalues_perm)
total_df_filt <- total_df %>%
  drop_na()
#X cuts

b <- seq(-0.9, 1, 0.1)
total_df_filt$Freq_top_cuts <- cut(as.numeric(as.character(total_df_filt$Freq_top)), 
                                   breaks = b)

p2 <- ggplot(total_df_filt[which(total_df_filt$pvalue < 0.05),],
       aes(x=as.character(Freq_top_cuts)))+ geom_histogram( aes(fill=scenario),
                                                            stat="count",
                                                            position = "dodge") +
  theme_minimal() +   xlab("Frequency AA TOP in mammalian gene alignments (Zoonomia)") + ylab("Total density") +
  scale_fill_manual(values = c("gold", "lightblue")) + theme(legend.title = element_blank(),
                                                             axis.text.x = element_text(angle=90))

#facet_wrap(~Trait, scales = "free") + scale_fill_manual(values = c("gold", "lightblue"))

saveRDS(p2, "FigureS35")

