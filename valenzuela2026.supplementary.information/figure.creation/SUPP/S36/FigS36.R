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
#primates_tree <- read.newick("../../data/Tree/FINAL_TREE_BL.newick")
#primates_tree <- root(primates_tree, outgroup = "Mus_musculus")

#A) READ LAST TRAITS FILE and metadata

#all_traits <- read.csv("../../data/Phenome/Traits_dataset_JUL2021.tsv", header=T, sep = "\t")

#Now let's bind the Group info


#B) READ CONFIG_FILES

#READ ALL CAAS PHYLOPS
#Here your path for files
path <- "PHYLOP_files/caas_phyloP/"


#nrow(SCENARIO1)
#nrow(SCENARIO2)

ALL_JOINED <- readRDS("PHYLOP_files/ALL_JOINED")
SCENARIO1 <- readRDS("PHYLOP_files/SCENARIO1")
SCENARIO2 <- readRDS("PHYLOP_files/SCENARIO2")
Sim <- readRDS("PHYLOP_files/Sim")

#ALL_PHYLOPs <- read.csv("PHYLOP_files/all_phyloP/ALL_CDS_PHYLOP.out", sep ="\t", header=FALSE)

ALL_JOINED_FILT <- distinct(ALL_JOINED[,which(colnames(ALL_JOINED) %in% c("Chrom", "Coord", "PhyloP", "strand"))])
SCENARIO1_FILT <- distinct(SCENARIO1[,which(colnames(SCENARIO1) %in% c("Chrom", "Coord", "PhyloP", "strand"))])
SCENARIO2_FILT <- distinct(SCENARIO2[,which(colnames(SCENARIO2) %in% c("Chrom", "Coord", "PhyloP", "strand"))])
SIM_FILT <- Sim[,c(1,2,4,3)]
colnames(SIM_FILT) <- c("Chrom", "Coord", "strand", "PhyloP")

#head(SCENARIO2_FILT)



ALL <- dplyr::bind_rows(list(CAAS=ALL_JOINED_FILT,
                             CDS=SIM_FILT,
                             SCENARIO1=SCENARIO1_FILT, 
                             SCENARIO2=SCENARIO2_FILT), .id = 'source')

#Transvar1 <- read.csv("../Tests_PhyloP_Coord_FEB22/zoonomia_transvar/ALL_PHYLOPs_TRANSVAR_scenario1.txt", header=FALSE, sep="\t")
#Transvar2 <- read.csv("../Tests_PhyloP_Coord_FEB22/zoonomia_transvar/ALL_PHYLOPs_TRANSVAR_scenario2.txt", header=FALSE, sep="\t")
#BOTH <- read.csv("../Tests_PhyloP_Coord_FEB22/zoonomia_transvar/ALL_PHYLOPs_TRANSVAR_BOTH.txt", header=FALSE, sep="\t")

#ALL <- dplyr::bind_rows(list(CAAS=BOTH, CDS=Sim,
 #                            SCENARIO1=Transvar1, SCENARIO2=Transvar2), .id = 'source')
#ALL <- dplyr::bind_rows(list(CAAS=Transvar, Background=Transvar), .id = 'source')

ALL$source[ALL$source == "CDS"] <- "1e+7 Background CDS" 
ALL$source[ALL$source == "SCENARIO1"] <- "Scenario1" 
ALL$source[ALL$source == "SCENARIO2"] <- "Scenario2" 

summary <- ALL%>% 
  group_by(source) %>% 
  tally()




library(ggpubr)
library(rstatix)


stat.test <- ALL %>% t_test(PhyloP ~ source) %>%
  adjust_pvalue(method = "fdr")# %>%
 # add_significance("p.adj.signif")
stat.test <- stat.test %>% add_xy_position(x = "source")

bxp <- ggplot(ALL, aes(x=source, y=PhyloP, fill=source)) + 
  geom_boxplot() +
  geom_text(data = summary,
            aes(source, Inf, label = paste("N= ",n, sep="")), vjust = 1) + theme_minimal() +
  xlab("Dataset") + ylab("Phylo-P score")

ALL$source<- factor(ALL$source , levels=c("1e+7 Background CDS", 
                                           "CAAS", 
                                           "Scenario1", 
                                           "Scenario2"))

bxp <- ggboxplot(ALL, x = "source", y = "PhyloP", fill = "source", 
                 palette = c("#00AFBB", "#E7B800", "#FC4E07", "#23764A")) +
  ylab("Phylo-P score") +
  #geom_text(data = summary,
   #         aes(source, Inf, label = paste("N= ",n, sep="")), vjust = 1)
  theme(axis.text.x = element_blank(), 
        axis.ticks.x =  element_blank(),
        axis.title.x = element_blank())
p2 <- bxp + stat_pvalue_manual(stat.test, label = "p.adj.signif", 
                         y.position = c(18, 12, 20, 10, 16, 14),
                         bracket.shorten = 0.05)

  #geom_signif(comparisons = list(c("SCENARIO1", "Background"),
        #                         c("SCENARIO2", "Background"),
         #                        c("SCENARIO2", "SCENARIO1")),
          #    test = "wilcox.test",
            #  map_signif_level = TRUE) + theme_minimal()

saveRDS(p2, "FigureS36")


