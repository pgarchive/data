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
library('rphylopic')
library(ggimage)

#I) READ STUDY CASES CAASES

pvalues_perm <- read.csv("CAAS_global_and_altern.txt", sep ="\t", header=FALSE)
colnames(pvalues_perm) <- c("Gene", "Trait", "Pos", "AA", "pvalue", "scenario", "Top_n",
                            "Bottom_n")
pvalues_perm$pvalue[which(pvalues_perm$pvalue == 0)] <- 0.001
pvalues_fdr <- pvalues_perm %>%
  drop_na() %>%
  group_by(Trait) %>%
  mutate(pval_adj = p.adjust(pvalue, method='fdr'))
pvalues_perm_sign <- pvalues_perm[which(pvalues_perm$pvalue <0.05),]

#unique(pvalues_perm_sign$Trait)

#II) OBTAIN SELECTION BARPLOTS

library("ggvenn")
pvalues_perm_analysis <- pvalues_perm[which(as.numeric(as.character(pvalues_perm$pvalue)) <0.05),]
levels(pvalues_perm_analysis$Trait) <- c(levels(pvalues_perm_analysis$Trait),
                                         c("Maximum_lifespan__y__G",
                                           "Maximum_lifespan__y__A", "Maximum_lifespan__y__B",
                                           "TrophicGuild_Insectivore_G",
                                           "TrophicGuild_Insectivore_A",
                                           "TrophicGuild_Insectivore_B", 
                                           "White_Blood_Cells_G", "White_Blood_Cells_A", "White_Blood_Cells_B"))
pvalues_perm_analysis$Trait[pvalues_perm_analysis$Trait == "Maximum_lifespan__y__altern1_config.txt"] <- "Maximum_lifespan__y__A"
pvalues_perm_analysis$Trait[pvalues_perm_analysis$Trait == "Maximum_lifespan__y__altern2_config.txt"] <- "Maximum_lifespan__y__B"
pvalues_perm_analysis$Trait[pvalues_perm_analysis$Trait == "Maximum.lifespan..y."] <- "Maximum_lifespan__y__G"
pvalues_perm_analysis$Trait[pvalues_perm_analysis$Trait == "TrophicGuild_Insectivore_altern1_config.txt"] <- "TrophicGuild_Insectivore_A"
pvalues_perm_analysis$Trait[pvalues_perm_analysis$Trait == "TrophicGuild_Insectivore_altern2_config.txt"] <- "TrophicGuild_Insectivore_B"
pvalues_perm_analysis$Trait[pvalues_perm_analysis$Trait == "TrophicGuild_Insectivore"] <- "TrophicGuild_Insectivore_G"
pvalues_perm_analysis$Trait[pvalues_perm_analysis$Trait == "White_Blood_Cells_altern1_config.txt"] <- "White_Blood_Cells_A"
pvalues_perm_analysis$Trait[pvalues_perm_analysis$Trait == "White_Blood_Cells_altern2_config.txt"] <- "White_Blood_Cells_B"
pvalues_perm_analysis$Trait[pvalues_perm_analysis$Trait == "White.Blood.Cells"] <- "White_Blood_Cells_G"

pvalues_perm_analysis$Pos <- as.numeric(as.character(pvalues_perm_analysis$Pos))+1

#Read selection files
SELECTION <- read.csv("ALL_SELECTION_STUDY_CASES_FILT.txt", sep ="\t", header=FALSE)
colnames(SELECTION) <- c("Trait", "Gene", "Pos",	"Partition", "target", "rate", 
                         "bias", "Bayes Factor", "Amino acid composition of site",	
"Substitution history on selected branches")

#SELECTION %>% drop_na() %>% count(Trait)
#View(head(SELECTION))

test <- SELECTION %>% drop_na()
nrow(SEL_PLOT[which(!is.na(SEL_PLOT$Selection)),])
nrow(test)

#JOIN
JOINED <- SELECTION %>% drop_na() %>%
  full_join(pvalues_perm_analysis, SELECTION,
          by=c("Gene" = "Gene",
               "Trait" = "Trait",
               "Pos" = "Pos"))

#View(JOINED)
#View(JOINED)

#TABLE_TO_RETRIEVE_OUNT
SEL_PLOT <- JOINED %>% 
  mutate(Selection = case_when(!is.na(pvalue) & !is.na(target) ~ "FADE positive & CAAS",
                               is.na(pvalue) & !is.na(target) ~ "FADE positive",
                               TRUE ~ NA
  )) %>% group_by(Trait, Selection) %>%
  count(Trait, Selection) %>%
  drop_na()

#View(SEL_PLOT)


#PATTERNS
#SEL_PLOT 

Max_Lifespan <- SEL_PLOT[which(SEL_PLOT$Trait %in% c(
  "Maximum_lifespan__y__G", "Maximum_lifespan__y__A",
  "Maximum_lifespan__y__B"
)),]
Max_Lifespan$x <- c("A", "A", "B", "B", "S")
traits_palette <- readRDS("traits_palette")
colfunc_ANC <- colorRampPalette(c(traits_palette[3], "white"))

library(ggpattern)
ggSEL2 <- ggplot(Max_Lifespan) + 
  geom_bar_pattern(aes(x = factor(x, level=c('S', 'A', 'B')),
                       y=n, fill=factor(x, level=c('S', 'A', 'B')),
                       pattern = Selection),
                   stat="identity")  +
 # theme_classic() +
  scale_fill_manual(values = 
                      c(colfunc_ANC(10)[c(9)], 
                        colfunc_ANC(10)[c(6)],
                        colfunc_ANC(10)[c(1)]))+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank()) +
  scale_pattern_manual(values = c("FADE positive & CAAS" = "stripe",
                                  "FADE positive" = "none")) +
  ylab("FADE positive") +
  xlab("CAAS selection")


saveRDS(ggSEL2, "ggSEL2_FEB24")

#Insect
Insect <- SEL_PLOT[which(SEL_PLOT$Trait %in% c(
  "TrophicGuild_Insectivore_A", "TrophicGuild_Insectivore_B",
  "TrophicGuild_Insectivore_G"
)),]
Insect$x <- c("A", "A", "B", "B", "S", "S")
traits_palette <- readRDS("traits_palette")
colfunc_ANC <- colorRampPalette(c(traits_palette[1], "white"))


ggSEL1 <- ggplot(Insect) + 
  geom_bar_pattern(aes(x = factor(x, level=c('S', 'A', 'B')),
                       y=n, fill=factor(x, level=c('S', 'A', 'B')),
                       pattern = Selection),
                   stat="identity")  +
  #theme_classic() +
  scale_fill_manual(values = 
                      c(colfunc_ANC(10)[c(9)], 
                        colfunc_ANC(10)[c(6)],
                        colfunc_ANC(10)[c(1)]))+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank()) +
  scale_pattern_manual(values = c("FADE positive & CAAS" = "stripe",
                                  "FADE positive" = "none")) +
  ylab("FADE positive") +
  xlab("CAAS selection")


saveRDS(ggSEL1, "ggSEL1_FEB24")

#View(JOINED)
#White Blood Cells
WBC <- SEL_PLOT[which(SEL_PLOT$Trait %in% c(
  "White_Blood_Cells_A", "White_Blood_Cells_B",
  "White_Blood_Cells_G"
)),]
WBC$x <- c("A", "A", "B", "B", "S", "S")
traits_palette <- readRDS("traits_palette")
colfunc_ANC <- colorRampPalette(c("grey", "white"))


ggSEL3 <- ggplot(WBC) + 
  geom_bar_pattern(aes(x = factor(x, level=c('S', 'A', 'B')),
                       y=n, fill=factor(x, level=c('S', 'A', 'B')),
                       pattern = Selection),
                   stat="identity")  +
 # theme_classic() +
  scale_fill_manual(values = 
                      c(colfunc_ANC(10)[c(9)], 
                        colfunc_ANC(10)[c(6)],
                        colfunc_ANC(10)[c(1)]))+
  theme(legend.position = "none",
        panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.background = element_blank()) +
  scale_pattern_manual(values = c("FADE positive & CAAS" = "stripe",
                                  "FADE positive" = "none")) +
  ylab("FADE positive") +
  xlab("CAAS selection")


saveRDS(ggSEL3, "ggSEL3_FEB24")




"""
OVERLAP_INSECT <- SEL_PLOT[which(SEL_PLOT$Selection == "FADE positive & CAAS" &
(SEL_PLOT$Trait == "TrophicGuild_Insectivore_G" |
                 SEL_PLOT$Trait == "TrophicGuild_Insectivore_B" |
                 SEL_PLOT$Trait == "TrophicGuild_Insectivore_A")),]
OVERLAP_INSECT$GENEPOS <- paste(OVERLAP_INSECT$Gene, OVERLAP_INSECT$Pos)

#length(unique(OVERLAP_INSECT$Gene))
#length(unique(OVERLAP_INSECT$GENEPOS))

#INSECTIVORE OVERLAP WITH SELECTION
INSECTIVORE <- pvalues_perm[which(as.character(pvalues_perm$Trait) == "TrophicGuild_Insectivore" |
                                   as.character(pvalues_perm$Trait) == "TrophicGuild_Insectivore_altern1_config.txt" |
                                   as.character(pvalues_perm$Trait) == "TrophicGuild_Insectivore_altern2_config.txt" ),]

pvalues_perm_analysis <- INSECTIVORE[which(as.numeric(as.character(INSECTIVORE$pvalue)) <0.05  ),] 
pvalues_perm_analysis$GENEPOS <- paste(pvalues_perm_analysis$Gene, pvalues_perm_analysis$Pos)

SEL_DATA <- SELECTION %>% drop_na()
SELECTION_INSECT <- SEL_DATA[which(SEL_DATA$Trait == "TrophicGuild_Insectivore_G" |
                                     SEL_DATA$Trait == "TrophicGuild_Insectivore_B" |
                                     SEL_DATA$Trait == "TrophicGuild_Insectivore_A"),]
SELECTION_INSECT$GENEPOS <- paste(SELECTION_INSECT$Gene, SELECTION_INSECT$Pos)

#VISUALIZE OVERLAP INSECT

pvalues_perm_analysis_GENES <- pvalues_perm_analysis[,c(1,2)] %>%
  distinct()

JOINED_INSECT <- pvalues_perm_analysis_GENES[which(pvalues_perm_analysis_GENES$Gene %in% SELECTION_INSECT$Gene),] 
x_gene <- list(
  "S-CAAS" = unique(JOINED_INSECT$Gene[as.character(JOINED_INSECT$Trait) == "TrophicGuild_Insectivore"]),
  "A-CAAS" = unique(JOINED_INSECT$Gene[as.character(JOINED_INSECT$Trait) == "TrophicGuild_Insectivore_altern1_config.txt"]), 
  "B-CAAS" = unique(JOINED_INSECT$Gene[as.character(JOINED_INSECT$Trait) == "TrophicGuild_Insectivore_altern2_config.txt"])
)
names(x_gene) <- c("S-CAAS", "A-CAAS", "B-CAAS")
traits_palette <- readRDS("traits_palette")
colfunc_ANC <- colorRampPalette(c(traits_palette[1], "white"))
ggvenn(
  x_gene, 
  fill_color = c(colfunc_ANC(10)[c(9)], colfunc_ANC(10)[c(6)], colfunc_ANC(10)[c(1)]),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = TRUE
)

#TEST INTERSECT
genes_intersecting <- intersect(unique(SELECTION_INSECT$Gene),unique(pvalues_perm_analysis$Gene))
intersection <- length(intersect(unique(SELECTION_INSECT$Gene),unique(pvalues_perm_analysis$Gene)))
background <- 16120 #ALLGENES, SELECTION RUN ON ALL ALIGNMENTS (FLEXIBILITY OF SPECIES)
phyper(intersection-1, length(unique(SELECTION_INSECT$Gene)), 
       background-length(unique(SELECTION_INSECT$Gene)), length(unique(pvalues_perm_analysis$Gene)),
       lower.tail=F)

OVERLAP_INSECT$GENEPOS <- paste(OVERLAP_INSECT$Gene, OVERLAP_INSECT$Pos)

length(unique(OVERLAP_INSECT$Gene))
length(unique(OVERLAP_INSECT$GENEPOS))

#################################
#LIFESPAN OVERLAP WITH SELECTION#
#################################

pvalues_perm_analysis <- LIFESPAN
pvalues_perm_analysis$GENEPOS <- paste(pvalues_perm_analysis$Gene, pvalues_perm_analysis$Pos)

SEL_DATA <- SELECTION %>% drop_na()
#View(SELECTION_LIFESPAN)
SELECTION_LIFESPAN <- SEL_DATA[which(SEL_DATA$Trait == "Maximum_lifespan__y__B" |
                                     SEL_DATA$Trait == "Maximum_lifespan__y__A" |
                                     SEL_DATA$Trait == "Maximum_lifespan__y__G"),]
SELECTION_INSECT$GENEPOS <- paste(SELECTION_INSECT$Gene, SELECTION_INSECT$Pos)
#TEST OVERLAP GENES
intersection <- length(intersect(unique(SELECTION_INSECT$Gene),unique(pvalues_perm_analysis$Gene)))
background <- 16120 #ALLGENES, SELECTION RUN ON ALL ALIGNMENTS (FLEXIBILITY OF SPECIES)
phyper(intersection-1, length(unique(SELECTION_INSECT$Gene)), 
       background-length(unique(SELECTION_INSECT$Gene)), length(unique(pvalues_perm_analysis$Gene)),
       lower.tail=F)

#POS
intersection <- length(intersect(unique(SELECTION_INSECT$GENEPOS),unique(pvalues_perm_analysis$GENEPOS)))


OVERLAP_INSECT$GENEPOS <- paste(OVERLAP_INSECT$Gene, OVERLAP_INSECT$Pos)

length(unique(OVERLAP_INSECT$Gene))
length(unique(OVERLAP_INSECT$GENEPOS))


#################################
#WBC OVERLAP WITH SELECTION#
#################################

pvalues_perm_analysis <- WBC_SUBSET
pvalues_perm_analysis$GENEPOS <- paste(pvalues_perm_analysis$Gene, pvalues_perm_analysis$Pos)

SEL_DATA <- SELECTION %>% drop_na()
#View(SELECTION_LIFESPAN)
unique(SEL_DATA$Trait)
SELECTION_INSECT <- SEL_DATA[which(SEL_DATA$Trait == "White_Blood_Cells_B" |
                                     SEL_DATA$Trait == "White_Blood_Cells_A" |
                                     SEL_DATA$Trait == "White_Blood_Cells_G"),]
SELECTION_INSECT$GENEPOS <- paste(SELECTION_INSECT$Gene, SELECTION_INSECT$Pos)

#VISUALIZE OVERLAP WBC

pvalues_perm_analysis_GENES <- pvalues_perm_analysis[,c(1,2)] %>%
  distinct()

JOINED_WBC <- pvalues_perm_analysis_GENES[which(pvalues_perm_analysis_GENES$Gene %in% SELECTION_INSECT$Gene),] 
x_gene <- list(
  "S-CAAS" = unique(JOINED_WBC$Gene[as.character(JOINED_WBC$Trait) == "White.Blood.Cells"]),
  "A-CAAS" = unique(JOINED_WBC$Gene[as.character(JOINED_WBC$Trait) == "White_Blood_Cells_altern1_config.txt"]), 
  "B-CAAS" = unique(JOINED_WBC$Gene[as.character(JOINED_WBC$Trait) == "White_Blood_Cells_altern2_config.txt"])
)
names(x_gene) <- c("S-CAAS", "A-CAAS", "B-CAAS")
traits_palette <- readRDS("traits_palette")
colfunc_ANC <- colorRampPalette(c("grey", "white"))
ggvenn(
  x_gene, 
  fill_color = c(colfunc_ANC(10)[c(9)], colfunc_ANC(10)[c(6)], colfunc_ANC(10)[c(1)]),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = TRUE
)

#TEST OVERLAP GENES

#TEST INTERSECT
intersection <- length(intersect(unique(SELECTION_INSECT$Gene),unique(pvalues_perm_analysis$Gene)))
background <- 16120 #ALLGENES, SELECTION RUN ON ALL ALIGNMENTS (FLEXIBILITY OF SPECIES)
phyper(intersection-1, length(unique(SELECTION_INSECT$Gene)), 
       background-length(unique(SELECTION_INSECT$Gene)), length(unique(pvalues_perm_analysis$Gene)),
       lower.tail=F)

#POS
intersection <- length(intersect(unique(SELECTION_INSECT$GENEPOS),unique(pvalues_perm_analysis$GENEPOS)))


OVERLAP_INSECT$GENEPOS <- paste(OVERLAP_INSECT$Gene, OVERLAP_INSECT$Pos)

length(unique(OVERLAP_INSECT$Gene))
length(unique(OVERLAP_INSECT$GENEPOS))




SEL_PLOT <- JOINED %>% 
  mutate(Selection = case_when(!is.na(pvalue) & !is.na(target) ~ "FADE positive & CAAS",
                               is.na(pvalue) & !is.na(target) ~ "FADE positive",
                               TRUE ~ NA
                               ))  %>% group_by(Trait, Selection) %>%
  count(Trait, Selection) %>%
  drop_na()

#test <- JOINED %>%
       drop_na()
#View(test)
"""

  