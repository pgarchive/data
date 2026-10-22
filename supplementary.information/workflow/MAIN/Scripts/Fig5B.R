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

pvalues_perm_sign <- pvalues_perm[which(pvalues_perm$pvalue <0.05),]
pvalues_perm_sign <- pvalues_perm[which(pvalues_perm$pvalue <0.05),]

#II) OBTAIN VENN DIAGRAM PLOT

library("ggvenn")
library(ggVennDiagram)

#View(FILT_CASES_OVERLAP)

#pvalues_perm_analysis <- pvalues_perm[which(as.numeric(as.character(pvalues_perm$pvalue)) <0.05),]
#ALL, NOT ONLY PVAL below 0.5
pvalues_perm_analysis <- pvalues_perm_sign
pvalues_perm_analysis$GENEPOS <- paste(pvalues_perm_analysis$Gene, pvalues_perm_analysis$Pos)
x <- list(
  "S-CAAS" = pvalues_perm_analysis$Gene[as.character(pvalues_perm_analysis$Trait) == "Maximum.lifespan..y."],
  "A-CAAS" = pvalues_perm_analysis$Gene[as.character(pvalues_perm_analysis$Trait) == "Maximum_lifespan__y__altern1_config.txt"], 
  "B-CAAS" = pvalues_perm_analysis$Gene[as.character(pvalues_perm_analysis$Trait) == "Maximum_lifespan__y__altern2_config.txt"]
)

names(x) <- c("S-CAAS", "A-CAAS", "B-CAAS")

#VERSION FOR STUDYING OVERLAP WITH POSITIONS WITH SIGNIFICANT PGLS
#x <- list(
# "S-PGLS<0.05" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "Maximum.lifespan..y." &
#                                           as.numeric(as.character(pvalues_perm_analysis$pval)) < 0.05],
#"A-PGLS<0.05" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "Maximum_lifespan__y__altern1_config.txt" &
#                                          as.numeric(as.character(pvalues_perm_analysis$pval)) < 0.05],
#"B-PGLS<0.05" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "Maximum_lifespan__y__altern2_config.txt" &
#                                          as.numeric(as.character(pvalues_perm_analysis$pval)) < 0.05]
#)

#y <- list(
# "S-PGLS>0.05" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "Maximum.lifespan..y." &
#                                                as.numeric(as.character(pvalues_perm_analysis$pval)) > 0.05],
#"A-PGLS>0.05" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "Maximum_lifespan__y__altern1_config.txt" &
#                                               as.numeric(as.character(pvalues_perm_analysis$pval)) > 0.05],
#"B-PGLS>0.05" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "Maximum_lifespan__y__altern2_config.txt" &
#                                               as.numeric(as.character(pvalues_perm_analysis$pval)) > 0.05]
#)
#forget - nothing to extract here #para el carro si cuesta mucho algo
#View(intersect(intersect(x$A,x$B),x$G))

#names(x) <- c("S-PGLS<0.05", "A-PGLS<0.05", "B-PGLS<0.05")
#names(y) <- c("S-PGLS>0.05", "A-PGLS>0.05", "B-PGLS>0.05")
traits_palette <- readRDS("traits_palette")
colfunc_ANC <- colorRampPalette(c(traits_palette[3], "white"))

#draw plot
library(ggvenn)
gg5B2 <- ggvenn(
  x, 
  fill_color = c(colfunc_ANC(10)[c(9)], colfunc_ANC(10)[c(6)], colfunc_ANC(10)[c(1)]),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = TRUE
)

#gg5B2_PGLSsign <- ggVennDiagram(
# x, 
#  category.names = c("S-PGLS<0.05", "A-PGLS<0.05", "B-PGLS<0.05"),
# label = "both",
#label_percent_digit = 2,
#label_alpha = 0,
#edge_size = 0
#) +
# scale_fill_gradient(low="white", high=traits_palette[3]) +
#  scale_x_continuous(expand = expansion(mult = .2)) +
# theme(legend.position = "none")

#gg5B2_PGLSnonsign <- ggVennDiagram(
# y, 
#category.names = c("S-PGLS>0.05", "A-PGLS>0.05", "B-PGLS>0.05"),
#label = "both",
#label_percent_digit = 2,
#label_alpha = 0,
#edge_size = 0
#) +
# scale_fill_gradient(low="white", high=traits_palette[3]) +
#scale_x_continuous(expand = expansion(mult = .2)) +
#theme(legend.position = "none")

saveRDS(gg5B2, "gg5B2_FEB24")

#Insect
x <- list(
  G = pvalues_perm_analysis$Gene[as.character(pvalues_perm_analysis$Trait) == "TrophicGuild_Insectivore"],
  A = pvalues_perm_analysis$Gene[as.character(pvalues_perm_analysis$Trait) == "TrophicGuild_Insectivore_altern1_config.txt"], 
  B = pvalues_perm_analysis$Gene[as.character(pvalues_perm_analysis$Trait) == "TrophicGuild_Insectivore_altern2_config.txt"]
)

names(x) <- c("S-CAAS","A-CAAS","B-CAAS")
traits_palette <- readRDS("traits_palette")
colfunc_ANC <- colorRampPalette(c(traits_palette[1], "white"))

#draw plot
gg5B1 <- ggvenn(
  x, 
  fill_color = c(colfunc_ANC(10)[c(9)], colfunc_ANC(10)[c(6)], colfunc_ANC(10)[c(1)]),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = TRUE
)

saveRDS(gg5B1, "gg5B1_FEB24")

#VERSION FOR STUDYING OVERLAP WITH POSITIONS WITH SIGNIFICANT PGLS
#x <- list(
# "S-PhyloLM<0.05" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "TrophicGuild_Insectivore" &
#                                                as.numeric(as.character(pvalues_perm_analysis$pval)) < 0.05],
#"A-PhyloLM<0.05" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "TrophicGuild_Insectivore_altern1_config.txt" &
#                                               as.numeric(as.character(pvalues_perm_analysis$pval)) < 0.05],
#"B-PhyloLM<0.05" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "TrophicGuild_Insectivore_altern2_config.txt" &
#                                               as.numeric(as.character(pvalues_perm_analysis$pval)) < 0.05]
#)

#y <- list(
# "S-PhyloLM>0.05" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "TrophicGuild_Insectivore" &
#                                                as.numeric(as.character(pvalues_perm_analysis$pval)) > 0.05],
#"A-PhyloLM>0.05" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "TrophicGuild_Insectivore_altern1_config.txt" &
#                                               as.numeric(as.character(pvalues_perm_analysis$pval)) > 0.05],
#"B-PhyloLM>0.05" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "TrophicGuild_Insectivore_altern2_config.txt" &
#                                               as.numeric(as.character(pvalues_perm_analysis$pval)) > 0.05]
#)
#forget - nothing to extract here #para el carro si cuesta mucho algo
#View(intersect(intersect(x$A,x$B),x$G))

#names(x) <- c("S-PhyloLM<0.05", "A-PhyloLM<0.05", "B-PhyloLM<0.05")
#names(y) <- c("S-PhyloLM>0.05", "A-PhyloLM>0.05", "B-PhyloLM>0.05")
#traits_palette <- readRDS("traits_palette")
#colfunc_ANC <- colorRampPalette(c(traits_palette[3], "white"))

#draw plot
#gg5B1_PGLSsign <- ggVennDiagram(
# x, 
#category.names = c("S-PhyloLM<0.05", "A-PhyloLM<0.05", "B-PhyloLM<0.05"),
#label = "both",
#label_percent_digit = 2,
#label_alpha = 0,
#edge_size = 0
#) +
# scale_fill_gradient(low="white", high=traits_palette[1]) +
#  scale_x_continuous(expand = expansion(mult = .2)) +
# theme(legend.position = "none")

#gg5B1_PGLSnonsign <- ggVennDiagram(
# y, 
#  category.names = c("S-PhyloLM>0.05", "A-PhyloLM>0.05", "B-PhyloLM>0.05"),
# label = "both",
#label_percent_digit = 2,
#label_alpha = 0,
#edge_size = 0
#) +
# scale_fill_gradient(low="white", high=traits_palette[1]) +
#scale_x_continuous(expand = expansion(mult = .2)) +
#theme(legend.position = "none")




#WBC
x <- list(
  G = pvalues_perm_analysis$Gene[as.character(pvalues_perm_analysis$Trait) == "White.Blood.Cells"],
  A = pvalues_perm_analysis$Gene[as.character(pvalues_perm_analysis$Trait) == "White_Blood_Cells_altern1_config.txt"], 
  B = pvalues_perm_analysis$Gene[as.character(pvalues_perm_analysis$Trait) == "White_Blood_Cells_altern2_config.txt"]
)

names(x) <- c("S-CAAS","A-CAAS","B-CAAS")
traits_palette <- readRDS("traits_palette")
colfunc_ANC <- colorRampPalette(c("grey", "white"))

#draw plot
gg5B3 <- ggvenn(
  x, 
  fill_color = c(colfunc_ANC(10)[c(9)], colfunc_ANC(10)[c(6)], colfunc_ANC(10)[c(1)]),
  stroke_size = 0.5, set_name_size = 4,
  show_percentage = TRUE
)

saveRDS(gg5B3, "gg5B3_FEB24")

#get p-values

#Internal_lambda_adj <- Internal_lambda %>%
# group_by(trait) %>%
#mutate(pval.adj = p.adjust(pval, method='fdr', 
#                               n = length(pval[!is.na(pval)])))

#View(Internal_lambda)









#################################
################################
#POSITIONS NOT SEGREGATING ######
#################################

#filt segregating cases
#FILT_CASES_OVERLAP_NOTSEGR <- joined_PGLS[which(joined_PGLS$N_top <= 5 |
#                                         joined_PGLS$N_bottom <= 5),]
#FILT_CASES_OVERLAP_SEGR <- joined_PGLS[which(as.numeric(as.character(joined_PGLS$N_top)) > 5 &
#                                     as.numeric(as.character(joined_PGLS$N_bottom)) > 5),] %>%
#  drop_na()

#View(FILT_CASES_OVERLAP_SEGR)
#nrow(joined_PGLS)
#nrow(FILT_CASES_OVERLAP_NOTSEGR)
#nrow(FILT_CASES_OVERLAP_SEGR)

#pvalues_perm_analysis <- FILT_CASES_OVERLAP_NOTSEGR
#pvalues_perm_analysis$GENEPOS <- paste(pvalues_perm_analysis$Gene, pvalues_perm_analysis$Pos)

#z <- list(
# "S-Cons" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "Maximum.lifespan..y."],
#  "A-Cons" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "Maximum_lifespan__y__altern1_config.txt"],
# "B-Cons" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "Maximum_lifespan__y__altern2_config.txt"]
#)


#names(z) <- c("S-Cons", "A-Cons", "B-Cons")
#traits_palette <- readRDS("traits_palette")
#colfunc_ANC <- colorRampPalette(c(traits_palette[3], "white"))

#draw plot
#gg5B2_Cons <- ggVennDiagram(
#  z, 
#category.names = c("S-Cons", "A-Cons", "B-Cons"),
#label = "both",
#label_percent_digit = 2,
#label_alpha = 0,
#edge_size = 0
#) +
# scale_fill_gradient(low="white", high=traits_palette[3]) +
#  scale_x_continuous(expand = expansion(mult = .2)) +
# theme(legend.position = "none")

#Insectivority
#z <- list(
# "S-Cons" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "TrophicGuild_Insectivore"],
#"A-Cons" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "TrophicGuild_Insectivore_altern1_config.txt"],
#"B-Cons" = pvalues_perm_analysis$GENEPOS[as.character(pvalues_perm_analysis$Trait) == "TrophicGuild_Insectivore_altern2_config.txt"]
#)


#names(z) <- c("S-Cons", "A-Cons", "B-Cons")
#traits_palette <- readRDS("traits_palette")
#colfunc_ANC <- colorRampPalette(c(traits_palette[3], "white"))

#draw plot
#gg5B2_Cons <- ggVennDiagram(
# z, 
#  category.names = c("S-Cons", "A-Cons", "B-Cons"),
# label = "both",
#  label_percent_digit = 2,
# label_alpha = 0,
#  edge_size = 0
#) +
# scale_fill_gradient(low="white", high=traits_palette[1]) +
#scale_x_continuous(expand = expansion(mult = .2)) +
#theme(legend.position = "none")


#PRUEBA CON SELECTION
#Read selection files
#SELECTION <- read.csv("ALL_SELECTION_STUDY_CASES_FILT.txt", sep ="\t", header=FALSE)
#colnames(SELECTION) <- c("Trait", "Gene", "Pos",	"Partition", "target", "rate", 
#                        "bias", "Bayes Factor", "Amino acid composition of site",	
#                       "Substitution history on selected branches")

#SELECTION$Pos0 <- SELECTION$Pos-1

#levels(joined_PGLS$Trait) <- c(levels(joined_PGLS$Trait),
#                                        c("Maximum_lifespan__y__G",
#                                         "Maximum_lifespan__y__A", "Maximum_lifespan__y__B",
#                                        "TrophicGuild_Insectivore_G",
#"TrophicGuild_Insectivore_A",
#                                     "TrophicGuild_Insectivore_B", 
#                                    "White_Blood_Cells_G", "White_Blood_Cells_A", "White_Blood_Cells_B"))
#joined_PGLS$Trait[joined_PGLS$Trait == "Maximum_lifespan__y__altern1_config.txt"] <- "Maximum_lifespan__y__A"
#joined_PGLS$Trait[joined_PGLS$Trait == "Maximum_lifespan__y__altern2_config.txt"] <- "Maximum_lifespan__y__B"
#joined_PGLS$Trait[joined_PGLS$Trait == "Maximum.lifespan..y."] <- "Maximum_lifespan__y__G"
#joined_PGLS$Trait[joined_PGLS$Trait == "TrophicGuild_Insectivore_altern1_config.txt"] <- "TrophicGuild_Insectivore_A"
#joined_PGLS$Trait[joined_PGLS$Trait == "TrophicGuild_Insectivore_altern2_config.txt"] <- "TrophicGuild_Insectivore_B"
#joined_PGLS$Trait[joined_PGLS$Trait == "TrophicGuild_Insectivore"] <- "TrophicGuild_Insectivore_G"
#joined_PGLS$Trait[joined_PGLS$Trait == "White_Blood_Cells_altern2_config.txt"] <- "White_Blood_Cells_B"
#joined_PGLS$Trait[joined_PGLS$Trait == "White.Blood.Cells"] <- "White_Blood_Cells_G"



#JOINED_SELECTION <-   left_join(SELECTION, joined_PGLS,
#                               by=c("Gene" = "Gene",
#                                   "Trait" = "Trait",
#                                   "Pos0" = "Pos")) 
#View(JOINED_SELECTION)

#ggplot(data = JOINED_SELECTION) +
# geom_histogram(aes(x=N_top, 
#                   fill=Trait)) +
#scale_fill_manual(values=c(colorRampPalette(c(traits_palette[3], "white"))(10)[1],
#                 colorRampPalette(c(traits_palette[3], "white"))(10)[3],
#                  colorRampPalette(c(traits_palette[1], "white"))(10)[9],
#                  colorRampPalette(c(traits_palette[1], "white"))(10)[1],
#                  colorRampPalette(c(traits_palette[1], "white"))(10)[3],
#                  colorRampPalette(c(traits_palette[5], "white"))(10)[9],
#                  colorRampPalette(c(traits_palette[5], "white"))(10)[1]))
#nrow(joined_PGLS[which(is.na(joined_PGLS$pval)),])
#nrow(joined_PGLS)

#ggplot(data = joined_PGLS[which(is.na(joined_PGLS$pval)),]) +
# geom_histogram(aes(x=N_top, ,
#                   fill=Trait)) +
#facet_wrap(~Trait, scales = "free") +
#scale_fill_manual(values=c(colorRampPalette(c(traits_palette[5], "white"))(10)[3],
#                           colorRampPalette(c(traits_palette[3], "white"))(10)[1],
#                           colorRampPalette(c(traits_palette[3], "white"))(10)[3],
#                           colorRampPalette(c(traits_palette[3], "white"))(10)[9],
#                          colorRampPalette(c(traits_palette[1], "white"))(10)[9],
#                         colorRampPalette(c(traits_palette[1], "white"))(10)[1],
#                        colorRampPalette(c(traits_palette[1], "white"))(10)[3],
#                       colorRampPalette(c(traits_palette[5], "white"))(10)[9],
#                      colorRampPalette(c(traits_palette[5], "white"))(10)[1]))


#ggplot(data = joined_PGLS[which(is.na(joined_PGLS$pval)),]) +
# geom_count(aes(x=N_top, y = N_bottom,
#                   fill=Trait)) +
#facet_wrap(~Trait, scales = "free") +
#scale_fill_manual(values=c(colorRampPalette(c(traits_palette[5], "white"))(10)[3],
#                           colorRampPalette(c(traits_palette[3], "white"))(10)[1],
#                           colorRampPalette(c(traits_palette[3], "white"))(10)[3],
#                           colorRampPalette(c(traits_palette[3], "white"))(10)[9],
#                           colorRampPalette(c(traits_palette[1], "white"))(10)[9],
#                           colorRampPalette(c(traits_palette[1], "white"))(10)[1],
#                           colorRampPalette(c(traits_palette[1], "white"))(10)[3],
#                           colorRampPalette(c(traits_palette[5], "white"))(10)[9],
#                           colorRampPalette(c(traits_palette[5], "white"))(10)[1]))


