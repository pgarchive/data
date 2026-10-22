##########################################################
##########################################################
##########################################################
##########################################################
##########SORTED HOMOGENEITY PLOT ANALYSIS ###############
##########################################################

library(ggplot2)
library(gridExtra)
library(tibble)
library(data.table)
library(ggrepel)
library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)

TABLE_AAs <- read.csv("Table_AAs_TOP.txt", 
                       header=F, sep = "\t")
#head(Species_power)
colnames(TABLE_AAs) <- c("ID", "Trait", "Pos", "AATop", "pvalue", "scenario")


traits_with_no_contrast <- c("Basophils.Abs",
                             "Basophils..",
                             "Blood.Urea.Nitrogen",
                             "Bun.Creat.Ratio",
                             "Cholesterol",
                             "Diencephalon",
                             "GGT",
                             "Globulin",
                             "Hippocampus",
                             "Hypothalamus",
                             "Maximum.lifespan..y.",
                             "M_Max_M_NeonateWeight_g",
                             "MOB",
                             "Piriform.Lobe",
                             "Schizocortex",
                             "Ser.Glut.Oxaloacetic.Transamin.",
                             "Striatum",
                             "Striatum..incl..NAcc.",
                             "Thalamus",
                             "Uric.Acid",
                             "XII")
TABLE_AAs_contr <- TABLE_AAs[which(!TABLE_AAs$Trait %in% traits_with_no_contrast),]

#AAs_palette <- randomcoloR::distinctColorPalette(k = 21)
#saveRDS(AAs_palette, "AAs_colors")

AAs_palette <- readRDS("AAs_colors")
ggplot(data=TABLE_AAs_contr) +
  geom_histogram(aes(x=AATop, fill=AATop),
                 stat="count") +
  scale_fill_manual(values = AAs_palette) +
  facet_wrap(~Trait ,scales="free")

#We have the species with available AA there and the species available for the trait

#READ INTERNAL P_VALUES
ggplot(data=TABLE_AAs_contr) +
  geom_histogram(aes(x=AATop, fill=AATop),
                 stat="count") +
  scale_fill_manual(values = AAs_palette) +
  theme_minimal() +
  coord_flip() +
  theme(legend.position = "none")

p2 <- ggplot(data=TABLE_AAs_contr[which(TABLE_AAs_contr$pvalue <0.05),]) +
  geom_histogram(aes(x=AATop, fill=AATop),
                 stat="count") +
  scale_fill_manual(values = AAs_palette) +
  theme_minimal() +
  coord_flip() +
  theme(legend.position = "none")

saveRDS(p2, "FigureS32")


