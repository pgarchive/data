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

TABLE_AAs_contr <- TABLE_AAs[which(!TABLE_AAs$Trait %in% traits_with_no_contrast),]

#AAs_palette <- randomcoloR::distinctColorPalette(k = 21)
#saveRDS(AAs_palette, "AAs_colors")

AAs_palette <- readRDS("AAs_colors")
TABLE_AAs <- read.csv("Table_CAAs_SCENARIO1_2.txt", 
                      header=F, sep = "\t")

TABLE_AAs_CONTR <- TABLE_AAs[which(!TABLE_AAs$V2 %in% traits_with_no_contrast),]
nrow(TABLE_AAs_CONTR[which(TABLE_AAs_CONTR$V6<0.05),])

#head(Species_power)
colnames(TABLE_AAs_CONTR) <- c("ID", "Trait", "Pos", "AATop", "AABottom", "pvalue", "scenario")
TABLE_AAs_CONTR$transition <- paste(TABLE_AAs_CONTR$AATop, TABLE_AAs_CONTR$AABottom, sep="/")

Counts_trans <- as.data.frame(table(TABLE_AAs_CONTR$transition))
#ggplot(data=Counts_trans) +
#  geom_bar(aes(x=reorder(Var1, -Freq), y=Freq),
#                 stat="identity") + coord_flip() 

SUST_MATRIX <- read.csv("BLOSUM62.txt", 
                      header=T, sep = "")

matrix_AAs <- melt(as.matrix(SUST_MATRIX))

MERGED <- TABLE_AAs_CONTR %>% left_join(matrix_AAs, by=c("AATop" = "Var1", "AABottom" = "Var2")) %>%
  drop_na()

#nrow(MERGED[which(MERGED$pvalue <0.05),])
p2 <- ggplot(data=MERGED[which(MERGED$pvalue <0.05),]) +
  geom_histogram(aes(x=value, fill=value),
                 stat="count") +
  scale_fill_manual(values = AAs_palette)

saveRDS(p2, "FigureS33")
