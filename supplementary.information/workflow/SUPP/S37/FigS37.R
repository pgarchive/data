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
library(ggplot2)
library(gridExtra)
library(tibble)
library(data.table)
library(ggrepel)
library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)


##################
##CAAS POSITIONS##
##################

path <- "VARIATION_files/"
list.directories <- list.files(path, pattern = "APR23.txt$")
#fields <- strsplit(list.directories, fixed = TRUE)
genes <- c()

ldf_CAAS <- list()
count <- 0
for (j in 1:length(list.directories)){
  count <- count + 1
  filepath <- paste(path, list.directories[j], sep="/")
  spc_df  <- read.csv(filepath, sep ="\t", header = FALSE)
  spc_df$V4 <- as.numeric(as.character(spc_df$V4))
  ldf_CAAS[[count]] <- as.data.frame(spc_df)}

CAAS_DF <- do.call(rbind, ldf_CAAS)
#nrow(CAAS_DF)
#CAAS_DF[which(CAAS_DF$V2 == "X"),]

CAAS_BIALLELIC <- CAAS_DF[!grepl(",", CAAS_DF$V5),]
#nrow(CAAS_BIALLELIC)

#tail(CAAS_BIALLELIC)
colnames(CAAS_BIALLELIC) <- c("GENE", "Chr", "SNP", "Pos", "AF")


#tail(CAAS_BIALLELIC$Pos)
#readTRANSVAR PHYLOPs

transvar <- read.csv("../S36/PHYLOP_files/caas_phyloP/ALL_CAAS_PHYLOPs_APR23.txt", sep="\t", header=F)
phylops_caas <- transvar[,1:5]
colnames(phylops_caas) <- c("Gene", "Chr", "Pos", "strand", "PhyloP")

CAAS_with_phyloP_df <- left_join(CAAS_BIALLELIC, phylops_caas, by = c("Chr", "Pos")) %>%
  drop_na()


#MERGE WITH ONLY CONTRAST PHENOTYPES
transvar_info <- read.csv("../S36/ALL_TRANSVAR_RESULTS_APR23_FORMATTED_FILT.txt", sep ="\t", header=FALSE)
colnames(transvar_info) <- c("Chrom", "Coord", "SNP", "Gene", "mutation", "annot", "strand")
#Input transvar
transvar_extra_info <- read.csv("../S35/ALL_TRANSVAR_INFO_APR23.txt", sep ="\t", header=FALSE)
colnames(transvar_extra_info) <- c("Gene", "Trait", "Pos", "mutation")
ALL_TRANSVAR_INFO_JOINED <- left_join(transvar_extra_info, transvar_info, by =c("Gene" = "Gene", 
                                                                                "mutation" = "mutation"))
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
ALL_TRANSVAR_INFO_JOINED_contr <- ALL_TRANSVAR_INFO_JOINED[which(!ALL_TRANSVAR_INFO_JOINED$Trait %in% traits_with_no_contrast),]

CAAS_with_phyloP_df_contr <- left_join(ALL_TRANSVAR_INFO_JOINED_contr, CAAS_with_phyloP_df, 
                                       by =c("Gene" = "Gene", 
                                             "Chrom" = "Chr",
                                             "Coord"="Pos")) %>%
  drop_na() %>%
  distinct(GENE, Chrom, SNP.y, Coord, AF, strand.y, PhyloP)

FIXED <- CAAS_with_phyloP_df_contr %>%
  mutate(fixed = ifelse(AF == "NULL", "fixed", "variant")) %>%
  drop_na() %>%
  filter(fixed == "fixed") %>%
  distinct(GENE, Chrom, SNP.y, Coord, AF, strand.y, PhyloP)

INVARIANT <- CAAS_with_phyloP_df_contr %>%
  mutate(fixed = ifelse(AF == "NULL", "fixed", "variant")) %>%
  drop_na()

INVARIANT$intervals <- cut(INVARIANT$PhyloP, breaks = seq(-22, 14, 2))

VARIANT <- INVARIANT[which(INVARIANT$fixed == "variant"),]

palette <- randomcoloR::distinctColorPalette(2)

#p1 <- ggplot(data=VARIANT[which(as.numeric(as.character(VARIANT$AF)) > 0.01),], aes(x=intervals, y=as.numeric(as.character(AF)))) +
#  geom_violin() +
#  facet_wrap(~Name, scales="free")

INVARIANT$intervals <- droplevels(INVARIANT$intervals)
p2 <- ggplot(data=INVARIANT) +
  geom_histogram(aes(x=intervals, fill=fixed),
                 stat='count',position="dodge2") +
  scale_fill_manual(values=palette[c(2,1)]) +
  labs(fill="Human position gnomAD") + xlab("Phylo-P intervals") +
  theme_minimal() +
  theme(legend.position = "top",
        axis.text.x = element_text(angle=90))
saveRDS(p2, "FigureS37")

