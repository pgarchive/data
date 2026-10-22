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

#READ TREE
primates_tree <- read.newick("../Fig1/FINAL_TREE_BL_WITHOUTNODELAB.newick")

#A) READ LAST TRAITS FILE and metadata

all_traits <- read.csv("../Fig1/PRIMATES_PHENOMIC_DATASET.tsv", header=T, sep = "\t")

#Now let's bind the Group info


#B) SELECT DESIRED SPECIES (sequenced_ones)

primates_palette <- readRDS("../Fig1/primate_families_palette.rds")
seq_species <- read.csv("../Fig1/names_FINAL.txt",
                        sep = "\t", header = TRUE)
colnames(seq_species)<- c("Tree_Species", "Genome_species", "Trait_species")
covered_species_all <- as.character(seq_species$Trait_species[!is.na(seq_species$Tree_Species)])
covered_species <- as.character(covered_species_all[!is.na(covered_species_all)])
sequenced_traits <- all_traits[which(all_traits$SpeciesBROAD %in% covered_species),]
levels(sequenced_traits$SpeciesBROAD) <- c(levels(sequenced_traits$SpeciesBROAD), "Tarsius_bancanus", "Leontocebus_fuscicolis",
                                           "Plecturocebus_hoffmansi") 
sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Cephalopachus_bancanus"] <- "Tarsius_bancanus"
sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Leontocebus_fuscicollis"] <- "Leontocebus_fuscicolis"
sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Plecturocebus_hoffmannsi"] <- "Plecturocebus_hoffmansi"


groups_list <- strsplit(as.character(sequenced_traits$GroupName), "_")
groups <- c()
for(i in 1:length(groups_list)){
  group <- paste(groups_list[[i]][1], "_", groups_list[[i]][2], sep="")
  groups[i] <- group}
groups
sequenced_traits <- cbind(sequenced_traits, groups)

#ADD QUALITATIVE HUMAN MEASUREMENTS

library(janitor)
human_measurements <- read.csv("../Fig1/human_direct_measurements.tsv", header=TRUE, sep = "\t")
human_transposed <- t(human_measurements) %>% row_to_names(row_number = 1) %>% head(1)
#convert to numeric columns
human_transposed <- as.data.frame(human_transposed) %>% 
  dplyr::select(Social_organization,
                Mating_system, 
                DielActivity,
                Locomotion,
                TrophicGuild,
                Paternal_care,
                Female_Ranging_patterns,
                Infanticide,
                Seasonal_breeders,
                Colour_vision)

sequenced_traits<- dplyr::bind_rows(sequenced_traits, human_transposed)
levels(sequenced_traits$SpeciesBROAD) <- c(levels(sequenced_traits$SpeciesBROAD),"Homo_sapiens")
sequenced_traits[nrow(sequenced_traits), "SpeciesBROAD"] <- "Homo_sapiens"
sequenced_traits[nrow(sequenced_traits), "groups"] <- "APE_Hominidae"

#remove ursinus
sequenced_traits <- sequenced_traits[which(sequenced_traits$SpeciesBROAD != "Papio_ursinus"),]

Traits_top <- read.csv("../Fig1/Qualitative_extremes", header=T, sep="\t")
INSECTIVORE_DATASET <- Traits_top[which(Traits_top$Trait == "TrophicGuild_Insectivore"),]
pruned_primates.tree<-drop.tip(primates_tree ,primates_tree$tip.label[-match(INSECTIVORE_DATASET$label[!is.na(INSECTIVORE_DATASET$label)], 
                                                                             primates_tree$tip.label)])

#pruned_primates.tree <- force.ultrametric(pruned_primates.tree)
insect <- as.numeric(as.character(INSECTIVORE_DATASET$value))
names(insect) <- INSECTIVORE_DATASET$label


#ALTERNATIVE WITH GGTREE

library(ggtree)

#VII) CREATE PLOT WITH Bayesian values for ancestral states in maximum lifespan phylogeny
# estimate continuous ASRs with Bayesian Brownian motion model
asr.bi <- data.frame(rerootingMethod(pruned_primates.tree, 
                                     insect)$marginal.anc)

asr.bi$node <- row.names(asr.bi)

library(magick)

p12 <- ggtree(force.ultrametric(pruned_primates.tree),
              layout='circular')

traits_palette <- readRDS("../Fig1/traits_palette")
colfunc_ANC <- colorRampPalette(c(traits_palette[1], "white"))
pies <- nodepie(asr.bi, cols=1:2, color=c(colfunc_ANC(10)[c(5,1)]), 
                alpha=0.8)
df_pies <- tibble::tibble(node=as.numeric(asr.bi$node), pies=pies)

#p12 <- inset(p12, pies, width=0.004, height =0.006)
#plot(p12)

#Annotate alternative Tops and Bottoms

Global_sel <- read.csv("TrophicGuild_Insectivore_config.txt",
                       sep="\t",
                       header=F)

AlterA <- read.csv("TrophicGuild_Insectivore_altern1_config.txt",
                   sep="\t",
                   header=F)

AlterB <- read.csv("TrophicGuild_Insectivore_altern2_config.txt",
                   sep="\t",
                   header=F)

selGlobal <- c()
selA <- c()
selB <- c()
for (i in 1:nrow(p12$data)){
  species <- p12$data$label[i]
  if (species %in% as.character(Global_sel$V1)){
    curr_global <- as.numeric(as.character(Global_sel$V2[as.character(Global_sel$V1) == species]))
    selGlobal <- c(selGlobal, curr_global)}
  else{
    selGlobal <- c(selGlobal, "other")
  }
  if (species %in% as.character(AlterA$V1)){
    currA <- as.numeric(as.character(AlterA$V2[as.character(AlterA$V1) == species]))
    selA <- c(selA, currA)}
  else{
    selA <- c(selA, "other")
  }
  if (species %in% as.character(AlterB$V1)){ 
    currB <- as.numeric(as.character(AlterB$V2[as.character(AlterB$V1) == species]))
    selB <- c(selB, currB)}
  else{
    selB <- c(selB, "other")
  }
}

p12$data$selGlobalValues <- selGlobal
p12$data$selAValues <- selA
p12$data$selBValues <- selB

library(ggtree)
library(ggimage)
library(ggtreeExtra)
library(ggplot2)


species_df <- p12$data[which(!is.na(p12$data$label)),]

SPECIES_DF <- species_df %>% select(label, selGlobalValues,
                                    selAValues, selBValues) %>%
  mutate(value=abs(rnorm(207)))
colnames(SPECIES_DF) <- c("label", "Global", "A", "B", "value")
SPECIES_DF <- as.data.frame(SPECIES_DF)
rownames(SPECIES_DF) <- SPECIES_DF$label

#tree
p11 <- ggtree(force.ultrametric(pruned_primates.tree),
              layout='circular')

#add families to tree

families <- c()
for(i in 1:length(p11$data$label)){
  if (as.character(p11$data$label[i]) %in% as.character(sequenced_traits$SpeciesBROAD)){
    family <- as.character(sequenced_traits$groups[sequenced_traits$SpeciesBROAD == as.character(p11$data$label[i])])
    families <- c(families, family)}
  else{
    family <- NA
    families <- c(families, family)
  }}

p11$data$Family <- families
primate_families_palette<- readRDS("primate_families_palette.rds")

#ADD FAMILIES SATURATED TO PHYLOGENY

desat <- function(cols, sat=0.5) {
  X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}

primate_families_palette_des <- desat(primate_families_palette, 0.25)


p12 <-p11 + geom_cladelabel(node=269, label="", 
                            offset=.001, align=TRUE,
                            barsize = 1.5, color=primate_families_palette_des[3]) +
  geom_cladelabel(node=354, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[2]) +
  geom_cladelabel(node=348, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[1]) +
  geom_cladelabel(node=217, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[6]) +
  geom_cladelabel(node=215, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[4]) +
  geom_cladelabel(node=237, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[7]) +
  geom_cladelabel(node=244, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[5]) +
  geom_cladelabel(node=256, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[8]) +
  geom_cladelabel(node=362, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[16]) +
  geom_cladelabel(node=385, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[13]) +
  geom_cladelabel(node=369, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[12]) +
  geom_cladelabel(node=382, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[14]) +
  geom_cladelabel(node=379, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[9])+
  geom_cladelabel(node=158, label="", extend=0.5,
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[10])+
  geom_cladelabel(node=407, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[15])+
  geom_cladelabel(node=403, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[11])


insect_df <- as.data.frame(cbind(names(insect), 
                                 insect))
colnames(insect_df) <- c("label", "insect")

p13 <- p12 %<+% insect_df

#View(p13$data)

OUTLIERS <- SPECIES_DF[which(SPECIES_DF$Global == 1),]
NON_OUTLIERS <- SPECIES_DF[which(!SPECIES_DF$Global == 1),]
#p13 <- p13 + 
# aes(color=as.numeric(as.character(gummi))) +
#scale_color_gradientn(colors=c("black",
#                              colfunc_ANC(10)[c(5)]))

p12 +
  geom_text(aes(label=node))

p13 <- p12 + 
  geom_fruit(data = NON_OUTLIERS,
             geom = geom_bar,
             stat="identity",
             mapping = aes(x=1,
                           y =label,
                           color=Global,
                           shape=Global),
             offset=0.13,
             pwidth=0.002,
             size=0.5,
             axis.params=list(
               axis="x", title="S",
               text.angle = 0,
               text.size = 0,
               line.size = 0,
               vjust = 0, 
               family="Helvetica", fontface="bold"
             )) +
  #scale_color_manual(values = "black") +
  scale_color_manual(values = c(colfunc_ANC(10)[c(7)],
                                "grey40")) + 
  geom_fruit(data = OUTLIERS,
             geom = geom_point,
             mapping = aes(y =label,
                           fill=Global,
                           shape=Global),
             offset=0,
             pwidth=0.002,
             size=2
  ) +
  scale_shape_manual(values = c(22, 24,1)) +
  #scale_color_manual(values = "black") +
  scale_fill_manual(values = c(colfunc_ANC(10)[c(5)])) 


p13$scales$scales <- lapply(p16$scales$scales, function(s) {
  if (!is.null(s$guide) && !is.atomic(s$guide) && !is.list(s$guide)) {
    s$guide <- "legend"
  }
  s
})

OUTLIERS_A <- SPECIES_DF[which(SPECIES_DF$A == 1),]
NON_OUTLIERS_A <- SPECIES_DF[which(!SPECIES_DF$A == 1),]
p14 <- p13 +  new_scale("shape") +
  geom_fruit(data = NON_OUTLIERS_A,
             geom = geom_bar,
             stat="identity",
             mapping = aes(x=1,
                           y =label,
                           color=A,
                           shape=A),
             offset=0.07,
             pwidth=0.002,
             size=0.5,
             axis.params=list(
               axis="x", title="A",
               text.angle = 0,
               text.size = 0,
               line.size = 0,
               vjust = 0, 
               family="Helvetica", fontface="bold"
             )) +
  #scale_color_manual(values = "black") +
  scale_color_manual(values = c(colfunc_ANC(10)[c(7)],
                                "grey40")) +
  geom_fruit(data = OUTLIERS_A,
             geom = geom_point,
             mapping = aes(y=label,
                           fill=A,
                           shape =A),
             pwidth=0.002,
             offset=0,
             size=2) +  scale_shape_manual(values = c(22, 24,1)) +
  scale_fill_manual(values = c(colfunc_ANC(10)[c(5)]))

library(ggpp)
p15 <- p14 %<+% df_pies[which(df_pies$node %in% c(208,
                                                  209,
                                                  362,
                                                  365, 366, 402,
                                                  210, 268,
                                                  408, 412, 407, 348,
                                                  347, 409)),]

p16 <- p15 +
  geom_plot(data=td_filter(!isTip), 
            mapping=aes(x=x,y=y, label=pies), 
            vp.width=0.03, vp.height=0.03, 
            hjust=0.5, vjust=0.5) 
OUTLIERS_B <- SPECIES_DF[which(SPECIES_DF$B == 1),]
NON_OUTLIERS_B <- SPECIES_DF[which(!SPECIES_DF$B == 1),]

p17 <- p16 + 
  geom_fruit(data = NON_OUTLIERS_B,
             geom = geom_bar,
             stat="identity",
             mapping = aes(x=1,
                           y =label,
                           color=B,
                           shape=B),
             offset=0.07,
             pwidth=0.002,
             size=0.5,
             axis.params=list(
               axis="x", title="B",
               text.angle = 0,
               text.size = 0,
               line.size = 0,
               vjust = 0, 
               family="Helvetica", fontface="bold"
             )) +
  #scale_color_manual(values = "black") +
  scale_color_manual(values = c(colfunc_ANC(10)[c(7)],
                                "grey40")) +
  geom_fruit(data = OUTLIERS_B,
             geom = geom_point,
             mapping = aes(y=label,
                           fill=B,
                           shape =B),
             pwidth=0.002,
             offset=0,
             size=2) +  scale_shape_manual(values = c(22, 24,1)) +
  scale_fill_manual(values = c(colfunc_ANC(10)[c(5)])) +
  xlim(0,0.07) +
  theme(legend.position = "none") +
  ggtitle(expression(bolditalic("Insectivorous diet"))) +
  theme(plot.title = element_text(size = 15, hjust = 0.5,
                                  color = colfunc_ANC(10)[c(3)]))


saveRDS(p17, "Fig5A1_JAN24")

