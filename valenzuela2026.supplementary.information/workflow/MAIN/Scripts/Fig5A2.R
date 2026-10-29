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
library(ggnewscale)

#I) READ LAST TRAITS FILE and metadata

all_traits <- read.csv("../Fig1/PRIMATES_PHENOMIC_DATASET.tsv", header=T, sep = "\t")

#Now let's bind the Group info


#II) SELECT DESIRED SPECIES (sequenced_ones)

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

#III) READ PRIMATE FAMILIES (sequenced_ones)

groups_list <- strsplit(as.character(sequenced_traits$GroupName), "_")
groups <- c()
for(i in 1:length(groups_list)){
  group <- paste(groups_list[[i]][1], "_", groups_list[[i]][2], sep="")
  groups[i] <- group}
groups
sequenced_traits <- cbind(sequenced_traits, groups)

#IV) ADD HUMAN QUANTITATIVE MEASUREMENTS

library(janitor)
human_measurements <- read.csv("../Fig1/human_direct_measurements.tsv", header=TRUE, sep = "\t")
human_transposed <- t(human_measurements) %>% row_to_names(row_number = 1) %>% head(1)
#convert to numeric columns
human_transposed <- as.data.frame(human_transposed) %>% 
  mutate_at(vars(BodyMass_kg,
                 Body.mass..g., 
                 Body_mass,
                 Brain_mass,
                 Neonatal.mass..g.,
                 LitterSz,
                 Interbirth.interval..days.,
                 Gestation..days.,
                 Lactation..days.,
                 Age.of.female.at.first.reproduction..y.,
                 Maximum.lifespan..y.,
                 Maximum.possible.reproductive.success.under.optimal.conditions..rmax.,
                 Albumin,
                 Calcium,
                 Cholesterol,
                 Phosphorus,
                 Potassium,
                 Ser.Glut.Pyruvic.Transaminase,
                 Ser.Glut.Oxaloacetic.Transamin.,
                 Total.Protein,
                 BodyMassMale_kg,
                 BodyMassFemale_kg,
                 MassDimorphism,
                 Offspring.per.litter,
                 Weaning_prop,
                 Alkaline.Phosphatase,
                 Basophils..,
                 Chloride,
                 Hematocrit,
                 Hemoglobin,
                 Lactic.Dehydrogenase,
                 Lymphocytes..,
                 Magnesium,
                 Monocytes..,
                 Neutrophil.Seg..,
                 Red.Blood.Cells,
                 Sodium,
                 Triglycerides,
                 Uric.Acid,
                 White.Blood.Cells,
                 High.Density.Lipoprotein,
                 Iron,
                 BasalMetabolicRate_ml_O2.h,
                 RestingMetabolicRate_kcal.day), ~as.numeric(as.character(.))) 
#convert to integer
human_transposed <- as.data.frame(human_transposed) %>% 
  mutate_at(vars(Relative_tail_length), ~as.integer(as.character(.))) 

sequenced_traits <- dplyr::bind_rows(sequenced_traits, human_transposed)
levels(sequenced_traits$SpeciesBROAD) <- c(levels(sequenced_traits$SpeciesBROAD),"Homo_sapiens")
sequenced_traits[nrow(sequenced_traits), "SpeciesBROAD"] <- "Homo_sapiens"
sequenced_traits[nrow(sequenced_traits), "groups"] <- "APE_Hominidae"
sequenced_traits <- sequenced_traits[which(sequenced_traits$SpeciesBROAD != "Papio_ursinus"),]


#A) MAXIMUM LIFESPAN
Traits_top <- read.csv("../Fig1/Quantitative_extremes", header=T, sep="\t")
LIFESPAN_DATASET <- Traits_top[which(Traits_top$variable == "Maximum.lifespan..y."),]
#sort(primates_tree$tip.label)
#sort(LIFESPAN_DATASET$label[!is.na(LIFESPAN_DATASET$label)])
pruned_primates.tree<-drop.tip(primates_tree ,primates_tree$tip.label[-match(LIFESPAN_DATASET$label[!is.na(LIFESPAN_DATASET$label)], 
                                                                             primates_tree$tip.label)])

#pruned_primates.tree <- force.ultrametric(pruned_primates.tree)
lifespan <- LIFESPAN_DATASET$value
names(lifespan) <- LIFESPAN_DATASET$label


#ANCESTRAL RECONSTRUCTION
asr.cont <- anc.Bayes(pruned_primates.tree, lifespan)
# define vector with node numbers as character strings
nodes <- c("92","109","119","153")

# add continuous ASR posterior distributions for selected nodes
pd <- asr.cont$mcmc[,nodes]
traits_palette <- readRDS("traits_palette")
colfunc_ANC <- colorRampPalette(c(traits_palette[3], "white"))

library(ggimage)
pdplots <- apply(pd, 2, function(y) {
  ggplot(data.frame(y=y), aes(y, fill=..x..)) +
    geom_histogram( binwidth=0.01) +
    xlim(range(unlist(pd))) +
    scale_fill_gradient(low=colfunc_ANC(10)[c(7)], high=colfunc_ANC(10)[c(1)]) +
    ggtree::theme_inset()
})


df_pdplots <- tibble::tibble(node=c(92,109,119,153), 
                             pdplots=pdplots)

library(ggtree)
p12 <- ggtree(force.ultrametric(pruned_primates.tree),
              layout='circular')

#Annotate alternative Tops and Bottoms

Global_sel <- read.csv("Maximum_lifespan__y__config.txt",
                       sep="\t",
                       header=F)

AlterA <- read.csv("Maximum.lifespan..y._altern1_config.txt",
                   sep="\t",
                   header=F)

AlterB <- read.csv("Maximum.lifespan..y._altern2_config.txt",
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


#ELABORATE DF
species_df <- p12$data[which(!is.na(p12$data$label)),]

SPECIES_DF <- species_df %>% select(label, selGlobalValues,
                                    selAValues, selBValues) %>%
  mutate(value=abs(rnorm(90)))
colnames(SPECIES_DF) <- c("label", "Global", "A", "B", "value")
SPECIES_DF <- as.data.frame(SPECIES_DF)
rownames(SPECIES_DF) <- SPECIES_DF$label


#ELABORATE PLOT
p11 <- ggtree(force.ultrametric(pruned_primates.tree),
              layout='circular',
              size=1)

lifespan_df <- as.data.frame(cbind(names(lifespan), 
                                   as.numeric(as.character(lifespan))))
colnames(lifespan_df) <- c("label", "lifespan")
p12 <- p11 %<+% lifespan_df

p12$data$width <- ifelse(is.na(p12$data$lifespan), 0.1, 0.4)


#ADD FAMILIES

#add families to tree

families <- c()
for(i in 1:length(p12$data$label)){
  if (as.character(p12$data$label[i]) %in% as.character(sequenced_traits$SpeciesBROAD)){
    family <- as.character(sequenced_traits$groups[sequenced_traits$SpeciesBROAD == as.character(p12$data$label[i])])
    families <- c(families, family)}
  else{
    family <- NA
    families <- c(families, family)
  }}

p12$data$Family <- families
primate_families_palette<- readRDS("../Fig1/primate_families_palette.rds")

#ADD FAMILIES SATURATED TO PHYLOGENY

desat <- function(cols, sat=0.5) {
  X <- diag(c(1, sat, 1)) %*% rgb2hsv(col2rgb(cols))
  hsv(X[1,], X[2,], X[3,])
}

primate_families_palette_des <- desat(primate_families_palette, 0.25)

#Colour as

x <- as_tibble(force.ultrametric(pruned_primates.tree))
td <- data.frame(node = nodeid(x, names(lifespan)),
                 trait = lifespan)
fit <- phytools::fastAnc(pruned_primates.tree, 
                         lifespan, vars=TRUE, CI=TRUE)
nd <- data.frame(node = names(fit$ace), trait = fit$ace)
d <- rbind(td, nd)
d$node <- as.numeric(d$node)

library(tidytree)
tree_formatted <- tidytree::full_join(x, d, by = 'node') %>% as.treedata()

#ADD COLOR FAMILIES

p12 <-ggtree(tree_formatted, aes(color=trait), layout = 'circular', 
             ladderize = TRUE)  + 
  scale_color_gradientn(colours=rev(c(colfunc_ANC(10)[c(1)], 
                                      colfunc_ANC(10)[c(5)],
                                      "grey40", 'black'))) +
  geom_cladelabel(node=120, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[3]) +
  geom_cladelabel(node=151, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[1]) +
  geom_cladelabel(node=64, label="", extend=0.5,
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[2]) +
  geom_cladelabel(node=99, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[6]) +
  geom_cladelabel(node=98, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[4]) +
  geom_cladelabel(node=109, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[7]) +
  geom_cladelabel(node=112, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[5]) +
  geom_cladelabel(node=118, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[8]) +
  geom_cladelabel(node=65, label="", extend = 0.5,
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[16]) +
  geom_cladelabel(node=163, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[13]) +
  geom_cladelabel(node=161, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[14]) +
  geom_cladelabel(node=160, label="",
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[9])+
  geom_cladelabel(node=66, label="", extend=0.4,
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[10]) +
  geom_cladelabel(node=174, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[15])+
  geom_cladelabel(node=170, label="", 
                  offset=.001, align=TRUE,
                  barsize = 1.5, color=primate_families_palette_des[11]) 

library(ggtreeExtra)
OUTLIERS <- SPECIES_DF[which(SPECIES_DF$Global == 1),]
NON_OUTLIERS <- SPECIES_DF[which(!SPECIES_DF$Global == 1),]

p13 <- p12 + new_scale_color() + geom_fruit(data = NON_OUTLIERS,
                                            geom = geom_bar,
                                            mapping = aes(x = 1,
                                                          y =label,
                                                          color=Global,
                                                          shape=Global),
                                            offset=0.13,
                                            pwidth=0.002,
                                            stat = "identity",
                                            size=0.5,
                                            axis.params=list(
                                              axis="x", title="S",
                                              text.angle = 0,
                                              text.size = 0,
                                              line.size = 0,
                                              vjust = 0, 
                                              family="Helvetica", fontface="bold"
                                            )) +
  scale_color_manual(values = c(colfunc_ANC(10)[c(7)],
                                "grey40")) + 
  new_scale_fill() +
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


OUTLIERS <- SPECIES_DF[which(SPECIES_DF$A == 1),]
NON_OUTLIERS <- SPECIES_DF[which(!SPECIES_DF$A == 1),]

library(ggnewscale)

p13$scales$scales <- lapply(p13$scales$scales, function(s) {
  if (!is.null(s$guide) && !is.atomic(s$guide) && !is.list(s$guide)) {
    s$guide <- "legend"
  }
  s
})

p14 <- p13 + ggnewscale::new_scale_color() + geom_fruit(data = NON_OUTLIERS,
                                                        geom = geom_bar,
                                                        mapping = aes(x = 1,
                                                                      y =label,
                                                                      color=A,
                                                                      shape=A),
                                                        offset=0.13,
                                                        pwidth=0.002,
                                                        stat = "identity",
                                                        size=0.5,
                                                        axis.params=list(
                                                          axis="x", title="A",
                                                          text.angle = 0,
                                                          text.size = 0,
                                                          line.size = 0,
                                                          vjust = 0, 
                                                          family="Helvetica", fontface="bold"
                                                        )) +
  scale_color_manual(values = c(colfunc_ANC(10)[c(7)],
                                "grey40")) + 
  new_scale_fill() +
  geom_fruit(data = OUTLIERS,
             geom = geom_point,
             mapping = aes(y =label,
                           fill=A,
                           shape=A),
             offset=0,
             pwidth=0.002,
             size=2
  ) +
  scale_shape_manual(values = c(22, 24,1)) +
  #scale_color_manual(values = "black") +
  scale_fill_manual(values = c(colfunc_ANC(10)[c(5)])) 




library(ggpp)
p15 <- p14 %<+% df_pdplots

p16 <- p15 +
  geom_plot(data=td_filter(!isTip), 
            mapping=aes(x=x,y=y, label=pdplots), 
            vp.width=0.1, vp.height=0.1, 
            hjust=0.5, vjust=0.2)


p16$scales$scales <- lapply(p16$scales$scales, function(s) {
  if (!is.null(s$guide) && !is.atomic(s$guide) && !is.list(s$guide)) {
    s$guide <- "legend"
  }
  s
})
OUTLIERS <- SPECIES_DF[which(SPECIES_DF$B == 1),]
NON_OUTLIERS <- SPECIES_DF[which(!SPECIES_DF$B == 1),]
p17 <- p16 + new_scale_color() + geom_fruit(data = NON_OUTLIERS,
                                            geom = geom_bar,
                                            mapping = aes(x = 1,
                                                          y =label,
                                                          color=B,
                                                          shape=B),
                                            offset=0.13,
                                            pwidth=0.002,
                                            stat = "identity",
                                            size=0.5,
                                            axis.params=list(
                                              axis="x", title="B",
                                              text.angle = 0,
                                              text.size = 0,
                                              line.size = 0,
                                              vjust = 0, 
                                              family="Helvetica", fontface="bold"
                                            )) +
  scale_color_manual(values = c(colfunc_ANC(10)[c(7)],
                                "grey40")) + 
  new_scale_fill() +
  geom_fruit(data = OUTLIERS,
             geom = geom_point,
             mapping = aes(y =label,
                           fill=B,
                           shape=B),
             offset=0,
             pwidth=0.002,
             size=2
  ) +
  scale_shape_manual(values = c(22, 24,1)) +
  #scale_color_manual(values = "black") +
  scale_fill_manual(values = c(colfunc_ANC(10)[c(5)])) +
  theme(legend.position = "none") +
  ggtitle(expression(bolditalic("Maximum lifespan"))) +
  theme(plot.title = element_text(size = 15, hjust = 0.5,
                                  color = colfunc_ANC(10)[c(3)]))

#MAXIMUM LIFESPAN

#p13_only_one <- p13 + scale_shape_manual(values = c(22, 24,1)) +
# #scale_color_manual(values = "black") +
#  scale_fill_manual(values = c(colfunc_ANC(10)[c(5)])) +
#  theme(legend.position = "none")

saveRDS(p17, "Fig5A2_JAN24")