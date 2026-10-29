####################################################################
##########(FIGURE1A) HETEROGENEITY CLUSTERS FIGURE##################
#######Create Heatmap of traits and clusters consistency to show####
######comparable sets of traits for further analyses ###############
################research topic######################################
####################################################################

#I) READ LIBRARIES AND TRAITS DATASET
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
library(psych)
library(network)
library(dendextend)
library(ape)
library(randomcoloR)
library(magrittr)
library(cluster)
library(cluster.datasets)
library(cowplot)
library(NbClust)
library(clValid)
library(ggfortify)
library(clustree)
library(dendextend)
library(factoextra)
library(corrplot)



all_traits <- read.csv("PRIMATES_PHENOMIC_DATASET.tsv", 
                       header=T, sep = "\t")


#II) SELECT DESIRED SPECIES

seq_species <- read.csv("names_FINAL.txt",
                        sep = "\t", header = TRUE)
colnames(seq_species)<- c("Tree_Species", "Genome_species", "Trait_species")
covered_species_all <- as.character(seq_species$Trait_species[!is.na(seq_species$Tree_Species)])
covered_species <- as.character(covered_species_all[!is.na(covered_species_all)])
sequenced_traits <- all_traits[which(all_traits$SpeciesBROAD %in% covered_species),]



##III) ADD FAMILY TAXONOMICAL INFORMATION
groups_list <- strsplit(as.character(sequenced_traits$GroupName), "_")
groups <- c()
for(i in 1:length(groups_list)){
  group <- paste(groups_list[[i]][1], "_", groups_list[[i]][2], sep="")
  groups[i] <- group}
groups


#IV) Read biological domains of traits
trait_domains <- read.csv("Traits_domains_sorted_FIG1.txt", sep = "\t")

#V) Binarize/Dichotomize the trait dataset table
#before calling MCA algorithm for factor analysis of traitsxspecies,
#let's transform the data into presence/absence of traits

binarized_traits <- sequenced_traits[,-which(colnames(sequenced_traits) %in% c("SpeciesBROAD",  "GroupName"))] %>%
  mutate_all(funs(ifelse(is.na(.), 0, 1)))

#VI) Study the clusters of homogeneity and consistency with FactorMineR

binarized_traits <- cbind(binarized_traits, groups)
binarized_traits <- binarized_traits %>%
  mutate_all(as.character)
rownames(binarized_traits) <- sequenced_traits[,1]
binarized_res.mca <- MCA(binarized_traits[,-ncol(binarized_traits)], 
                         graph=FALSE, ncp = 10)

#VII) Check for inflektion point with "elbow" clustering method and uik criteria
library(inflection)
info_cluster_elbow <- fviz_nbclust(t(binarized_traits[,-ncol(binarized_traits)]), 
                                   hcut, method = "wss", k.max = 90, print.summary = TRUE)  + 
  geom_vline(xintercept = 13, linetype = 2)

uik(as.numeric(as.character(info_cluster_elbow$data$clusters)),
    as.numeric(as.character(info_cluster_elbow$data$y)))



#VIII) Trait clusters description
trans_binarized_res.mca <- MCA(t(binarized_traits[,-ncol(binarized_traits)]), 
                               graph=FALSE, ncp = 13)
trans_res.HCPC <- HCPC(trans_binarized_res.mca, nb.clust = 13, max = NULL, graph = TRUE)
fviz_cluster(trans_res.HCPC,
             show.clust.cent = TRUE, # Show cluster centers
             palette = distinctColorPalette(20),         # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             REPEL = TRUE,
             main = "Factor map"
)

#IX) Species cluster description
res.HCPC <- HCPC(binarized_res.mca, nb.clust = -1, min = 3, max = NULL, graph = FALSE)
fviz_cluster(res.HCPC,
             show.clust.cent = TRUE, # Show cluster centers
             palette = "jco", # Color palette see ?ggpubr::ggpar
             ggtheme = theme_minimal(),
             main = "Factor map"
)

#X) Prepare binarized matrix for plot, and sort the original matrix

binarized_traits_nongroup <- binarized_traits[,which(!colnames(binarized_traits) %in% c("SpeciesBROAD", "groups"))]
sorted_rows_bin <- binarized_traits_nongroup[order(rowSums(as.data.frame(sapply(binarized_traits_nongroup[,-ncol(binarized_traits_nongroup)], 
                                                                       as.numeric)))), ]

sorted_cols_bin <- sorted_rows_bin[, rev(order(colSums(as.data.frame(sapply(sorted_rows_bin,
                                                                            as.numeric)))))]


#XI) MELT THIS SORTED MATRIX
melted_tot_all <- sorted_cols_bin  %>%
  gather(var, val) 
melted_tot_all <- as.data.frame(cbind(rep(rownames(sorted_cols_bin), length(colnames(sorted_cols_bin))), 
                                      melted_tot_all))
colnames(melted_tot_all) <- c("species", "traits", "value")



#XII) Prepare dendogram, both for species and traits
#Reorder dendogram branches, both for species and for traits

sequenced_traits$groups <- groups

unloadNamespace("dendextend")
attachNamespace("dendextend")

colorsRainbow <- rainbow(4, alpha = 1)
dendogram_species <- as.dendrogram(res.HCPC$call$t$tree) %>% 
  set("branches_lwd", 0.7) %>% # hang the leaves (at some height)
  raise.dendrogram(0.02) %>%
  set("labels_cex", 10) %>%
  rotate(as.character(rownames(sorted_cols_bin)))
ggdend_species <- as.ggdend(dendogram_species)
dendogram_traits <- as.dendrogram(trans_res.HCPC$call$t$tree) %>% 
  set("branches_k_color", k = 13) %>% 
  set("branches_lwd", 0.7) 

#FLIP BRANCHES TRAITS

#View(as.data.frame(cbind(trans_res.HCPC$data.clust$clust,
 #                   rownames(trans_res.HCPC$data.clust))))
coverage_traits_clustered <- as.data.frame(cbind(rownames(trans_res.HCPC$data.clust),
                                                 as.numeric(as.character(trans_res.HCPC$data.clust$clust))))

colnames(coverage_traits_clustered) <- c("Trait", "Cluster")
coverage_traits_clustered <- coverage_traits_clustered[order(as.numeric(as.character(coverage_traits_clustered$Cluster)),
                                                             decreasing = TRUE),]

coverage_traits_clustered$Trait
dd2 <- rotate(dendogram_traits, coverage_traits_clustered$Trait)
#dend <- click_rotate(dendogram_traits)

# keep playing with the rotation:
#while (TRUE) dend <- click_rotate(dend)
# the same as
#dend <- click_rotate(dend, continue = TRUE)


#dend <- dend %>% 
 # raise.dendrogram(1) %>%
  #hang.dendrogram(hang_height = 0.01)

#dend %>% plot(main = "plot")
#ggdend_traits <- as.ggdend(dend)

#sort dend species
sorted_dend_species <- dendogram_species %>% 
  rotate(as.character(rownames(sorted_cols_bin)))
ggdend_species <- as.ggdend(sorted_dend_species)

#sort dend traits
#%>% rotate(as.factor(colnames(sorted_cols_bin)))
sorted_dend_traits <- dd2
ggdend_traits <- as.ggdend(sorted_dend_traits)



#XIII) Create dataset for species phylogeny

#levels(sequenced_traits$SpeciesBROAD) <- c(levels(sequenced_traits$SpeciesBROAD), "Tarsius_bancanus", "Leontocebus_fuscicolis",
#                                           "Plecturocebus_hoffmansi", "Lagothrix_lagothricha") 
#sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Cephalopachus_bancanus"] <- "Tarsius_bancanus"
#sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Leontocebus_fuscicollis"] <- "Leontocebus_fuscicolis"
#sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Plecturocebus_hoffmannsi"] <- "Plecturocebus_hoffmansi"
#sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Lagothrix_lagotricha"] <- "Lagothrix_lagothricha"


tree <- read.tree("FINAL_TREE_BL_WITHOUTNODELAB.newick")
library(ggtree)
library(phytools)

species_match <- match(sequenced_traits$SpeciesBROAD,
                       tree$tip.label)  

pruned_primates.tree<-drop.tip(tree,tree$tip.label[-species_match])
pruned_primates.tree <- force.ultrametric(pruned_primates.tree)
p <- ggtree(pruned_primates.tree)


families <- c()
for(i in 1:length(p$data$label)){
  if (as.character(p$data$label[i]) %in% as.character(sequenced_traits$SpeciesBROAD)){
    family <- as.character(sequenced_traits$groups[sequenced_traits$SpeciesBROAD == as.character(p$data$label[i])])
    families <- c(families, family)}
  else{
    family <- NA
    families <- c(families, family)
  }}

p$data$Family <- families
primate_families_palette<- readRDS("primate_families_palette.rds")

p1 <- p + geom_tippoint(aes(color=Family), size=1.5,
                        show.legend = FALSE) +
  scale_color_manual(values = primate_families_palette)  +
  theme_tree() +
  xlim(0,0.1)
#+
  #geom_text(aes(label=node), hjust=-.3) +
  #geom_cladelabel(node=301, label="Cercopithecidae", 
   #               offset=.005, align=TRUE) + 
  #geom_cladelabel(node=385, label="Hylobatidae", 
   #               offset=.005, align=TRUE) +
  #geom_cladelabel(node=380, label="Hominidae", 
   #               offset=.005, align=TRUE) +
  #geom_cladelabel(node=234, label="Callitrichidae", 
   #               offset=.005, align=TRUE) +
  #geom_cladelabel(node=231, label="Aotidae", 
   #               offset=.005, align=TRUE) +
  #geom_cladelabel(node=255, label="Cebidae", 
   #               offset=.005, align=TRUE) +
  #geom_cladelabel(node=264, label="Atelidae", 
   #               offset=.005, align=TRUE) +
  #geom_cladelabel(node=276, label="Pitheciidae", 
   #               offset=.005, align=TRUE) +
  #geom_cladelabel(node=394, label="Tarsiidae", 
   #               offset=.005, align=TRUE) +
  #geom_cladelabel(node=418, label="Lemuridae", 
   #               offset=.005, align=TRUE) +
  #geom_cladelabel(node=402, label="Indriidae", 
   #               offset=.005, align=TRUE) +
  #geom_cladelabel(node=414, label="Lepilemuridae", 
   #               offset=.005, align=TRUE, fontsize=2.5) +
  #geom_cladelabel(node=411, label="Cheirogaleidae", 
   #               offset=.005, align=TRUE, fontsize = 2.5) +
  #geom_cladelabel(node=174, label="Daubentonidae", 
   #               offset=.005, align=TRUE, fontsize=2.5) +
  #geom_cladelabel(node=439, label="Lorisidae", 
   #               offset=.005, align=TRUE) +
  #geom_cladelabel(node=435, label="Galagidae", 
   #               offset=.005, align=TRUE) 
  
#saveRDS(p1, "labelled_tree_JAN24")
#labelled_TREE <- readRDS("labelled_tree_JAN24")
#p1 <- labelled_TREE
  

#XIV) Create dataset of traits dendogram
segment_data_traits <- with(
  segment(ggdend_traits), 
  data.frame(x = y, y = x, xend = yend, yend = xend, col = col))
#View(ggdend_traits)

# Use the dendrogram label data to position the gene (species)
species_pos_table <- with(
  p$data[which(!is.na(p$data$label)),], 
  data.frame(y_center = y, 
             species = as.character(p$data$label[!is.na(p$data$label)]), 
             height = 1))
#levels(species_pos_table$species) <- c(levels(species_pos_table$species), "Cephalopachus_bancanus", "Leontocebus_fuscicollis",
 #                                          "Plecturocebus_hoffmannsi", "Lagothrix_lagotricha") 
#species_pos_table$species[species_pos_table$species == "Tarsius_bancanus"] <- "Cephalopachus_bancanus"
#species_pos_table$species[species_pos_table$species == "Leontocebus_fuscicolis"] <- "Leontocebus_fuscicollis"
#species_pos_table$species[species_pos_table$species == "Plecturocebus_hoffmansi"] <- "Plecturocebus_hoffmannsi"
#species_pos_table$species[species_pos_table$species == "Lagothrix_lagothricha"] <- "Lagothrix_lagotricha"


# Use the dendrogram label data to position the samples (traits)
traits_pos_table <- with(
  ggdend_traits$labels, 
  data.frame(x_center = x, traits = as.character(label), width = 1))


#XV) Join both dendogram datasets
all_hmap <- melted_tot_all %>% 
  left_join(species_pos_table) %>%
  left_join(traits_pos_table)


#XVI) Limits for the vertical axes
species_axis_limits <- with(
  species_pos_table, 
  c(min(y_center - 0.5 * height), max(y_center + 0.5 * height))
) + 
  0.1 * c(-1, 1) # extra spacing: 0.1

#XVII) Limits for the horizontal axes
traits_axis_limits <- with(
  traits_pos_table,
  c(min(x_center - 0.5 * width), max(x_center + 0.5 * width))) + 
  0.1 * c(-1, 1) # extra spacing: 0.1


#XVIII) Create plot combining all information of absence/presence into heatmap

#cluster8_lines<-data.frame(x=c(12,33), y=c(199,0), 
 #                    xend=c(33,33), yend=c(199,199))


#trait_order <- all_hmap[,c("traits", "x_center")] %>%
 # distinct() %>%
#  arrange(x_center)
#View(trait_order)
#traits_order_with_cluster <- left_join(trait_order, coverage_traits_clustered, 
#          by=c("traits" = "Trait"))
#View(traits_order_with_cluster)
#View(traits_pos_table)

#species_axis_limits <- with(
 # species_pos_table, 
  #c(max(y_center - 0.5 * height), min(y_center + 0.5 * height))
#) + 
 # 0.1 * c(-1, 1) # extra spacing: 0.1

#XIX) Add families colors for dendogram of species
#we add colours for groups each variable phylogenetically

#primates_palette <- readRDS("primate_families_palette.rds")
#selected_groups <- c()
#for (i in 1:length(species_pos_table$species)){
 # selected_groups[i] <- as.character(binarized_traits$groups[binarized_traits$SpeciesBROAD == 
  #                                                             as.character(species_pos_table$species[i])])
  #species_pos_table$families[i] <- as.character(binarized_traits$groups[binarized_traits$SpeciesBROAD == 
   #                                                                       as.character(species_pos_table$species[i])])
#}
#names(primates_palette) <- levels(factor(binarized_traits$groups))
#selected_colours <- c()
#for (i in 1:length(selected_groups)){
 # selected_colours[i] <- primates_palette[names(primates_palette) == selected_groups[i]]
  #}

#XX) Create legend for species dendogram plot

library(grid)
#get_legend <- ggplot() +
#  geom_segment(data = dendrogram_segments, 
#             aes(x=x, y=y, xend=xend, yend=yend)) +
#  geom_point(data = dendrogram_ends ,
#               aes(x=x, y=y.x, xend=xend, yend=yend, col = families), shape=19, size=3) + # test aes is for plotly
#  scale_color_manual(values = primates_palette) +
#  geom_point(data=species_pos_table, aes(y=0, x=y_center, col = families
   #                                      ),alpha =0) +
  #geom_label(data=species_pos_table, aes(y=0, x=y_center, col = families, 
                                        # label = species)) +
#  ylim(c(0.1,0)) +
#  theme(legend.key=element_blank(),
#        legend.position = "top",
#        legend.title = element_blank(),
#        legend.direction = "vertical",
#        legend.spacing.y = unit(0.0005, 'cm'))
#legend <- cowplot::get_legend(get_legend)
#grid.draw(legend)

#XXI) Get plot for species dendogram

#list_Fams <- strsplit(as.character(dendrogram_ends$families),
#                         "_", fixed = TRUE)
#for (i in 1:length(list_Fams)){
#  dendrogram_ends$Fams[i] <- list_Fams[[i]][2]
#}

#list_Fams_names <- strsplit(as.character(names(primates_palette)),
#                      "_", fixed = TRUE)

#for (i in 1:length(list_Fams_names)){
#  names(primates_palette)[i] <- list_Fams_names[[i]][2]
#}


#dendrogram_ends$Fams <- factor(dendrogram_ends$Fams, levels=rev(c("Cercopithecidae", "Hylobatidae", "Hominidae" ,
#                                                                      "Callitrichidae", "Aotidae", "Cebidae", "Atelidae",
#                                                                      "Pitheciidae", "Tarsiidae", "Lemuridae", 
#                                                                      "Indriidae",   "Lepilemuridae", "Cheirogaleidae", 
#                                                                      "Daubentoniidae", "Lorisidae", "Galagidae")))
#sorted_primates_palette <- primates_palette[rev(levels(dendrogram_ends$Fams))]


#XXII) Add biological Trait domains information and create dendogram plot
#for traits

library(wesanderson)

traits_palette <- readRDS("traits_palette")
traits_palette <- c(traits_palette[1:4], "#daa18f")
saveRDS(traits_palette, "traits_palette")
list_Primary <- strsplit(as.character(trait_domains$Domain_Secondary_Domain),
                                          "_", fixed = TRUE)
for (i in 1:length(list_Primary)){
  trait_domains$Primary_domain[i] <- list_Primary[[i]][1]
}

#selected_groups <- c()
domains <- c()
for (i in 1:length(traits_pos_table$traits)){
  if (!traits_pos_table$traits[i] %in% trait_domains$Measurement){
    print(traits_pos_table$traits[i])
  }
  domain_curr <- as.character(trait_domains$Primary_domain[trait_domains$Measurement ==  as.character(traits_pos_table$traits[i])])
 domains <- c(domains, domain_curr)
  }

traits_pos_table$domains <- as.character(domains)
traits_pos_table$domains[traits_pos_table$domains == "Life"] <- "Life-history"
#read traits palette
plt_traits_dendr_with_legend <- ggplot(segment_data_traits) + 
  geom_segment(aes(x = y, y = x, xend = yend, yend = xend)) + 
  geom_bar(data=traits_pos_table, aes(y=0, x=x_center, fill = domains ,
                                      col = domains),
           size=1.8, stat="identity", width = 1.2) +
  scale_fill_manual(values = traits_palette[c(1,2,3,4,6)]) +
  scale_color_manual(values = traits_palette[c(1,2,3,4,6)]) +
  scale_y_continuous(position = "right") + 
  scale_x_continuous( limits = traits_axis_limits, 
                      breaks = traits_pos_table$x_center,
                      labels = traits_pos_table$traits,
                      position = "top", 
                      expand = c(0, 0)) +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size = rel(0.8), angle = 45, hjust = 0),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.text=element_text(size=10),
        legend.position = "top",
        title = element_blank()
  ) +  guides(color = guide_legend(override.aes = list(shape = 2))) 

legend <- get_legend(plt_traits_dendr_with_legend)
legend <- ggdraw(legend)

plt_traits_dendr <- ggplot(segment_data_traits) + 
  geom_segment(aes(x = y, y = x, xend = yend, yend = xend)) + 
  geom_bar(data=traits_pos_table, aes(y=0, x=x_center, fill = domains ,
                                      col = domains),
            size=1.8, stat="identity", width = 1.2) +
  scale_fill_manual(values = traits_palette[c(1,2,3,4,6)]) +
  scale_color_manual(values = traits_palette[c(1,2,3,4,6)]) +
  scale_y_continuous(position = "right") + 
  scale_x_continuous( limits = traits_axis_limits, 
                      breaks = traits_pos_table$x_center,
                      labels = traits_pos_table$traits,
                      position = "top", 
                      expand = c(0, 0)) +
  theme(panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size = rel(0.8), angle = 45, hjust = 0),
        axis.text.y = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        legend.position = "none",
        title = element_blank()
  ) +  guides(color = guide_legend(override.aes = list(shape = 2))) +
  theme(legend.text=element_text(size=10))

#READ CLUSTERS INFO

coverage_traits_clustered <- read.csv("Traits_clusters_FINAL_JAN23.tsv", sep = "\t")
#View(coverage_traits_clustered1)
#View(coverage_traits_clustered)
coverage_traits_clustered$Annot[coverage_traits_clustered$Cluster %in% c(7:13) ] <- "Widely distributed"
coverage_traits_clustered$y[coverage_traits_clustered$Cluster %in% c(7:13) ] <- 0.3
coverage_traits_clustered$Annot[coverage_traits_clustered$Cluster %in% c(2,3,5:6) ] <- "Sexual selection, Brain and Biometabolites"
coverage_traits_clustered$y[coverage_traits_clustered$Cluster %in% c(2,3,5:6) ] <- 0.2
coverage_traits_clustered$Annot[coverage_traits_clustered$Cluster %in% c(1,4) ] <- "Strepshirrines clades"
coverage_traits_clustered$y[coverage_traits_clustered$Cluster %in% c(1,4) ] <- 0.1

coverage_traits_clustered <- coverage_traits_clustered[order(match(coverage_traits_clustered$Cluster,
                                      c(13,12,11,10,9,8,7,6,5,3,2,4,1))), ]

coverage_traits_clustered$x <- c(1:263)

info <- coverage_traits_clustered %>%
  group_by(Annot) %>%
  summarize(xend = max(x), x = min(x),
            yend = max(y), y = min(y))


#View(traits_pos_table)

plt_labels <- ggplot(data = info) +
  geom_segment(aes(x = as.numeric(as.character(x)),
                   xend = as.numeric(as.character(xend)),
                   y=as.numeric(as.character(y)),
                   yend = as.numeric(as.character(y))))+
  annotate("text",
            x=30, y=0.4, label="WD", fontface="italic",
           size = 3) +
  annotate("text",
           x=125, y=0.3, label="ID (SS+Br+BMT)", fontface="italic",
           size = 3) +
  annotate("text",
           x=225, y=0.2, label="STR", fontface="italic",
           size = 3)  +
    theme_classic() +
  scale_x_continuous(expand = c(0, 0)) + 
  theme(axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.line=element_blank(),
        axis.ticks = element_blank()) +
  ylim(0,0.5)


#GET roman numbers and plot heatmap

coverage_traits_clustered$Cluster <- as.roman(coverage_traits_clustered$Cluster)

all_hmap_labelled <- left_join(all_hmap,coverage_traits_clustered,
                               by=c("traits" = "Trait"))

#View(all_hmap_labelled)

#Cluster 8 mins and maxs
min(as.numeric(as.character(all_hmap_labelled$x_center[as.character(all_hmap_labelled$Cluster) == "VIII"])))
max(as.numeric(as.character(all_hmap_labelled$x_center[as.character(all_hmap_labelled$Cluster) == "VIII"])))


breaks_roman_groups <- all_hmap_labelled %>%
  group_by(Cluster) %>%
  summarize(median_x = median(x))
plot_all_hmap <- ggplot(unique(all_hmap_labelled), 
                        aes(x = x_center, y = y_center, fill = value, 
                            height = height, width = width)) + 
  geom_tile() +
  geom_rect(data = all_hmap_labelled, aes(xmin=37,ymin=0.4,xmax=56, ymax=223),
           fill="grey", color="black", alpha = 0) +
  scale_x_continuous(limits = traits_axis_limits, 
                     expand = c(0, 0),
                     labels = breaks_roman_groups$Cluster,
                     breaks = breaks_roman_groups$median_x) + 
  # For the y axis, alternatively set the labels as: gene_position_table$gene
  scale_y_continuous(breaks = species_pos_table[, "y_center"], 
                     #labels = species_pos_table[, "species"],
                     limits = species_axis_limits) + 
  scale_fill_manual(values = c("white", "darkseagreen")) +
  geom_rect(data = all_hmap, aes(xmin=12,ymin=0,xmax=33, ymax=223),
            fill="grey", alpha=0.2) +
  theme_classic() +
  theme(axis.text.x = element_text(size=7,
                                   angle = 90), 
        # margin: top, right, bottom, and left
        strip.placement = "outside",
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(),
        axis.line.y = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "None",
        panel.spacing = unit(0, "lines"),
        strip.background = element_blank()) 



  

#XXIII) Create final plot merging all information obtained previously
library(cowplot)
library(ggpubr)

#blank <- grid.rect(gp=gpar(col="white"))
p2 <- ggdraw() +
  draw_plot(p1, x = 0, y = 0.17, width = .55, height = .52, hjust=0, vjust = 0) +
  draw_plot(plot_all_hmap, x = 0.3, y = 0.1, width = .7, height = .6, hjust=0, 
            vjust = 0) +
  draw_plot(plt_traits_dendr, x = 0.3, y = 0.66, width = .7, height = .35, hjust=0,
            vjust = 0)  +
  draw_plot(plt_labels, x = 0.3, y = 0, width = .7, height = .1, hjust=0,
            vjust = 0)
  

p <- legend
#ggdraw(p)
#saveRDS(p, "ggPHYLOleg")
saveRDS(p2, "ggPHYLO_JAN24")
