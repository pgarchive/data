######################################
##COMPONENTS STUDY FOR TRAITS GROUPS##
###, BASED ON MERGED MULTIPLE GROUPS##
###IN THE ANALYSIS OF CATEGORIES #####
###AND VARIABLES #####################
######################################

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
library(ape)
library(ggtree)
library(randomcoloR)
library(devtools)
library(phytools)
library(fastDummies)
options(scipen=0)



#A) READ LAST TRAITS FILE and metadata

all_traits <- read.csv("../../Figures/Fig1/PRIMATES_PHENOMIC_DATASET.tsv", header=T, sep = "\t")

#Now let's bind the Group info


#B) SELECT DESIRED SPECIES (sequenced_ones)


primates_palette <- readRDS("../../Figures/Fig1/primate_families_palette.rds")
seq_species <- read.csv("../../Figures/Fig1/names_FINAL.txt",
                        sep = "\t", header = TRUE)
colnames(seq_species)<- c("Tree_Species", "Genome_species", "Trait_species")
covered_species_all <- as.character(seq_species$Trait_species[!is.na(seq_species$Tree_Species)])
covered_species <- as.character(covered_species_all[!is.na(covered_species_all)])
sequenced_traits <- all_traits[which(all_traits$SpeciesBROAD %in% covered_species),]

groups_list <- strsplit(as.character(sequenced_traits$GroupName), "_")
groups <- c()
for(i in 1:length(groups_list)){
  group <- paste(groups_list[[i]][1], "_", groups_list[[i]][2], sep="")
  groups[i] <- group}
groups
sequenced_traits <- cbind(sequenced_traits, groups)
sequenced_traits <- sequenced_traits %>%
  mutate_at(vars(`IUCN_ER`), funs(case_when(
    . %in% c("5", "4", "3") ~ "1", 
    . %in% c("2", "1") ~ "0")))


#read redundant measurements

redundant <- read.csv("Redundant_measurements.tsv", header=FALSE, sep = "\t")

#D) BINARIZE THE MATRIX

binarized_traits <- sequenced_traits[,-which(colnames(sequenced_traits) %in% c("SpeciesBROAD",  "GroupName",
                                                                               "groups") |
                                               colnames(sequenced_traits) %in% redundant$V1)] %>%
  mutate_all(funs(ifelse(is.na(.), 0, 1)))


#E) STUDY THE CLUSTERS OF HOMOGENEITY AND CONSISTENCY

#binarized_traits <- cbind(binarized_traits, groups)

rownames(binarized_traits) <- sequenced_traits[,1]


##############################################
#A) ANALYSIS OF CORRELATIONS WITH CLUSTERS###
##############################################

#SELECT SPECIES AND TRAITS OF INTEREST
#for cluster 7 (brain), need to low down the nº of correlating species

coverage_traits_clustered <- read.csv("../../Figures/Fig1/Traits_clusters_FINAL_JAN23.tsv", sep = "\t")

traits_clust <- coverage_traits_clustered[which(coverage_traits_clustered$Cluster ==4),]

data_clust <- binarized_traits[,which(colnames(binarized_traits) %in% traits_clust$Trait)]

num_data_clust <- sapply(data_clust, as.numeric)
rownames(num_data_clust) <- rownames(data_clust)
##Subset as usual
best_species_clust <- rownames(num_data_clust[which(rowSums(num_data_clust) >= 0.90*ncol(num_data_clust)),])
best_species_clust
all_traits_clust <- binarized_traits[which(rownames(binarized_traits) %in% best_species_clust),]



#totals <- colSums(sapply(all_traits_clust, as.numeric))
best_traits_clust <- all_traits_clust[, which(colSums(all_traits_clust) >= 0.90*nrow(all_traits_clust))]
best_traits_clust <- colnames(best_traits_clust)

#MULTIVARIATE ANALYSIS
trues <- sequenced_traits$SpeciesBROAD %in% best_species_clust
#length(trues[trues==TRUE])

selected_subset <- sequenced_traits[which(sequenced_traits$SpeciesBROAD %in% best_species_clust),
                                    which(colnames(sequenced_traits) %in% best_traits_clust)] 

rownames(selected_subset) <- sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD %in% best_species_clust]
selected_subset <- selected_subset[,which(!colnames(selected_subset) %in% as.character(redundant$V1))]



library(missMDA)
no_missing <- selected_subset %>% drop_na()

res.famd <- FAMD(no_missing) #cluster6 | cluster7 | cluster 8

nrow(selected_subset)
#sequenced_traits[,1]
fviz_famd_var(res.famd, "quali.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)
options(ggrepel.max.overlaps = Inf)

p1 <- fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)

#View(res.famd$quanti.var$contrib)



#2) make correlations for quantitative sets of data

var_cor <- as.data.frame(res.famd$quanti.var)

primates_tree <- read.newick("../../Figures/Fig1/FINAL_TREE_BL_WITHOUTNODELAB.newick")

#CLUSTS
quanti_clust <- selected_subset[,which(colnames(selected_subset) %in% rownames(var_cor)
                                         & colnames(selected_subset) != "IUCN_ER"
                                       & colnames(selected_subset) != "R_Min_LitterSize")] %>% drop_na()
matched_rows <- match(rownames(quanti_clust),
                      primates_tree$tip.label)
matched_rows <- matched_rows[!is.na(matched_rows)]
clust.tree<-drop.tip(primates_tree,primates_tree$tip.label[-matched_rows])
quanti_clust_phylo <- quanti_clust[which(rownames(quanti_clust) %in% primates_tree$tip.label[matched_rows] ),] 

#imputed_DAT <- imputeFAMD(quanti_clust_phylo)
correl_struct <- phyl.vcv(as.matrix(quanti_clust_phylo),vcv(clust.tree),1)
correl <- cov2cor(correl_struct$R)
#correl_struct <- phyl.vcv(as.matrix(imputed_DAT$completeObs),vcv(clust.tree),1)
correl <- correl[rowSums(is.na(correl)) != ncol(correl)-1,]
correl <- correl[,colSums(is.na(correl)) != nrow(correl)]



#r.xy<-cov2cor(correl_struct$R)["BodyMass_kg","White.Blood.Cells"]
#t.xy<-r.xy*sqrt((Ntip(clust.tree)-2)/(1-r.xy^2))
#P.xy<-2*pt(abs(t.xy),df=Ntip(clust.tree)-2,lower.tail=F)

t.xy<-correl*sqrt((Ntip(clust.tree)-2)/(1-correl^2))
#View(t.xy)
P.xy<-2*pt(abs(t.xy),df=Ntip(clust.tree)-2,lower.tail=F)
P.xy

# Obtener los índices del triángulo superior (sin la diagonal)
upper.tri.idx <- which(upper.tri(P.xy), arr.ind = TRUE)

# Extraer p-values únicos (solo del triángulo superior)
p_values_unique <- P.xy[upper.tri(P.xy)]

# Aplicar corrección FDR
fdr_values <- p.adjust(p_values_unique, method = "fdr")

#View(p_values_unique)
# Crear dataframe con las coordenadas y los valores corregidos
df_results <- data.frame(
  X = rownames(P.xy)[upper.tri.idx[, 1]],
  Y = colnames(P.xy)[upper.tri.idx[, 2]],
  P_value = p_values_unique,
  FDR = fdr_values,
  Index = seq_along(fdr_values) # Genera un índice único por fila
)

# Duplicar los valores para reflejar la simetría en el dataframe final
df_final <- rbind(
  df_results, 
  df_results %>% rename(X = Y, Y = X) # Intercambia X e Y para reflejar el triángulo inferior
)

# Agregar filas con la diagonal como NA
diag_data <- data.frame(
  X = rownames(P.xy),
  Y = rownames(P.xy),
  P_value = NA,
  FDR = NA,
  Index = NA
)

# Unir todo en un solo dataframe
df_final <- rbind(df_final, diag_data) %>% arrange(FDR)
colnames(df_final) <- c("Var1", "Var2", "p-value", "fdr")

#MELTED <- melt(P.xy)
#MELTED$fdr <- p.adjust(MELTED$value)

FDR <- df_final[,c(1,2,4)]
#FDR <- MELTED[,c(1,2,4)]
FDR_MATRIX <- acast(FDR, Var1~Var2, value.var='fdr', fill=0)


library(GGally)
#Correlation matrices

#View(correl)
p1 <- ggcorr(data = NULL, 
       cor_matrix = correl,
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_size = 2,
       label_round = 2,
       color = "black",
       low = "darkred",
       mid = "#EEEEEE", high = "darkgreen",
       layout.exp = 0.4,
       size=3) 


#p-value matrices
DF <-melt(FDR_MATRIX)
library(ComplexHeatmap)

ggcorr(data = NULL, 
       cor_matrix = FDR_MATRIX,
       midpoint = 0.5,
       nbreaks = 6,
       label = TRUE,
       label_size = 3,
     #  layout.exp	= 2,
       label_round = 3,
       limits = FALSE,
       color = "black",
       low = "lightblue",
       mid = "white", high = "white",
     layout.exp=3,
       size=3,
     hjust = 1) +
  anno_zoom(align_to = DF$Var1)


#PRUEBA TOTAL:)
library(devtools)

layout.exp=2


library("ggcor")

trait_domains <- read.csv("../../Figures/Fig1/Traits_domains_sorted_FIG1.txt", sep = "\t")
traits_palette <- readRDS("../../Figures/Fig1/traits_palette")

DF_domains <- left_join(DF, trait_domains, by=c("Var1"= "Measurement"))

SPLITTED <- strsplit(as.character(DF_domains$Domain_Secondary_Domain), split="_", fixed = TRUE)
names <- c()
for (i in 1:length(SPLITTED)){
  name <- SPLITTED[[i]][1]
  if (name == "Life"){
  name <- gsub(name, "Life", "Life-history")
  }
  names <- c(names, name)
}

DF_domains$Domains <- names


DF_summ <- DF_domains %>%
  distinct(Var1, Domains)

DF_summ <- as.data.frame(cbind(as.character(DF_summ$Domains),
                               as.character(DF_summ$Var1)))
rownames(DF_summ) <- DF_summ$Var1

SORTED_MATRIX <- FDR_MATRIX[order(apply(FDR_MATRIX, 1, function(row) sum(row)), decreasing = F), ]
SORTED_MATRIX_COLUMNS <- SORTED_MATRIX[,rownames(SORTED_MATRIX)]


DF <-melt(FDR_MATRIX)

# Ordenar por la media de FDR en vez de la suma, ignorando NA
row_order <- order(rowMeans(FDR_MATRIX, na.rm = TRUE), decreasing = FALSE)

# Aplicar el orden a las filas y columnas
SORTED_MATRIX <- FDR_MATRIX[row_order, row_order]
SORTED_MATRIX_COLUMNS <- SORTED_MATRIX[,rownames(SORTED_MATRIX)]


trait_domains <- read.csv("Traits_domains_SPLIT_FEB22.txt", sep = "\t")


get_upper_tri <- function(CorMat){
  CorMat[lower.tri(CorMat)]<- NA
  return(CorMat)
}


upper_tri <- get_upper_tri(SORTED_MATRIX_COLUMNS)
diag(upper_tri) <- NA
SORTED_DF <- melt(upper_tri) %>%
  drop_na()

#COLORS DOMAINS
#View(DF_summ)
dmc <- data.frame(x = factor(DF_summ$V2), colour = factor(DF_summ$V1))
rownames(dmc) <- dmc$x
SORTED_DMC <- dmc[match(rownames(SORTED_MATRIX_COLUMNS), rownames(dmc)), ]
SORTED_DMC$x <- factor(SORTED_DMC$x, levels=rev(rownames(SORTED_MATRIX_COLUMNS)))

SORTED_MERGED <- left_join(SORTED_DF, SORTED_DMC, by=c("Var2" = "x"))
SORTED_MERGED2 <- left_join(SORTED_MERGED, SORTED_DMC, by=c("Var1" = "x"))

#View(SORTED_MERGED)
#GENERATE PLOT
#SORTED_DF$Var2
p1 <- ggplot(data = SORTED_MERGED2) +
  labs(x = NULL, y = NULL) + 
  geom_tile(mapping = aes(Var1, Var2, 
                          fill = value),
            drop=FALSE) +
  scale_y_discrete(
    limits = rev
  ) +
  labs(fill = "FDR p-value") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 9, angle = 30,
                                   hjust =1),
        axis.text.y = element_text(size = 9),
        axis.ticks.x =  element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        legend.position = "right",
        text = element_text(family = "Helvetica"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)) +
  new_scale_fill() +
  labs(x = NULL, y = NULL) +
  geom_tile(
    mapping = aes(x = -0.45, y = Var2, fill = colour.x, width = 0.5)
  ) +
  
  geom_tile(
    aes(y = -0.45, x = Var1, fill = colour.y, width = 1)
  ) +
  labs(fill = "Primary domain") +

  scale_fill_manual(values = traits_palette[c(1,2,4,5)]) +
  theme(legend.position = "right")
  

saveRDS(p1, "FigureS15")


















#####################################################
#3) MAKE QUALITATIVE CORRELATIONS WITH PHYLIP (THRESHOLD MODEL)##
#####################################################

sessionInfo()
library(mvMORPH)

#CLUSTER7-13
library(Rphylip)
primates_tree <- read.newick("/home/upf/Documents/LInux_DIC_2020/Documents/PhD_EvoGenomics/3rd_year/GENOME-PHENOME-PAPER_NOV2021/data/Tree/FINAL_TREE_BL.newick")
metadata_numeric <- read.csv("../../results/Metadata_numeric_traits_2.0.tsv", header=T, sep = "\t")


select_dummy <- dummy_cols(selected_subset, 
                           select_columns = c("DielActivity", "TrophicGuild",
                                              "Locomotion", "NatalDispersal",
                                              "Habitat", "IUCN_ER",
                                              "Seasonal_breeders",
                                              "Mating_system", "Social_organization",
                                              "Colour_vision", "Paternal_care", 
                                              "Female_Ranging_patterns", "Infanticide",
                                              "Monopolization", "StableGroups",
                                              "SocialMatingSystem"))
#select_dummy<- select_dummy[,-c(1:4)]
#select_dummy <- select_dummy[,-c(4)]
#select_dummy <- select_dummy[,-c(5:10)]
#select_dummy <- select_dummy[,-c(6:29)]

colnames(select_dummy)
traits_filtered <- c("IUCN_ER", "DielActivity_CA",
                     "DielActivity_N", "DielActivity_D",
                     "TrophicGuild_Folivore", "TrophicGuild_Folivore_frugivore",                              
                     "TrophicGuild_Frugivore", "TrophicGuild_Gummivore",
                     "TrophicGuild_Insectivore", "TrophicGuild_Omnivore",
                     "Locomotion_AR", "Locomotion_BOTH","Locomotion_T",
                     "NatalDispersal_Bisexual", "NatalDispersal_Sexbiased",
                     "Habitat_Forest", "Habitat_Forest/Savanna",
                     "Habitat_Forest/Savanna/Shrubland", 
                     "Habitat_Forest/Shrubland",
                     "Habitat_Forest/Wetlands",
                     "Seasonal_breeders_Yes",
                     "Mating_system_Monogamous",
"Mating_system_Polygynandrous",
"Mating_system_Polygynandrous/Monogamous",
"Mating_system_Polygynous", 
"Mating_system_Polygynous/Monogamous",                          
"Mating_system_Polygynous/Polygynandrous",
"Social_organization_Multimale",                                
"Social_organization_Pair",                                     
"Social_organization_Pair/Multimale",                           
"Social_organization_Pair/Unimale",                             
"Social_organization_Pair/Unimale/Multimale",                   
"Social_organization_Solitary",                                 
"Social_organization_Unimale",                                  
"Social_organization_Unimale/Multimale", 
"Colour_vision_Dichromatic",                                    
"Colour_vision_Monochromatic",                                  
"Colour_vision_Polymorphic",  
"Paternal_care_Present",
"Colour_vision_Trichromatic",
"Female_Ranging_patterns_Discrete",
"Seasonal_breeders_Yes")

selected_qual_dat <- select_dummy[, which(colnames(select_dummy) %in%
                                       traits_filtered)] %>% select(-contains("_NA"))


rownames(selected_qual_dat) <- rownames(selected_subset)


write.table(selected_qual_dat, "qual_dummy",
            quote = FALSE, sep = "\t", row.names = TRUE,
            col.names=NA)


#SELECT SPECIES FROM TREE
species_to_keep <- match(rownames(selected_qual_dat), primates_tree$tip.label)
species_to_keep <- species_to_keep[!is.na(species_to_keep)]
pruned.tree<-drop.tip(primates_tree,primates_tree$tip.label[-species_to_keep])
write.tree(pruned.tree, file="pruned_tree.newick")
select_dummy <- select_dummy[rownames(select_dummy) %in% pruned.tree$tip.label,]
X <- select_dummy

library(parallel)
library(foreach)
cores=detectCores()
cl <- makeCluster(4, outfile = "prueba.txt")
clusterEvalQ(cl, library("Rphylip"))
clusterExport(cl, c("threshBayes"), 
              envir=environment())
doParallel::registerDoParallel(cl = cl)

corrs <- c()
corrs_final <- foreach(i=1:(ncol(X)-1)) %dopar% sample_fun(i, corrs = corrs)

#n.obs <- sapply(corrs_final, length)
#seq.max <- seq_len(max(n.obs))

mat <- read.csv("all_cors_pasted_JAN23.txt", sep="\t", header = FALSE)
df2 <- mutate_all(mat, function(x) as.numeric(as.character(x)))


#View(mat)
rownames(selected_qual_dat)
rownames(df2)
length(colnames(selected_qual_dat))
length(rownames(selected_qual_dat))
colnames(select_dummy)
rownames(df2) <- colnames(selected_qual_dat)
colnames(df2) <- colnames(selected_qual_dat)
#mat <- do.call(rbind, corrs_final)
#rownames(mat) <- colnames(X)[1:(ncol(X)-1)]
#colnames(mat) <- colnames(X)[2:(ncol(X))]

#View(df2)
library(GGally)
diag(df2)<-1
df2[upper.tri(df2)] <- t(df2)[upper.tri(t(df2))]
df2 <- df2[which(!rownames(df2) %in% c("TrophicGuild_Insectivore")),
    which(!colnames(df2) %in% c("TrophicGuild_Insectivore"))]
#diag(m) <- diag
ggcorr(data = NULL, 
       cor_matrix = df2,
       nbreaks = 6,
       hjust = 0.8,
       label = TRUE,
       label_round = 2,
       label_size = 2,
       color = "black",
       low = "darkred",
       mid = "#EEEEEE", high = "darkgreen",
       layout.exp = 0.4,
       size=3) 

#View(df2)
#FDR qualitative analysis
#View(t.xy)
t.xy<-df2*sqrt((Ntip(pruned.tree)-2)/(1-df2^2))
#View(t.xy)
t.xy <- as.matrix(t.xy)
t.xy <- mutate_all(t.xy, function(x) as.numeric(as.character(t.xy)))
P.xy<-2*pt(abs(t.xy),df=Ntip(pruned.tree)-2,lower.tail=F)
P.xy

MELTED <- melt(P.xy)
MELTED$fdr <- p.adjust(MELTED$value)
FDR <- MELTED[,c(1,2,4)]
FDR_MATRIX <- acast(FDR, Var1~Var2, value.var='fdr', fill=0)
DF <-melt(FDR_MATRIX)

SORTED_MATRIX <- FDR_MATRIX[order(apply(FDR_MATRIX, 1, function(row) sum(row)), decreasing = F), ]
SORTED_MATRIX_COLUMNS <- SORTED_MATRIX[,rownames(SORTED_MATRIX)]

trait_domains <- read.csv("../FIGURES_PAPER_2022/Figure3/Traits_domains_SPLIT_FEB22.txt", sep = "\t")

qualitative_traits <- c("AgilityScore", "BodyFlesh",
                        "BodyHair", "BodyPelageDichrom", 
                        "Colour_vision",
                        "DielActivity",
                        "GenitalHair",
                        "GenitalSkin",
                        "Habitat",
                        "HeadFlesh",
                        "HeadHair",
                        "HeadSkin",
                        "IUCN_Pop_Trend",
                        "Locomotion",
                        "L_Pct_ND_InfMort_lt30d",
                        "Mating_System",
                        "Paternal_care",
                        "Female_Ranging_Patterns",
                        "Positional.behaviour",
                        "Prehensile_tail",
                        "R_Max_LitterSize",
                        "R_MostCommon_LitterSize_1",
                        "R_MostCommon_LitterSize_2",
                        "R_MostCommon_LitterSize_3",
                        "R_Pattern_Breeding",
                        "Seasonal_breeders",
                        "SexualOrnamentDimorphism",
                        "SexualOrnamentDimorphism..ExcludingDichromatism.",
                        "SocialMatingSystem",
                        "Social_organization",
                        "SocialSystem",
                        "TrophicGuild",
                        "Female_Ranging_patterns",
                        "StableGroups")

DF$TRAIT <- DF$Var1
levels(DF$TRAIT) <- c(levels(DF$TRAIT) ,
                                        c("Colour_vision", "Social_organization",
                                          "Mating_system", "DielActivity",
                                          "TrophicGuild", "Locomotion",
                                          "Habitat", "SocialMatingSystem", "Monopolization", "NatalDispersal",
                                          "Paternal_care", "Female_Ranging_patterns",
                                          "Seasonal_breeders"))

DF$TRAIT[startsWith(as.character(DF$Var1), "Colour_vision_")] <- "Colour_vision"
DF$TRAIT[startsWith(as.character(DF$Var1), "Social_organization_")] <- "Social_organization"
DF$TRAIT[startsWith(as.character(DF$Var1), "Mating_system_")] <- "Mating_system"
DF$TRAIT[startsWith(as.character(DF$Var1), "DielActivity_")] <- "DielActivity"
DF$TRAIT[startsWith(as.character(DF$Var1), "TrophicGuild_")] <- "TrophicGuild"
DF$TRAIT[startsWith(as.character(DF$Var1), "Locomotion_")] <- "Locomotion"
DF$TRAIT[startsWith(as.character(DF$Var1), "Habitat_")] <- "Habitat"
DF$TRAIT[startsWith(as.character(DF$Var1), "SocialMatingSystem_")] <- "SocialMatingSystem"
DF$TRAIT[startsWith(as.character(DF$Var1), "Monopolization_")] <- "Monopolization"
DF$TRAIT[startsWith(as.character(DF$Var1), "NatalDispersal_")] <- "NatalDispersal"
DF$TRAIT[startsWith(as.character(DF$Var1), "Paternal_care_")] <- "Paternal_care"
DF$TRAIT[startsWith(as.character(DF$Var1), "Female_Ranging_patterns_")] <- "Female_Ranging_patterns"
DF$TRAIT[startsWith(as.character(DF$Var1), "Seasonal_breeders_")] <- "Seasonal_breeders"




DF_domains <- dplyr::left_join(DF,
                                                 trait_domains, by=c("TRAIT" = "Measurement"))

DF_summ <- DF_domains %>%
  distinct(Var1, Domain)

DF_summ_SUMM <- as.data.frame(cbind(as.character(DF_summ$Domain),
                                    as.character(DF_summ$Var1)))


get_upper_tri <- function(CorMat){
  CorMat[lower.tri(CorMat)]<- NA
  return(CorMat)
}

#quickcor(SORTED_MATRIX_COLUMNS, type = "upper", show.diag=FALSE,
 #        cor.test = FALSE) +

upper_tri <- get_upper_tri(SORTED_MATRIX_COLUMNS)
diag(upper_tri) <- NA
SORTED_DF <- melt(upper_tri) %>%
  drop_na()

#COLORS DOMAINS
dmc <- data.frame(x = factor(DF_summ_SUMM$V2), colour = factor(DF_summ_SUMM$V1))
rownames(dmc) <- dmc$x
SORTED_DMC <- dmc[match(rownames(SORTED_MATRIX_COLUMNS), rownames(dmc)), ]
SORTED_DMC$x <- factor(SORTED_DMC$x, levels=rev(rownames(SORTED_MATRIX_COLUMNS)))

SORTED_MERGED <- left_join(SORTED_DF, SORTED_DMC, by=c("Var2" = "x"))
SORTED_MERGED2 <- left_join(SORTED_MERGED, SORTED_DMC, by=c("Var1" = "x"))

#View(SORTED_MERGED)
#GENERATE PLOT
#SORTED_DF$Var2
ggplot(data = SORTED_MERGED2) +
  labs(x = NULL, y = NULL) + 
  geom_tile(mapping = aes(Var1, Var2, 
                          fill = value),
            drop=FALSE) +
  scale_y_discrete(
    limits = rev
  ) +
  labs(fill = "FDR p-value") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 9, angle = 30,
                                   hjust =1),
        axis.text.y = element_text(size = 9),
        axis.ticks.x =  element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(0.4, 'cm'),
        text = element_text(family = "Helvetica"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)) +
  new_scale_fill() +
  labs(x = NULL, y = NULL) +
  geom_tile(
    mapping = aes(x = -0.45, y = Var2, fill = colour.x, width = 0.5)
  ) +

  geom_tile(
    aes(y = -0.45, x = Var1, fill = colour.y, width = 1)
  ) +
  labs(fill = "Primary domain") +
  
  scale_fill_manual(values = traits_palette[c(1,2,6)])































#stop cluster
stopCluster(cl)


sample_fun <- function(i, corrs){
  corrs <- c(corrs, rep(NA, i-1))
  library(Rphylip)
  st <- i + 1
  for (j in st:ncol(X)){
    #k <- k + 1
    print(j)
    subset <- X[,(c(i,j))]
    complete_subset <- subset[complete.cases(subset),]
    species_to_keep <- match(rownames(complete_subset), pruned.tree$tip.label)
    species_to_keep <- species_to_keep[!is.na(species_to_keep)]
    new_pruned.tree<-drop.tip(pruned.tree,pruned.tree$tip.label[-species_to_keep])
    complete_subset  <- complete_subset[rownames(complete_subset ) %in% new_pruned.tree$tip.label,]
    if (class(try(threshBayes(new_pruned.tree, complete_subset))) == "try-error"){
      corr <- NA
      corrs <- c(corrs, corr)
    }
    else{
      if (typeof(complete_subset[,1]) == "double"){
        type1 = "continuous"
      }
      else{
        type1 = "discrete"
      }
      if (typeof(complete_subset[,2]) == "double"){
        type2 = "continuous"
      }
      else{
        type2 = "discrete"
      }
      fit <- threshBayes(new_pruned.tree, complete_subset, types=c(type1, type2))
      corr <- mean(fit$par[,"r"])
      #fits[[k]] <- fit
      corrs <- c(corrs, corr)
    }
    
  }
  return(corrs)
}


#parallelize process traits

library(parallel)

n_cores <- detectCores(logical = TRUE)
cl <- makeCluster(4, outfile = "prueba.txt")

clusterEvalQ(cl, c(library(caper)))
clusterExport(cl, c("pgls"), 
              envir=environment())

init <- Sys.time()

cor_list <- parLapply(cl, as.data.frame(permulations_df[,1:1000]), species_cancer = species_cancer,
                      ldf_subset = ldf_subset, pruned_tree = pruned_tree, 
                      PGLS_fun)
compare$runtime[4] <- Sys.time() - init

stopCluster(cl)









#MULTIVARIATE ANALYSIS
res.famd <- FAMD(selected_subset)


fviz_famd_var(res.famd, "quali.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07")
)

options(ggrepel.max.overlaps = Inf)
fviz_famd_var(res.famd, "quanti.var", col.var = "contrib", 
              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
              repel = TRUE)


var_cor <- as.data.frame(var_quanti$var$cor)

cluster_biomet <- log(quanti+0.001)
write.table(cluster_biomet, "Cluster_biometabolites", quote=FALSE, sep="\t")


library(Hmisc)
library(corrplot)
Cluster_biometabolites <- read.csv("Cluster_biometabolites", sep ="\t")
results <- cor(Cluster_biometabolites,
           use="pairwise.complete.obs")
corrplot(results,
         tl.col = "black")




"""
redundant_morphometric1_traits <- rownames(var_quanti_cor[as.numeric(as.character(var_quanti_cor$Dim.1)) >= 0.75,])
redundant_dimorphism_traits <- rownames(var_quanti_cor[as.numeric(as.character(var_quanti_cor$Dim.1)) >= 0.75,])
redundant_morphometric2_traits <- rownames(var_quanti_cor[as.numeric(as.character(var_quanti_cor$Dim.1)) >= 0.75,])
redundant_clust4_traits <- rownames(var_quanti_cor[as.numeric(as.character(var_quanti_cor$Dim.1)) >= 0.75,])
redundant_biometabolite_traits <- rownames(var_quanti_cor[as.numeric(as.character(var_quanti_cor$Dim.1)) >= 0.75,])
redudant_brain_traits <- rownames(var_quanti_cor[as.numeric(as.character(var_quanti_cor$Dim.1)) >= 0.75,])
redudant_LHT_traits <- rownames(var_quanti_cor[as.numeric(as.character(var_quanti_cor$Dim.1)) >= 0.75,])
redudant_ecological1_traits <- rownames(var_quanti_cor[as.numeric(as.character(var_quanti_cor$Dim.1)) >= 0.75,])
redudant_LHT2_traits <-rownames(var_quanti_cor[as.numeric(as.character(var_quanti_cor$Dim.1)) >= 0.75,])
redudant_ecological2_traits <-rownames(var_quanti_cor[as.numeric(as.character(var_quanti_cor$Dim.1)) >= 0.75,])
"""

#LOAD METADATA TO KEEP ONLY categorical variables

metadata_categorical <- read.csv("../../results/Metadata/Metadata_categorical_traits_FINAL.tsv", header=T, sep = "\t")
metadata_binary <- read.csv("../../results/Metadata/Metadata_binary_traits_FINAL.tsv", header=T, sep = "\t")

best_traits_clust_QUAL <- best_traits_clust[which(best_traits_clust %in% as.character(metadata_categorical$trait) |
                                                           best_traits_clust %in% as.character(metadata_binary$Trait))]
quali <- sequenced_traits[which(rownames(sequenced_traits) %in% best_species_clust),
                           which(colnames(sequenced_traits) %in% best_traits_clust_categorical)]

var_quali<- MCA(quali)
var_quali_cor <- as.data.frame(var_quali$var$cos2)

fviz_mca_var(var_quali, col.var="cos2", alpha.var="contrib", 
              gradient.cols = c("blue", "green", "red"), repel = TRUE)



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("Biostrings","ggtree"))















test <- sequenced_traits[,which(colnames(sequenced_traits) %in% c("SpeciesBROAD",
                                                                  "BV",
                                                                  "MOB",
                                                                  "AOB"))] %>%
  drop_na()
View(test)

labels <- match(rownames(test), primates_tree$tip.label)
labels_good <- labels[!is.na(labels)]
pruned_tree<-drop.tip(primates_tree,primates_tree$tip.label[-labels_good])
corr <- vcv.phylo(pruned_tree)
rownames(test) <- as.character(test$SpeciesBROAD)
SORTED_DF <- test[match(rownames(corr), rownames(test)), ]
library(caper)
primate_data <- comparative.data(phy = pruned_tree, data = SORTED_DF, names.col = "SpeciesBROAD", vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
pgls(as.numeric(as.character(MOB))~as.numeric(as.character(BV)),
     data= primate_data, lambda="ML")

MOB <- SORTED_DF[, "AOB"]
BV <- SORTED_DF[, "BV"]

# Give them names
names(MOB) <- names(BV) <- rownames(SORTED_DF)

# Calculate PICs
hPic <- pic(MOB, pruned_tree)
aPic <- pic(BV, pruned_tree)
picModel <- lm(hPic ~ aPic - 1)
plot(hPic ~ aPic)
abline(a = 0, b = coef(picModel))




library(cluster)
gap_stat <- clusGap(var_quanti_cor, FUN = kmeans, nstart = 1, K.max = 2, B = 50)
fviz_gap_stat(gap_stat) + theme_minimal() + ggtitle("fviz_gap_stat: Gap Statistic") #valida el cluster único
res.km <- kmeans(var_quanti_cor, centers = 1, nstart = 1)
grp <- as.factor(res.km$cluster)
plot(var_quanti_cor, col = res.km$cluster)
fviz_cluster(res.km, col.var = grp,
             data = var_quanti_cor,
             show.clust.cent = TRUE,
             repel=TRUE)
fviz_nbclust(var_quanti, kmeans, method = "wss", k.max = 2) + theme_minimal() + ggtitle("the Elbow Method")

