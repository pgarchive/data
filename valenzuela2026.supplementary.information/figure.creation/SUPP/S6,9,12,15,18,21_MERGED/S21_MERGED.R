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


rownames(binarized_traits) <- sequenced_traits[,1]


##############################################
#A) ANALYSIS OF CORRELATIONS WITH CLUSTERS###
##############################################


coverage_traits_clustered <- read.csv("../../Figures/Fig1/Traits_clusters_FINAL_JAN23.tsv", sep = "\t")

traits_clust <- coverage_traits_clustered[which(coverage_traits_clustered$Cluster ==10),]

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




#####################################################
#3) MAKE QUALITATIVE CORRELATIONS WITH PHYLIP (THRESHOLD MODEL)##
#####################################################

sessionInfo()
library(mvMORPH)

#CLUSTER7-13
library(Rphylip)


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




#SELECT SPECIES FROM TREE
species_to_keep <- match(rownames(selected_qual_dat), primates_tree$tip.label)
species_to_keep <- species_to_keep[!is.na(species_to_keep)]
pruned.tree<-drop.tip(primates_tree,primates_tree$tip.label[-species_to_keep])
write.tree(pruned.tree, file="pruned_tree.newick")
select_dummy <- select_dummy[rownames(select_dummy) %in% pruned.tree$tip.label,]
X <- select_dummy

#COMPUTED FROM CLUSTER

selected_qual_dat <- read.csv("../S21/qual_dummy_JAN23", sep = "\t", 
                              header=T, row.names = 1)
mat <- read.csv("../S21/all_cors_pasted_JAN23.txt", sep="\t", header = FALSE)
df2 <- mutate_all(mat, function(x) as.numeric(as.character(x)))


rownames(df2) <- colnames(selected_qual_dat)
colnames(df2) <- colnames(selected_qual_dat)


#View(df2)
library(GGally)
diag(df2)<-1
df2[upper.tri(df2)] <- t(df2)[upper.tri(t(df2))]
rownames(df2)
colnames(df2)
p1 <- ggcorr(data = NULL, 
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


pruned.tree <- read.newick("../S21/pruned_tree.newick")
t.xy<-df2*sqrt((Ntip(pruned.tree)-2)/(1-df2^2))
#View(t.xy)
t.xy <- as.matrix(t.xy)
t.xy <- mutate_all(t.xy, function(x) as.numeric(as.character(t.xy)))
P.xy<-2*pt(abs(t.xy),df=Ntip(pruned.tree)-2,lower.tail=F)

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

#View(df_final)
# Ver resultado
#(df_final)

#MELTED <- melt(P.xy)
#MELTED$fdr <- p.adjust(MELTED$value)
#FDR <- MELTED[,c(1,2,4)]

FDR <- df_final[,c(1,2,4)]
#FDR_MATRIX <- acast(FDR, Var1~Var2, value.var='fdr', fill=0)
FDR_MATRIX <- acast(FDR, Var1~Var2, value.var='fdr')
#View(FDR_MATRIX)
DF <-melt(FDR_MATRIX)
#View(DF)

# Ordenar por la media de FDR en vez de la suma, ignorando NA
row_order <- order(rowMeans(FDR_MATRIX, na.rm = TRUE), decreasing = FALSE)

# Aplicar el orden a las filas y columnas
SORTED_MATRIX <- FDR_MATRIX[row_order, row_order]
SORTED_MATRIX_COLUMNS <- SORTED_MATRIX[,rownames(SORTED_MATRIX)]

trait_domains <- read.csv("Traits_domains_SPLIT_FEB22.txt", sep = "\t")

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




DF_domains <- dplyr::left_join(DF,trait_domains, by=c("TRAIT" = "Measurement"))

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
SORTED_MERGED3 <- left_join(SORTED_MERGED, SORTED_DMC, by=c("Var1" = "x"))


#View(SORTED_MERGED3)
p1 <- ggplot(data = SORTED_MERGED3) +
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
  
  scale_fill_manual(values = traits_palette[c(1,2,5)])

############################################
#MERGING STYEPS /COMBINE BOTH DIAGONAL INFO#
############################################
# 1. Paso: convertir la matriz de correlaciones en formato largo
# 1. Preparar df de correlaciones en formato largo

cor_df <- as.data.frame(as.table(as.matrix(df2))) %>%
  rename(Var1 = Var1, Var2 = Var2, cor = Freq) %>%
  filter(as.character(Var1) < as.character(Var2))

# 2. Preparar df de FDR desde SORTED_MERGED2
fdr_df <- SORTED_MERGED3 %>%
  rename(Var1 = 1, Var2 = 2, fdr = value)

# Función para normalizar pares Var1 y Var2 (orden alfabético)
normalize_pairs <- function(df, var1_col = "Var1", var2_col = "Var2") {
  v1 <- as.character(df[[var1_col]])
  v2 <- as.character(df[[var2_col]])
  
  df[[var1_col]] <- ifelse(v1 < v2, v1, v2)
  df[[var2_col]] <- ifelse(v1 < v2, v2, v1)
  
  return(df)
}

# 3. Normalizar pares en ambos datasets
cor_df <- normalize_pairs(cor_df, "Var1", "Var2")
fdr_df <- normalize_pairs(fdr_df, "Var1", "Var2")

# 4. Duplicar filas de FDR para tener espejo (triángulo inferior)
fdr_df_mirror <- fdr_df %>%
  mutate(orig_Var1 = Var1, orig_Var2 = Var2) %>%
  transmute(
    Var1 = orig_Var2,
    Var2 = orig_Var1,
    fdr = fdr
  )


# 5. Combinar FDR normal + espejo
fdr_full <- bind_rows(fdr_df, fdr_df_mirror) %>%
  distinct(Var1, Var2, .keep_all = TRUE)

#cor mirror
cor_df_mirror <- cor_df %>%
  rename(orig_Var1 = Var1, orig_Var2 = Var2) %>%
  transmute(
    Var1 = orig_Var2,
    Var2 = orig_Var1,
    cor = cor
  )

cor_df_full <- bind_rows(cor_df, cor_df_mirror) %>%
  distinct(Var1, Var2, .keep_all = TRUE)

# 6. Matriz combinada de correlación (solo triángulo superior) y FDR (completo)
combined_full <- full_join(
  cor_df_full %>% select(Var1, Var2, cor),
  fdr_full %>% select(Var1, Var2, fdr),
  by = c("Var1", "Var2")
)


# 7. Obtener mínimo FDR para cada trait en Var1
min_fdr_var1 <- combined_full %>%
  filter(!is.na(fdr)) %>%
  group_by(Var1) %>%
  summarise(min_fdr = min(fdr, na.rm = TRUE))

# 8. Obtener mínimo FDR para cada trait en Var2
min_fdr_var2 <- combined_full %>%
  filter(!is.na(fdr)) %>%
  group_by(Var2) %>%
  summarise(min_fdr = min(fdr, na.rm = TRUE))

# 9. Renombrar para juntar
min_fdr_var2_renamed <- min_fdr_var2 %>%
  rename(Var1 = Var2)

# 10. Combinar mínimos y calcular mínimo combinado
min_fdr_traits <- full_join(min_fdr_var1, min_fdr_var2_renamed, by = "Var1") %>%
  mutate(
    trait = Var1,
    min_fdr_combined = pmin(min_fdr.x, min_fdr.y, na.rm = TRUE)
  ) %>%
  select(trait, min_fdr_combined) %>%
  distinct()

# 11. Ordenar traits según min_fdr_combined ascendente
ordered_traits <- min_fdr_traits %>%
  arrange(min_fdr_combined) %>%
  pull(trait)

# 12. Reordenar factores en combined_full con el orden calculado
combined_full <- combined_full %>%
  mutate(
    Var1 = factor(Var1, levels = ordered_traits),
    Var2 = factor(Var2, levels = ordered_traits),
    x = as.numeric(Var1),
    y = as.numeric(Var2)
  ) %>%
  filter(as.character(Var1) != as.character(Var2))  # Eliminamos diagonal

combined_full %>% filter(as.character(Var1) == as.character(Var2))

# 13. Gráfico final
cor_upper <- filter(combined_full, !is.na(cor), y > x)
fdr_lower <- filter(combined_full, !is.na(fdr), y < x)

#COLOR
# Crear dataframe para las barras laterales de color por dominio

trait_colors_df1 <- SORTED_MERGED3 %>%
  select(trait = Var2, colour = colour.x) %>%
  distinct() %>%                    # Para eliminar duplicados, si los hubiera
  filter(trait %in% ordered_traits) %>%  # Solo los traits que quieres
  mutate(pos = match(trait, ordered_traits)) %>%
  arrange(pos)
trait_colors_df2 <- SORTED_MERGED3 %>%
  select(trait = Var1, colour = colour.y) %>%
  distinct() %>%                    # Para eliminar duplicados, si los hubiera
  filter(trait %in% ordered_traits) %>%  # Solo los traits que quieres
  mutate(pos = match(trait, ordered_traits)) %>%
  arrange(pos)
trait_colors_df <- rbind(trait_colors_df1, trait_colors_df2)

library(ggnewscale)
traits_palette <- c(
  "Behavior" = "#FF0000",
  "Ecology" = "#00A08A",
  "Life-history" = "#F2AD00",
  "Morphology" = "navyblue",
  "Physiology" = "#daa18f"
)

# Plot

library(ggnewscale)
p1 <- ggplot() +
  # ─── Primary domain (Marginal bars) ─────────────────────────
  geom_tile(
    data = trait_colors_df,
    aes(x = -0.45, y = pos, fill = colour),
    width = 0.5
  ) +
  geom_tile(
    data = trait_colors_df,
    aes(x = pos, y = -0.45, fill = colour),
    height = 0.5
  ) +
  scale_fill_manual(
    values = traits_palette[c("Morphology", "Life-history", "Behavior", "Ecology")],
    name = "Primary domain"
  ) 


p2 <- p1 + 
  # ─── Reset fill scale ───────────────────────────────────────
  new_scale_fill() +
  
  # ─── Correlation (Upper triangle) ───────────────────────────
  geom_tile(data = filter(cor_upper, x != y), aes(x = x, y = y, fill = cor)) +
  scale_fill_gradient2(
    low = "darkred", mid = "white", high = "darkgreen",
    midpoint = 0, name = "Correlation (r)"
  ) 

clean_guides <- function(plot_obj) {
  plot_obj$scales$scales <- lapply(plot_obj$scales$scales, function(s) {
    if (!is.null(s$guide) && !is.atomic(s$guide) && !is.list(s$guide)) {
      s$guide <- "legend"
    }
    s
  })
  plot_obj
}

p2 <- clean_guides(p2)

p3 <- p2 +
  
  new_scale_fill() +
  
  # ─── FDR (Lower triangle) ───────────────────────────────────
  geom_tile(data = filter(fdr_lower, x != y), aes(x = x, y = y, fill = fdr)) +
  scale_fill_gradient(
    low = "darkblue", high = "white",
    limits = c(0, 1), na.value = "grey90",
    name = "FDR p-value"
  ) +
  
  # ─── Layout and theme ───────────────────────────────────────
  coord_fixed() +
  scale_x_continuous(
    breaks = 1:length(ordered_traits), labels = ordered_traits, expand = c(0, 0)
  ) +
  scale_y_continuous(
    breaks = 1:length(ordered_traits), labels = ordered_traits, expand = c(0, 0)
  ) +
  theme_minimal() +
  theme(
    legend.key.size = unit(0.4, "cm"),
    axis.text.x = element_text(angle = 90, hjust = 1),
    axis.text.y = element_text(size = 8),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank()
  ) 
  