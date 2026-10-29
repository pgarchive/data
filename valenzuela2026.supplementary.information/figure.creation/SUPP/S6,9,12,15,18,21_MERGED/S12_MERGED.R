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


#B) SELECT DESIRED SPECIES (sequenced_ones)


primates_palette <- readRDS("../../Figures/Fig1/primate_families_palette.rds")
seq_species <- read.csv("../../Figures/Fig1/names_FINAL.txt",
                        sep = "\t", header = TRUE)
colnames(seq_species)<- c("Tree_Species", "Genome_species", "Trait_species")
covered_species_all <- as.character(seq_species$Trait_species[!is.na(seq_species$Tree_Species)])
covered_species <- as.character(covered_species_all[!is.na(covered_species_all)])
sequenced_traits <- all_traits[which(all_traits$SpeciesBROAD %in% covered_species),]
#levels(sequenced_traits$SpeciesBROAD) <- c(levels(sequenced_traits$SpeciesBROAD), "Tarsius_bancanus", "Leontocebus_fuscicolis",
#                                           "Plecturocebus_hoffmansi", "Lagothrix_lagothricha") 
#sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Cephalopachus_bancanus"] <- "Tarsius_bancanus"
#sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Leontocebus_fuscicollis"] <- "Leontocebus_fuscicolis"
#sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Plecturocebus_hoffmannsi"] <- "Plecturocebus_hoffmansi"
#sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Lagothrix_lagotricha"] <- "Lagothrix_lagothricha"

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

redundant <- read.csv("../S6/Redundant_measurements.tsv", header=FALSE, sep = "\t")

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

traits_clust <- coverage_traits_clustered[which(coverage_traits_clustered$Cluster ==5),]

data_clust <- binarized_traits[,which(colnames(binarized_traits) %in% traits_clust$Trait)]

num_data_clust <- sapply(data_clust, as.numeric)
rownames(num_data_clust) <- rownames(data_clust)
##Subset as usual
best_species_clust <- rownames(num_data_clust[which(rowSums(num_data_clust) >= 0.90*ncol(num_data_clust)),])
best_species_clust
all_traits_clust <- binarized_traits[which(rownames(binarized_traits) %in% best_species_clust),]

#filter traits for correlation study
#all_traits_clust <- all_traits_clust[,-ncol(all_traits_clust)] %>%
#  mutate_all(as.numeric)


#Cluster6 SEXUAL DIMORPHISM[8] - SPECIES WITH >90% coverafe of traits in cluster | NO MISS-IMPUTE



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

#selected_subset <- selected_subset %>% 
# mutate(SocialSystem = ifelse(as.character(SocialSystem) == "Polygynandry", "Polygynandry_sec", as.character(SocialSystem))) %>%
#mutate(SocialSystem = ifelse(as.character(SocialSystem) == "Pair", "Pair_sec", as.character(SocialSystem))) %>%
#  mutate(SocialSystem = ifelse(as.character(SocialSystem) == "Solitary", "Solitary_sec", as.character(SocialSystem))) %>%
#mutate(SocialSystem = ifelse(as.character(SocialSystem) == "Multimale", "Multimale_sec", as.character(SocialSystem))) 


library(missMDA)
no_missing <- selected_subset %>% drop_na()
#imputed <- imputeFAMD(selected_subset)

res.famd <- FAMD(no_missing) #cluster6 | cluster7 | cluster 8
#res.famd <- FAMD(selected_subset, tab.disj = imputed$tab.disj)
#res.famd_raw <- FAMD(selected_subset)

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
#REMOVE missinf CORRs from the dataset
correl <- correl[rowSums(is.na(correl)) != ncol(correl)-1,]
correl <- correl[,colSums(is.na(correl)) != nrow(correl)]
#correl <- correl[rowSums(correl==0) != ncol(correl)-1,]
#correl <- correl[,colSums(correl==0) != nrow(correl)]
#nrow(correl)
#ncol(correl)


#adjust corr
#library(remotes)
#remotes::install_github("bcdudek/bcdstats")
#correl_struct$C

#correl <- correl[which(!rownames(correl)=="BodyFlesh"),  ]
#correl <- correl[,which(!colnames(correl)=="BodyFlesh") ]
#adjust.corr(
# correl,
#  type = c("pearson"),
# use = c("pairwise.complete.obs"),
#adjust = "holm"
#)
t.xy<-correl*sqrt((Ntip(clust.tree)-2)/(1-correl^2))
#View(t.xy)
P.xy<-2*pt(abs(t.xy),df=Ntip(clust.tree)-2,lower.tail=F)


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
View(df_final)
FDR <- df_final[,c(1,2,4)]
#FDR <- MELTED[,c(1,2,4)]
FDR_MATRIX <- acast(FDR, Var1~Var2, value.var='fdr', fill=0)

"""
library(TestCor)
ApplyFdrCor(
correl,
alpha = 0.05,
stat_test = "student",
method = "BH",
Nboot = 1000,
vect = FALSE,
arr.ind = FALSE
)
"""

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

DF <-melt(FDR_MATRIX)
library(ComplexHeatmap)
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
#colnames(DF) <- c("Var1", "Var2", "r")
#inverse <- 1- FDR_MATRIX
#algo raro pasa, no coge bien los p-values y tengo que cambiar los labels
#quickcor(SORTED_MATRIX_COLUMNS, type = "upper", show.diag=FALSE) +
# ggcor::anno_link(data = DF_summ_SUMM, aes(colour = V1, label = V1)) +
#scale_colour_manual(values = traits_palette[c(1,2,3,4,6)]) +
#scale_fill_gradient(high = "#132B43", low = "#56B1F7",
#                    labels=c('0', '0.5', '0.75'),
#                   breaks=c(1, 0.5, 0)) +
#geom_color() +
#remove_axis("x") +
#labs(fill='FDR p-value', colour='Biological domains') +
#guides(colour=FALSE)

# Ordenar por la media de FDR en vez de la suma, ignorando NA
row_order <- order(rowMeans(FDR_MATRIX, na.rm = TRUE), decreasing = FALSE)

# Aplicar el orden a las filas y columnas
SORTED_MATRIX <- FDR_MATRIX[row_order, row_order]
SORTED_MATRIX_COLUMNS <- SORTED_MATRIX[,rownames(SORTED_MATRIX)]

#SORTED_MATRIX <- FDR_MATRIX[order(apply(FDR_MATRIX, 1, function(row) sum(row)), decreasing = F), ]
#SORTED_MATRIX_COLUMNS <- SORTED_MATRIX[,rownames(SORTED_MATRIX)]


trait_domains <- read.csv("../S7/Traits_domains_SPLIT_FEB22.txt", sep = "\t")


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
p2 <- ggplot(data = SORTED_MERGED2) +
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
  
  scale_fill_manual(values = traits_palette[c(1,2,3,4,5)]) +
  theme(legend.position = "right")


############################################
#MERGING STYEPS /COMBINE BOTH DIAGONAL INFO#
############################################
# 1. Paso: convertir la matriz de correlaciones en formato largo
# 1. Preparar df de correlaciones en formato largo
cor_df <- as.data.frame(as.table(correl)) %>%
  rename(Var1 = Var1, Var2 = Var2, cor = Freq) %>%
  filter(as.character(Var1) < as.character(Var2))

# 2. Preparar df de FDR desde SORTED_MERGED2
fdr_df <- SORTED_MERGED2 %>%
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

trait_colors_df1 <- SORTED_MERGED2 %>%
  select(trait = Var2, colour = colour.x) %>%
  distinct() %>%                    # Para eliminar duplicados, si los hubiera
  filter(trait %in% ordered_traits) %>%  # Solo los traits que quieres
  mutate(pos = match(trait, ordered_traits)) %>%
  arrange(pos)
trait_colors_df2 <- SORTED_MERGED2 %>%
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
    values = traits_palette,
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
    axis.text.x = element_text(angle = 90, hjust = 1, size = 7),
    axis.text.y = element_text(size = 7),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank()
  ) 
  