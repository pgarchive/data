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

traits_clust <- coverage_traits_clustered[which(coverage_traits_clustered$Cluster ==6),]

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

#install.packages("~/Descargas/devtools_2.4.3.tar.gz", repos = NULL, type="source")
#install_github("shinepreventer/ggcor-hou")


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

#matrix(FDR_MATRIX[order(row(FDR_MATRIX), FDR_MATRIX)], ncol=ncol(FDR_MATRIX), byrow=TRUE)

SORTED_MATRIX <- FDR_MATRIX[order(apply(FDR_MATRIX, 1, function(row) sum(row)), decreasing = F), ]
SORTED_MATRIX_COLUMNS <- SORTED_MATRIX[,rownames(SORTED_MATRIX)]



DF <-melt(FDR_MATRIX)

# Ordenar por la media de FDR en vez de la suma, ignorando NA
row_order <- order(rowMeans(FDR_MATRIX, na.rm = TRUE), decreasing = FALSE)

# Aplicar el orden a las filas y columnas
SORTED_MATRIX <- FDR_MATRIX[row_order, row_order]
SORTED_MATRIX_COLUMNS <- SORTED_MATRIX[,rownames(SORTED_MATRIX)]

#SORTED_MATRIX <- FDR_MATRIX[order(apply(FDR_MATRIX, 1, function(row) sum(row)), decreasing = F), ]
#SORTED_MATRIX_COLUMNS <- SORTED_MATRIX[,rownames(SORTED_MATRIX)]

trait_domains <- read.csv("Traits_domains_SPLIT_FEB22.txt", sep = "\t")


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

  scale_fill_manual(values = traits_palette[c(2,3,4)]) +
  theme(legend.position = "right")
  

saveRDS(p1, "FigureS9")



