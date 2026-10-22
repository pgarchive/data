###########################
#VISUALIZE MANHATTAN PLOTS#
###########################

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
library(harmonicmeanp)
library(combinationpvalues)


#FROM HERE
joined <- readRDS("joined_FEB24_QUAL")
#write.table(joined, "joined_FEB24_QUAL.txt", sep="\t", quote=F, row.names=F,
 #           col.names=F)

genes_coordinates <- read.csv("../../Figures/Fig2/genes_coordinates_sorted.txt", header=T, sep = "\t")

colnames(genes_coordinates) <- c("gene", "Chrom", "bpIN", "bpEND")
GENES_DF <- dplyr::left_join(joined, 
                             genes_coordinates, by=c("gene"))  %>%
  drop_na()

#colnames(trait_domains)

#PLOT CHRX

Gen_Phen_dataset <- GENES_DF[which(!is.na(GENES_DF$Chrom) &
                                     GENES_DF$Chrom %in% c(1:22,"X")),]

levels(Gen_Phen_dataset$Chrom) <- c(levels(Gen_Phen_dataset$Chrom), 23)
Gen_Phen_dataset$Chrom[Gen_Phen_dataset$Chrom == "X"] <- 23
Gen_Phen_dataset$Chrom <- as.numeric(as.character(Gen_Phen_dataset$Chrom))

#III) ADJUST P-VALUES FOR CAAS

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

adjusted_pvalues_harmonic_CAAS_DATAFRAME <- Gen_Phen_dataset[which(Gen_Phen_dataset$Method == "CAAS" &
                                                                     !Gen_Phen_dataset$trait %in% traits_with_no_contrast),] %>%
  group_by(Method, trait) %>%
  mutate(adjusted_harm_fdr = p.adjust(mean_permpval, method = "fdr",
                                      n=length(mean_permpval[!is.na(mean_permpval)]))) %>%
  mutate(adjusted_harm_bonf = p.adjust(mean_permpval, method = "bonf",
                                       n=length(mean_permpval[!is.na(mean_permpval)]))) %>%
  ungroup(Method, trait) %>%
  group_by(Method, gene) %>%
  slice(which.min(mean_permpval))

adjusted_pvalues_harmonic_CAAS <- Gen_Phen_dataset[which(Gen_Phen_dataset$Method == "CAAS" &
                                                           !Gen_Phen_dataset$trait %in% traits_with_no_contrast),] %>%
  group_by(Method, trait) %>%
  mutate(adjusted_harm_fdr = p.adjust(mean_permpval, method = "fdr")) %>%
  mutate(adjusted_harm_bonf = p.adjust(mean_permpval, method = "bonf")) %>%
  ungroup(Method, trait) %>%
  group_by(Method, gene) %>%
  slice(which.min(mean_permpval)) %>%
  ungroup(Method, gene) %>%
  summarize(index_threeshold_fdr = which.min(abs(sort(adjusted_harm_fdr)-0.05)),
            index_threeshold_bonf = which.min(abs(sort(adjusted_harm_bonf)-0.05)))

#III) ADJUST P-VALUES FOR RERS
adjusted_pvalues_harmonic_RERS_DATAFRAME <- Gen_Phen_dataset[which(Gen_Phen_dataset$Method == "RERS"),] %>%
  group_by(Method, trait) %>%
  mutate(adjusted_harm_fdr = p.adjust(mean_permpval, method = "fdr",
                                      n=length(mean_permpval[!is.na(mean_permpval)]))) %>%
  mutate(adjusted_harm_bonf = p.adjust(mean_permpval, method = "bonf",
                                       n=length(mean_permpval[!is.na(mean_permpval)]))) %>%
  ungroup(Method, trait) %>%
  group_by(Method, gene) %>%
  slice(which.min(mean_permpval))

adjusted_pvalues_harmonic_RERS <- Gen_Phen_dataset[which(Gen_Phen_dataset$Method == "RERS"),] %>%
  group_by(Method, trait) %>%
  mutate(adjusted_harm_fdr = p.adjust(mean_permpval, method = "fdr")) %>%
  mutate(adjusted_harm_bonf = p.adjust(mean_permpval, method = "bonf")) %>%
  ungroup(Method, trait) %>%
  group_by(Method, gene) %>%
  slice(which.min(mean_permpval)) %>%
  ungroup(Method, gene) %>%
  summarize(index_threeshold_fdr = which.min(abs(sort(adjusted_harm_fdr)-0.05)),
            index_threeshold_bonf = which.min(abs(sort(adjusted_harm_bonf)-0.05)))


#IV) Compute thresholds for manhattan figure
#SUBSET_RERS <- Gen_Phen_dataset[which(Gen_Phen_dataset$Method == "RERS"),]
#SUBSET_CAAS <- Gen_Phen_dataset[which(Gen_Phen_dataset$Method == "CAAS"),]  %>% 
  #!Gen_Phen_dataset$trait %in% traits_with_no_contrast),] %>%
adjusted_pvalues <- bind_rows(adjusted_pvalues_harmonic_CAAS, 
                                adjusted_pvalues_harmonic_RERS, 
                                .id = 'Method')
adjusted_pvalues$Method <- ifelse(adjusted_pvalues$Method == 1, "CAAS", "RERS")
adjusted_pvalues_harmonic_RERS_DATAFRAME <- adjusted_pvalues_harmonic_RERS_DATAFRAME %>%
  arrange(as.numeric(as.character(adjusted_harm_fdr)))
adjusted_pvalues_harmonic_CAAS_DATAFRAME <- adjusted_pvalues_harmonic_CAAS_DATAFRAME %>%
  arrange(as.numeric(as.character(adjusted_harm_fdr)))
threshold_pvalue_fdr_RER <- adjusted_pvalues_harmonic_RERS_DATAFRAME[as.numeric(as.character(adjusted_pvalues_harmonic_RERS$index_threeshold_fdr)), ]
threshold_pvalue_bonf_RER <- adjusted_pvalues_harmonic_RERS_DATAFRAME[as.numeric(as.character(adjusted_pvalues_harmonic_RERS$index_threeshold_bonf)), ]
threshold_pvalue_fdr_CAAS <- adjusted_pvalues_harmonic_CAAS_DATAFRAME[as.numeric(as.character(adjusted_pvalues_harmonic_CAAS$index_threeshold_fdr)), ]
threshold_pvalue_bonf_CAAS <- adjusted_pvalues_harmonic_CAAS_DATAFRAME[as.numeric(as.character(adjusted_pvalues_harmonic_CAAS$index_threeshold_bonf)), ]
logthreshold_fdr_RER <- -log10(max(adjusted_pvalues_harmonic_RERS_DATAFRAME$mean_permpval[adjusted_pvalues_harmonic_RERS_DATAFRAME$adjusted_harm_fdr<0.05]))
logthreshold_bonf_RER <- -log10(max(adjusted_pvalues_harmonic_RERS_DATAFRAME$mean_permpval[adjusted_pvalues_harmonic_RERS_DATAFRAME$adjusted_harm_bonf<0.05]))
logthreshold_fdr_CAAS<- -log10(max(adjusted_pvalues_harmonic_CAAS_DATAFRAME$mean_permpval[adjusted_pvalues_harmonic_CAAS_DATAFRAME$adjusted_harm_fdr<0.05]))
logthreshold_bonf_CAAS <- -log10(max(adjusted_pvalues_harmonic_CAAS_DATAFRAME$mean_permpval[adjusted_pvalues_harmonic_CAAS_DATAFRAME$adjusted_harm_bonf<0.05]))

all_log <- as.data.frame(rbind(c(logthreshold_bonf_CAAS, logthreshold_fdr_CAAS),
                               c(logthreshold_bonf_RER, logthreshold_fdr_RER)))
all_log$Method <- c("CAAS", "RERS")
colnames(all_log) <-  c("bonf", "fdr", "Method")

#V) REFORMAT COORDINATE+P-VALUE INFORMATION FOR MANHATTAN PLOT 

#View(adjusted_pvalues_harmonic_RERS_DATAFRAME)
genes_df_rers <- adjusted_pvalues_harmonic_RERS_DATAFRAME %>% 
  group_by(Chrom) %>% 
  summarise(chr_len=max(bpEND)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(Chrom, tot) %>%
  
  # Add this info to the initial dataset
  left_join(adjusted_pvalues_harmonic_RERS_DATAFRAME, ., by=c("Chrom"="Chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chrom, bpEND) %>%
  mutate( BPcum=bpEND+tot) %>%
  drop_na() 
#  mutate( is_highlight=ifelse(-log10(mean_permpval)>logthreshold, "yes", "no")) %>%
#  mutate( is_annotate=ifelse(-log10(mean_permpval)>logthreshold, "yes", "no")) 

# unique(genes_df$Chrom)
axisdf <- genes_df_rers %>% group_by(Chrom) %>% 
  summarise(center=( max(BPcum, na.rm = TRUE) + min(BPcum, na.rm = TRUE) ) / 2 )

#CAAS

#unique(genes_df_caas$trait)
genes_df_caas <- adjusted_pvalues_harmonic_CAAS_DATAFRAME %>% 
  group_by(Chrom) %>% 
  summarise(chr_len=max(bpEND)) %>%
  # Calculate cumulative position of each chromosome
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  dplyr::select(Chrom, tot) %>%
  
  # Add this info to the initial dataset
  left_join(adjusted_pvalues_harmonic_CAAS_DATAFRAME, ., by=c("Chrom"="Chrom")) %>%
  
  # Add a cumulative position of each SNP
  arrange(Chrom, bpEND) %>%
  mutate( BPcum=bpEND+tot) %>%
  drop_na() 

#  mutate( is_highlight=ifelse(-log10(mean_permpval)>logthreshold, "yes", "no")) %>%
#  mutate( is_annotate=ifelse(-log10(mean_permpval)>logthreshold, "yes", "no")) 

# unique(genes_df$Chrom)
axisdf_caas <- genes_df_caas %>% group_by(Chrom) %>% 
  summarise(center=( max(BPcum, na.rm = TRUE) + min(BPcum, na.rm = TRUE) ) / 2 )


#VI) RETRIEVE LIST/SUBSET OF GENE-TRAIT ASSOCIATIONS above the THREESHOLD
best_CAAS <- genes_df_caas[which(genes_df_caas$mean_permpval<1e-6),]
#best_CAAS[which(best_CAAS$trait == "MOB"),]
best_RERS <- genes_df_rers[which(genes_df_rers$mean_permpval<1e-8),]
#View(best_RERS)
#View(best_RERS)
best <- bind_rows(best_CAAS, best_RERS, .id = 'Method')
best$Method <- ifelse(best$Method == 1, "CAAS", "RERS")

#VII) MAKE MAKE MANHATTAN PLOT FOR RERS
genes_df_fdr_caas <- genes_df_caas %>%
  group_by(Method, trait) %>%
  mutate(fdr = p.adjust(mean_permpval, method = "fdr",
                        n=length(mean_permpval[!is.na(mean_permpval)]))) %>%
  mutate(color = case_when(as.numeric(as.character(adjusted_harm_fdr)) < 0.05 ~ "FDR-trait<0.05",
                           as.numeric(as.character(mean_permpval)) <0.05 ~ "p-val<0.05",
                           TRUE ~ "p-val>0.05")) %>%
  group_by(Method, gene) %>%
  slice(which.min(mean_permpval))

genes_df_fdr_rers <- genes_df_rers %>%
  group_by(Method, trait) %>%
  mutate(fdr = p.adjust(mean_permpval, method = "fdr",
                        n=length(mean_permpval[!is.na(mean_permpval)]))) %>%
  mutate(color = case_when(as.numeric(as.character(adjusted_harm_fdr)) < 0.05 ~ "FDR-trait<0.05",
                           as.numeric(as.character(mean_permpval)) <0.05 ~ "p-val<0.05",
                           TRUE ~ "p-val>0.05")) %>%
  group_by(Method, gene) %>%
  slice(which.min(mean_permpval))

#View(genes_df_fdr_rers)

colfunc_CAAS <- colorRampPalette(c("darkgoldenrod3", "white"))
CAAS_palette <- colfunc_CAAS(10)[c(1,5,7)]
colfunc_RERs <- colorRampPalette(c("steelblue4", "white"))
RERs_palette <- colfunc_RERs(10)[c(1,5,7)]

genes_df_fdr_caas <- genes_df_fdr_caas %>%
  arrange(color)# %>%
# filter(color!="p-val>0.05")
genes_df_fdr_rers <- genes_df_fdr_rers %>%
  arrange(color) %>%
  filter(color!="p-val>0.05")


#VIII) MAKE MANHATTAN PLOT FOR CAAS

#View(genes_df_fdr_caas)
p_caas <- ggplot(genes_df_fdr_caas,
                 aes(x=BPcum, 
                     y=-log10(mean_permpval))) +
  
  # Show all points
  geom_point( aes(color=color, size=3)) +
  scale_x_continuous( label = c(axisdf_caas$Chrom[1:22],"X"), breaks= axisdf_caas$center ) +
  scale_color_manual(values = CAAS_palette ) +     # remove space between plot area and x axis
  
  geom_label_repel( data=best_CAAS, aes(label=as.character(gene), 
                                        y=-log10(mean_permpval)),
                    size=3,
                    max.overlaps = Inf,
                    fontface = 'italic') +
  geom_hline(data = all_log[which(all_log$Method == "CAAS"),], 
             aes(yintercept=bonf), lty=2, col="red") +
  geom_hline(data = all_log[which(all_log$Method == "CAAS"),], 
             aes(yintercept=fdr), lty=2, col="darkgreen") +
  # Custom the theme:
  theme_classic() +
  scale_size(guide = "none") +
  theme( 
    legend.position = "top",
    legend.title = element_blank(),
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    axis.text.x = element_text(size = 7,
                               angle = 90),
    axis.text.y = element_text(size = 10),
    strip.background = element_blank(),
    strip.text = element_blank()
  ) +
  ylab(expression("-Log10( AA-trait"~italic("p")["Phylolm"]~")")) +
  guides(color = guide_legend(override.aes = list(size=4))) +
  xlab("Chromosome") +
  facet_wrap(~Method, scales="free")

saveRDS(p_caas, "CAAS_MANHATTAN_FDR-TRAIT_FEB24_nompval_contrasts")



