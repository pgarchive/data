########################
##Functional Enrich#####
########################

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
library(WebGestaltR)
library(randomcoloR)

View(listGeneSet(organism = "hsapiens"))
listIdType()

trait_domains <- read.csv("../../Figures/Fig2/Traits_domains_SPLIT_FEB22_ALL_UNDERSCORES.txt", header=T, sep = "\t")

#READ ENRICHMENTS fdr<0.05

enrichments <- read.csv("../../Figures/Fig2/ALL_ENRICHMENTS_SEP23_FINAL.txt", header=F, sep = "\t")
colnames(enrichments) <- c("geneSet", "description", "link", "size", 
                           "overlap",	"expect",	"enrichmentRatio",	"pValue",
                           "FDR",	"overlapId", "userId", "DB", "Method", "Trait")

#FIGURE CAAS
traits_with_no_contrast <- c("Basophils_Abs",
                             "Basophils_",
                             "Blood_Urea_Nitrogen",
                             "Bun_Creat_Ratio",
                             "Cholesterol",
                             "Diencephalon",
                             "GGT",
                             "Globulin",
                             "Hippocampus",
                             "Hypothalamus",
                             "Maximum_lifespan__y",
                             "M_Max_M_NeonateWeight_g",
                             "MOB",
                             "Piriform_Lobe",
                             "Schizocortex",
                             "Ser_Glut_Oxaloacetic_Transamin",
                             "Striatum",
                             "Striatum__incl__NAcc",
                             "Thalamus",
                             "Uric_Acid",
                             "XII")


new_names <- c()
for (i in 1:length(enrichments$Trait)){
  Name <- enrichments$Trait[i]
  #print(Name)
  list_filtered <- strsplit(as.character(Name), split="_", fixed=TRUE)
  #print(list_filtered)
  set_of_characters <- list_filtered[[1]][2:length(list_filtered[[1]])]
  new_name <- paste(set_of_characters, collapse="_")
  new_names <- c(new_names, new_name)
}

enrichments$trait <- new_names


colnames(trait_domains) <- c("trait", "Domain", "Secondary_Domain")

enrichments$TRAIT <- enrichments$trait
levels(enrichments$TRAIT) <- c(levels(enrichments$TRAIT) ,
                               c("Colour_vision", "Social_organization",
                                 "Mating_system", "DielActivity",
                                 "TrophicGuild", "Locomotion",
                                 "Habitat", "SocialMatingSystem", "Monopolization", "NatalDispersal"))

enrichments$TRAIT[startsWith(as.character(enrichments$trait), "Colour_vision_")] <- "Colour_vision"
enrichments$TRAIT[startsWith(as.character(enrichments$trait), "Social_organization_")] <- "Social_organization"
enrichments$TRAIT[startsWith(as.character(enrichments$trait), "Mating_system_")] <- "Mating_system"
enrichments$TRAIT[startsWith(as.character(enrichments$trait), "DielActivity_")] <- "DielActivity"
enrichments$TRAIT[startsWith(as.character(enrichments$trait), "TrophicGuild_")] <- "TrophicGuild"
enrichments$TRAIT[startsWith(as.character(enrichments$trait), "Locomotion_")] <- "Locomotion"
enrichments$TRAIT[startsWith(as.character(enrichments$trait), "Habitat_")] <- "Habitat"
enrichments$TRAIT[startsWith(as.character(enrichments$trait), "SocialMatingSystem_")] <- "SocialMatingSystem"
enrichments$TRAIT[startsWith(as.character(enrichments$trait), "Monopolization_")] <- "Monopolization"
enrichments$TRAIT[startsWith(as.character(enrichments$trait), "NatalDispersal_")] <- "NatalDispersal"

enrichments_domains <- dplyr::left_join(enrichments,
                                        trait_domains, by=c("TRAIT" = "trait"))

#View(enrichments_domains)

CAAS_enrichments <- enrichments_domains[which(!enrichments_domains$trait %in% traits_with_no_contrast &
                                                enrichments_domains$Method == "CAAS"),]
RERS_enrichments <- enrichments_domains[which(enrichments_domains$Method == "RERS"),]

enrichments_domains <- rbind(CAAS_enrichments, RERS_enrichments)

counts_RERS <- RERS_DOMAIN %>%
  group_by(Domain) %>%
  summarize(distinct_trait = n_distinct(trait))

counts_CAAS <- CAAS_DOMAIN %>%
  group_by(Domain) %>%
  summarize(distinct_trait = n_distinct(trait))

joined <- bind_rows(counts_CAAS, counts_RERS, .id = 'Method')
joined$Method <- ifelse(joined$Method == 1, "CAAS", "RERS")


summary_domains <- enrichments_domains %>%
  group_by(DB, Secondary_Domain, Method) %>%
  summarize(distinct_trait = n_distinct(trait))

ggplot(data=enrichments_domains[which(enrichments_domains$DB == "Disgenet_disease"),]) +
  geom_histogram(aes(x=Trait, fill=Method), stat="count") +
  coord_flip() +
  facet_wrap(~Method, scales="free") +
  theme_minimal()

ggplot(data=enrichments[which(enrichments$DB == "pathway_KEGG"),]) +
  geom_histogram(aes(x=Trait, fill=Method), stat="count") +
  coord_flip() +
  scale_fill_manual(values = c("darkgoldenrod3", "steelblue4")) +
  facet_wrap(~Method, scales="free") +
  theme_minimal() +
  theme(strip.text = element_blank(),
        legend.position = "none") +
  ylab("Nº of enrichments") +
  xlab("Trait measurement")

library(ggmosaic)
test <- summary_domains[which(summary_domains$DB == "MF" |
                                summary_domains$DB == "CC" |
                                summary_domains$DB == "BP" ),]
test$Secondary_Domain <- factor(test$Secondary_Domain, 
                                           levels = rev(unique(test$Secondary_Domain)))
test$Secondary_Domain <- droplevels(test$Secondary_Domain)

ggplot(data=test) +
  geom_mosaic(aes(x = product(Secondary_Domain), fill = Method,
                  weight = as.numeric(as.character(distinct_trait))),
              na.rm = TRUE, show.legend = FALSE) +
  scale_fill_manual(values = c("darkgoldenrod3", "steelblue4")) +
  coord_flip() +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=10),
        plot.title = element_text(hjust = 0.5))


#SUPP PLOTS 50-53
t1 <- ggplot(data=enrichments_domains[which(enrichments_domains$DB == "pathway_Panther"),]) +
  geom_histogram(aes(x=trait, fill=Method), stat="count") +
  coord_flip() +
  scale_fill_manual(values = c("darkgoldenrod3", "steelblue4")) +
  # facet_wrap(~Method, scales="free") +
  theme_minimal() +
  theme(strip.text = element_blank(),
        legend.position = "none") +
  ylab("Nº of enrichments") +
  xlab("Trait measurement")

saveRDS(t1, "FigureS50")

