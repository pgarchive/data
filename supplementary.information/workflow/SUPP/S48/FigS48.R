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

library(ggmosaic)
test <- summary_domains[which(summary_domains$DB == "MF" |
                                summary_domains$DB == "CC" |
                                summary_domains$DB == "BP" ),]
test$Secondary_Domain <- factor(test$Secondary_Domain, 
                                           levels = rev(unique(test$Secondary_Domain)))
test$Secondary_Domain <- droplevels(test$Secondary_Domain)


#MOSAIC PLOTS GO ENRICHMENTS -FEB24


disease_raw_dataset <- enrichments_domains[which(enrichments_domains$DB == "BP" |
                                                   enrichments_domains$DB == "MF"  | 
                                                   enrichments_domains$DB == "CC"),]

disease_counts <- enrichments_domains[which(enrichments_domains$DB == "BP" |
                                              enrichments_domains$DB == "MF"  | 
                                              enrichments_domains$DB == "CC"),] %>%
  group_by(Secondary_Domain, Method) %>%
  count(Secondary_Domain, Method, description) 

#View(disease_counts[which(disease_counts$Method == "RERS"),])


#View(disease_counts)
#View(disease_merged)
disease_dataset <- left_join(disease_raw_dataset, disease_counts) %>%
  group_by(Secondary_Domain, Method) %>%
  slice_max(n, n=1) %>%
  slice_min(pValue, n= 1) %>%
  slice(1) %>%
  mutate(enrichment_median = median(enrichmentRatio))
library(ggmosaic)
disease_dataset$Secondary_Domain <- as.factor(gsub("_", " ", disease_dataset$Secondary_Domain))
disease_dataset$Secondary_Domain <- factor(disease_dataset$Secondary_Domain, levels=c("Vestibular/Visual Systems",rev(
  levels(disease_dataset$Secondary_Domain))))
disease_dataset$Secondary_Domain[which(disease_dataset$Secondary_Domain == "Vestibular and Visual Systems")] <-
  "Vestibular/Visual Systems"
disease_dataset$Secondary_Domain <- factor(disease_dataset$Secondary_Domain, levels = rev(unique(disease_dataset$Secondary_Domain)))

g <- ggplot(data =  disease_dataset) +
  geom_mosaic(aes(x = product(Secondary_Domain), fill = Method,
                  weight = as.numeric(as.character(enrichment_median))),
              na.rm = TRUE, show.legend = FALSE) +
  scale_fill_manual(values = c("darkgoldenrod3", "steelblue4")) +
  coord_flip() +
  theme_minimal() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size=8),
        plot.title = element_text(hjust = 0.5)) 

#merge x,y dataframe
temp <- ggplot_build(g)$data %>%
  as.data.frame %>% 
  mutate(Method = ifelse(fill == "darkgoldenrod3", "CAAS", "RERS")) %>%
  transmute(x__Secondary_Domain = x__Secondary_Domain,
            Method = Method,
            xmin = xmin,
            xmax = xmax,
            ymin = ymin,
            ymax = ymax,
            fill = fill,
            x.position = (xmax + xmin) / 2,
            y.position = (ymax + ymin) / 2) 
MERGED <- left_join(disease_dataset, temp, by=c("Secondary_Domain" = "x__Secondary_Domain",
                                                "Method" = "Method" ))

library(dplyr)
JITTERED_DATA <- MERGED %>%
  ungroup() %>%
  dplyr::select(Method, Secondary_Domain) %>%
  group_by(Method, Secondary_Domain) %>%
  count(Method, Secondary_Domain) %>%
  mutate(jitter = ifelse(n>1, "Yes", "No"))

JITTERED_MERGED <- left_join(MERGED, JITTERED_DATA,
                             by=c("Secondary_Domain" = "Secondary_Domain",
                                  "Method" = "Method"))

TO_JITTER <- JITTERED_MERGED[which(JITTERED_MERGED$jitter == "Yes"),]
NOT_TO_JITTER <- JITTERED_MERGED[which(JITTERED_MERGED$jitter == "No"),]

levels(TO_JITTER$description) <- c(levels(TO_JITTER$description),
                                   "Weight",
                                   "Neuronal\nmigration",
                                   "Respiratory infect",
                                   "Spleen",
                                   "Hepatomegaly",
                                   "Voice abnormal")
levels(NOT_TO_JITTER$description) <- c(levels(NOT_TO_JITTER$description),
                                       "Weight",
                                       "Neuronal\nmigration",
                                       "Resp infection",
                                       "Spleen",
                                       "Hepatomegaly",
                                       "Voice abnormal")

NOT_TO_JITTER$description[which(as.character(NOT_TO_JITTER$description) == "Failure to gain weight")] <- "Weight"
TO_JITTER$description[which(as.character(TO_JITTER$description)== "Failure to gain weight")] <- "Weight"
NOT_TO_JITTER$description[which(as.character(NOT_TO_JITTER$description) == "Abnormality of neuronal migration")] <- "Neuronal\nmigration"
TO_JITTER$description[which(as.character(TO_JITTER$description)== "Abnormality of neuronal migration")] <- "Neuronal\nmigration"
TO_JITTER$description[which(as.character(TO_JITTER$description)== "Abnormality of neuronal migration")] <- "Neuronal\nmigration"
NOT_TO_JITTER$description[which(as.character(NOT_TO_JITTER$description) == "Abnormality of neuronal migration")] <- "Neuronal\nmigration"
TO_JITTER$description[which(as.character(TO_JITTER$description)== "Recurrent respiratory infections")] <- "Resp infection"
NOT_TO_JITTER$description[which(as.character(NOT_TO_JITTER$description) == "Recurrent respiratory infections")] <- "Resp infection"
TO_JITTER$description[which(as.character(TO_JITTER$description)== "Splenomegaly")] <- "Spleen"
NOT_TO_JITTER$description[which(as.character(NOT_TO_JITTER$description) == "Splenomegaly")] <- "Spleen"
TO_JITTER$description[which(as.character(TO_JITTER$description)== "Abnormality of the voice")] <- "Voice abnormal"
NOT_TO_JITTER$description[which(as.character(NOT_TO_JITTER$description) == "Abnormality of the voice")] <- "Voice abnormal"
NOT_TO_JITTER$needs_rotating <- rep(FALSE, nrow(NOT_TO_JITTER))
NOT_TO_JITTER$label_size <- rep(2.5, nrow(NOT_TO_JITTER))
NOT_TO_JITTER$needs_rotating[NOT_TO_JITTER$size == 219] <- TRUE
NOT_TO_JITTER$label_size[NOT_TO_JITTER$size == 219] <- 2

#View(NOT_TO_JITTER)

gg_FIG1 <- g +  geom_text_repel(data = TO_JITTER , aes(x = TO_JITTER$x.position,
                                                       y = TO_JITTER$y.position, 
                                                       label = TO_JITTER$description), size=3,
                                segment.color = "transparent",
                                fontface="bold",
                                color="white")+
  geom_text(data = NOT_TO_JITTER , aes(x = NOT_TO_JITTER$x.position,
                                       y = NOT_TO_JITTER$y.position, 
                                       label = NOT_TO_JITTER$description,
                                       angle = ifelse(needs_rotating, 60,0),
                                       size = as.factor(label_size)),
            segment.line = "transparent",
            fontface="bold",
            color="white") +
  scale_size_manual(values = c(3, 3)) +
  theme(legend.position = "none")

saveRDS(gg_FIG1, "FigureS48")
