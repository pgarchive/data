##########################################
#MERGE 2 DATAFRAMES AND CONVERT TO MATRIX#
##########################################

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
library(tidyverse)
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
library(FactoMineR)
library(corrplot)

#I) READ OBJECTS
#FINAL_DATASET_CAAS_OVERLAPS <- readRDS("FINAL_DATASET_CAAS_OVERLAPS_MAR23")
#FINAL_DATASET_RERS_OVERLAPS <- readRDS("FINAL_DATASET_RERS_OVERLAPS_MAR23")

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

FINAL_DATASET_CAAS_OVERLAPS_RAW <- readRDS("DATASET_CAAS_OVERLAPS_MAR24")
FINAL_DATASET_RERS_OVERLAPS <- readRDS("DATASET_RERS_OVERLAPS_MAR24")
#test <- FINAL_DATASET_RERS_OVERLAPS[which(FINAL_DATASET_RERS_OVERLAPS$Trait1 != "SexualSizeDimorphism" &
#                                    FINAL_DATASET_RERS_OVERLAPS$Trait2 != "SexualSizeDimorphism"),]

#saveRDS(test, "DATASET_RERS_OVERLAPS_MAR24")
#FINAL_DATASET_RERS_OVERLAPS[which(FINAL_DATASET_RERS_OVERLAPS$Trait1 == "SexualSizeDimorphism"),]
#test[which(test$Trait2 == "SexualSizeDimorphism"),]

#CONTR
FINAL_DATASET_CAAS_OVERLAPS <- FINAL_DATASET_CAAS_OVERLAPS_RAW[which(!FINAL_DATASET_CAAS_OVERLAPS_RAW$Trait1 %in% traits_with_no_contrast &
                                                                       !FINAL_DATASET_CAAS_OVERLAPS_RAW$Trait2 %in% traits_with_no_contrast),]
unique(FINAL_DATASET_RERS_OVERLAPS$Trait2)
#II) READ Biological domains information

trait_domains <- read.csv("Traits_domains_SPLIT_FEB22.txt", sep = "\t")



traits_palette <- readRDS("traits_palette.rds")
traits_palette <- c(traits_palette, randomColor(4))
colnames(trait_domains) <- c("trait", "Domain", "Secondary_Domain", "Terciary_domain",
                             "Fourth_Domain")
levels(trait_domains$Secondary_Domain) <- c("Biometabolites", "Body Mass", "Brain", "Climate",
                                            "Coloration-Ornamentation", "Conservation",
                                            "Diet", "Energetics", "Habitat Use", "Longevity", "Organs", "Reproduction", "Sexual Selection", 
                                            "Social organization", "Vestibular/Visual Systems")

FINAL_DATASET_RERS_OVERLAPS$TRAIT1 <- FINAL_DATASET_RERS_OVERLAPS$Trait1
FINAL_DATASET_RERS_OVERLAPS$TRAIT2 <- FINAL_DATASET_RERS_OVERLAPS$Trait2

levels(FINAL_DATASET_RERS_OVERLAPS$TRAIT1) <- c(levels(FINAL_DATASET_RERS_OVERLAPS$TRAIT1) ,
                                                c("Colour_vision", "Social_organization",
                                                  "Mating_system", "DielActivity",
                                                  "TrophicGuild", "Locomotion",
                                                  "Habitat", "SocialMatingSystem", "Monopolization", "NatalDispersal"))

levels(FINAL_DATASET_RERS_OVERLAPS$TRAIT2) <- c(levels(FINAL_DATASET_RERS_OVERLAPS$TRAIT2) ,
                                                c("Colour_vision", "Social_organization",
                                                  "Mating_system", "DielActivity",
                                                  "TrophicGuild", "Locomotion",
                                                  "Habitat", "SocialMatingSystem", "Monopolization", "NatalDispersal"))

#III) READ CREATE MATRIX OF RERs overlaps with fdr<0.05
FINAL_DATASET_RERS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait1), "Colour_vision_")] <- "Colour_vision"
FINAL_DATASET_RERS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait1), "Social_organization_")] <- "Social_organization"
FINAL_DATASET_RERS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait1), "Mating_system_")] <- "Mating_system"
FINAL_DATASET_RERS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait1), "DielActivity_")] <- "DielActivity"
FINAL_DATASET_RERS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait1), "TrophicGuild_")] <- "TrophicGuild"
FINAL_DATASET_RERS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait1), "Locomotion_")] <- "Locomotion"
FINAL_DATASET_RERS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait1), "Habitat_")] <- "Habitat"
FINAL_DATASET_RERS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait1), "SocialMatingSystem_")] <- "SocialMatingSystem"
FINAL_DATASET_RERS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait1), "Monopolization_")] <- "Monopolization"
FINAL_DATASET_RERS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait1), "NatalDispersal_")] <- "NatalDispersal"

#trait2
FINAL_DATASET_RERS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait2), "Colour_vision_")] <- "Colour_vision"
FINAL_DATASET_RERS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait2), "Social_organization_")] <- "Social_organization"
FINAL_DATASET_RERS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait2), "Mating_system_")] <- "Mating_system"
FINAL_DATASET_RERS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait2), "DielActivity_")] <- "DielActivity"
FINAL_DATASET_RERS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait2), "TrophicGuild_")] <- "TrophicGuild"
FINAL_DATASET_RERS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait2), "Locomotion_")] <- "Locomotion"
FINAL_DATASET_RERS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait2), "Habitat_")] <- "Habitat"
FINAL_DATASET_RERS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait2), "SocialMatingSystem_")] <- "SocialMatingSystem"
FINAL_DATASET_RERS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait2), "Monopolization_")] <- "Monopolization"
FINAL_DATASET_RERS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_RERS_OVERLAPS$Trait2), "NatalDispersal_")] <- "NatalDispersal"

#BIND DOMAINS

FINAL_DATASET_RERS_OVERLAPS_domain1 <- dplyr::left_join(FINAL_DATASET_RERS_OVERLAPS,
                                                        trait_domains, by=c("TRAIT1" = "trait"))
FINAL_DATASET_RERS_OVERLAPS_domain2 <- dplyr::left_join(FINAL_DATASET_RERS_OVERLAPS_domain1,
                                                        trait_domains, by=c("TRAIT2" = "trait"))


#IV) ARRANGE LOG-PVALUE FOR MAIN TEXT DRAFT SET OF RERS OVERLAPS
FINAL_DATASET_RERS_OVERLAPS_domain2$logpvalue <- -log10(as.numeric(as.character(FINAL_DATASET_RERS_OVERLAPS_domain2$pvalue)))

#View(FINAL_DATASET_RERS_OVERLAPS_domain2)
#test_RERS <- FINAL_DATASET_RERS_OVERLAPS_domain2[which(as.character(FINAL_DATASET_RERS_OVERLAPS_domain2$Secondary_Domain.x) !=
#                                                         as.character(FINAL_DATASET_RERS_OVERLAPS_domain2$Secondary_Domain.y)),]
#
#test_numb <- test_RERS%>%
#  arrange(desc(logpvalue))
#View(test_RERS[,c(1,2,3,9,13,16)])

##############################################################################
#V) FILTER BY FDR AND OBTAIN MATRIX OF MEDIAN PROPORTION OF OVERLAPPING GENES#
##############################################################################

#View(FINAL_DATASET_RERS_OVERLAPS_domain2)
FINAL_DATASET_RERS_OVERLAPS_domain2$fdr <- p.adjust(as.numeric(as.character(FINAL_DATASET_RERS_OVERLAPS_domain2$pvalue)), 
                                                    method = p.adjust.methods, 
                                                    n = length(as.numeric(as.character(FINAL_DATASET_RERS_OVERLAPS_domain2$pvalue))))


#View(FINAL_DATASET_RERS_OVERLAPS_domain2)

totals_RERs <- FINAL_DATASET_RERS_OVERLAPS_domain2 %>%
  group_by(G1 = pmin(as.character(Secondary_Domain.x), 
                     as.character(Secondary_Domain.y)), 
           G2 = pmax(as.character(Secondary_Domain.x), 
                     as.character(Secondary_Domain.y))) %>%
  count() 

#View(FINAL_DATASET_RERS_OVERLAPS_domain2)

sign_RERs <-  FINAL_DATASET_RERS_OVERLAPS_domain2[which(as.numeric(as.character(FINAL_DATASET_RERS_OVERLAPS_domain2$fdr)) < 0.05),] %>%
  group_by(G1 = pmin(as.character(Secondary_Domain.x), 
                     as.character(Secondary_Domain.y)), 
           G2 = pmax(as.character(Secondary_Domain.x), 
                     as.character(Secondary_Domain.y))) %>%
  count()

prop_RERs <- left_join(totals_RERs, sign_RERs, by=c("G1",
                                                    "G2")) %>%
  replace_na(list(n.y = 0)) %>%
  mutate(Perc = n.y/n.x)

#View(prop_RERs)

FINAL_DATASET_RERS_OVERLAPS_domain2 <- FINAL_DATASET_RERS_OVERLAPS_domain2 %>%
  mutate(G1 = pmin(as.character(Secondary_Domain.x), 
            as.character(Secondary_Domain.y)), 
G2 = pmax(as.character(Secondary_Domain.x), 
          as.character(Secondary_Domain.y)))
  
  
FILTERED_RERS <- left_join(FINAL_DATASET_RERS_OVERLAPS_domain2, prop_RERs,
                           by=c("G1",
                                "G2"))
#View(FILTERED_RERS)
#View(FILTERED_RERS)


FILTERED_RERS_CURATED <- FILTERED_RERS[, c("Secondary_Domain.x", "Secondary_Domain.y", 
                                           "Perc.y")]

require(data.table)

#View(FILTERED_RERS_MAT)

#test <- FILTERED_RERS_CURATED %>% group_by(G1 = pmin(as.character(Secondary_Domain.x), 
  #                                                   as.character(Secondary_Domain.y)), 
   #                                        G2 = pmax(as.character(Secondary_Domain.x), 
    #                                                 as.character(Secondary_Domain.y)))
#View(test)

FILTERED_RERS_MAT <- FILTERED_RERS_CURATED %>% group_by(G1 = pmin(as.character(Secondary_Domain.x), 
                                                                  as.character(Secondary_Domain.y)), 
                                                        G2 = pmax(as.character(Secondary_Domain.x), 
                                                                  as.character(Secondary_Domain.y))) %>%
  summarise(Median_overlap = median(as.numeric(as.character(Perc.y))))

FILTERED_RERS_MAT <- as.data.frame(FILTERED_RERS_MAT)
FILTERED_RERS_MAT[, 'G1'] <- as.factor(FILTERED_RERS_MAT[, 'G1'])
FILTERED_RERS_MAT[, 'G2'] <- as.factor(FILTERED_RERS_MAT[, 'G2'])

FILTERED_RERS_MAT$G1 <- factor(FILTERED_RERS_MAT$G1, levels = rev(levels(trait_domains$Secondary_Domain))) 
FILTERED_RERS_MAT$G2 <- factor(FILTERED_RERS_MAT$G2, levels = rev(levels(trait_domains$Secondary_Domain)))
FILTERED_RERS_MAT_COMPLETED <- complete(FILTERED_RERS_MAT, G1, G2, fill=list(RERS=0))


MATRIX_RERS <- with(FILTERED_RERS_MAT_COMPLETED, 
                    tapply(as.numeric(as.character(Median_overlap)), 
                           list(G1, G2), median))
MATRIX_RERS[is.na(MATRIX_RERS)] <- 0

#VI) Convert matrix into dataframe
get_upper_tri <- function(CorMat){
  CorMat[lower.tri(CorMat)]<- NA
  return(CorMat)
}

get_lower_tri <- function(CorMat){
  CorMat[upper.tri(CorMat)]<- NA
  return(CorMat)
}

lower_tri <- get_lower_tri(MATRIX_RERS)
meltlower <- melt(lower_tri, na.rm = T)
#library(tidyverse)

#VII) READ CREATE MATRIX OF CAAS overlaps with fdr<0.05


FINAL_DATASET_CAAS_OVERLAPS$TRAIT1 <- FINAL_DATASET_CAAS_OVERLAPS$Trait1
FINAL_DATASET_CAAS_OVERLAPS$TRAIT2 <- FINAL_DATASET_CAAS_OVERLAPS$Trait2

levels(FINAL_DATASET_CAAS_OVERLAPS$TRAIT1) <- c(levels(FINAL_DATASET_CAAS_OVERLAPS$TRAIT1) ,
                                                c("Colour_vision", "Social_organization",
                                                  "Mating_system", "DielActivity",
                                                  "TrophicGuild", "Locomotion",
                                                  "Habitat", "SocialMatingSystem", "Monopolization", "NatalDispersal"))

levels(FINAL_DATASET_CAAS_OVERLAPS$TRAIT2) <- c(levels(FINAL_DATASET_CAAS_OVERLAPS$TRAIT2) ,
                                                c("Colour_vision", "Social_organization",
                                                  "Mating_system", "DielActivity",
                                                  "TrophicGuild", "Locomotion",
                                                  "Habitat", "SocialMatingSystem", "Monopolization", "NatalDispersal"))

#trait1
FINAL_DATASET_CAAS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait1), "Colour_vision_")] <- "Colour_vision"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait1), "Social_organization_")] <- "Social_organization"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait1), "Mating_system_")] <- "Mating_system"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait1), "DielActivity_")] <- "DielActivity"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait1), "TrophicGuild_")] <- "TrophicGuild"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait1), "Locomotion_")] <- "Locomotion"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait1), "Habitat_")] <- "Habitat"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait1), "SocialMatingSystem_")] <- "SocialMatingSystem"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait1), "Monopolization_")] <- "Monopolization"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT1[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait1), "NatalDispersal_")] <- "NatalDispersal"

#trait2
FINAL_DATASET_CAAS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait2), "Colour_vision_")] <- "Colour_vision"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait2), "Social_organization_")] <- "Social_organization"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait2), "Mating_system_")] <- "Mating_system"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait2), "DielActivity_")] <- "DielActivity"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait2), "TrophicGuild_")] <- "TrophicGuild"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait2), "Locomotion_")] <- "Locomotion"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait2), "Habitat_")] <- "Habitat"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait2), "SocialMatingSystem_")] <- "SocialMatingSystem"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait2), "Monopolization_")] <- "Monopolization"
FINAL_DATASET_CAAS_OVERLAPS$TRAIT2[startsWith(as.character(FINAL_DATASET_CAAS_OVERLAPS$Trait2), "NatalDispersal_")] <- "NatalDispersal"



#BIND DOMAINS

FINAL_DATASET_CAAS_OVERLAPS_domain1 <- dplyr::left_join(FINAL_DATASET_CAAS_OVERLAPS,
                                                        trait_domains, by=c("TRAIT1" = "trait"))
FINAL_DATASET_CAAS_OVERLAPS_domain2 <- dplyr::left_join(FINAL_DATASET_CAAS_OVERLAPS_domain1,
                                                        trait_domains, by=c("TRAIT2" = "trait"))

FINAL_DATASET_CAAS_OVERLAPS_domain2$logpvalue <- -log10(as.numeric(as.character(FINAL_DATASET_CAAS_OVERLAPS_domain2$pvalue)))

#VIII) ARRANGE LOG-PVALUE FOR MAIN TEXT DRAFT SET OF CAAS OVERLAPS

#test_CAAS <- FINAL_DATASET_CAAS_OVERLAPS_domain2[which(as.character(FINAL_DATASET_CAAS_OVERLAPS_domain2$Secondary_Domain.x) !=
#                                                         as.character(FINAL_DATASET_CAAS_OVERLAPS_domain2$Secondary_Domain.y)),]
#traits_with_no_contrast <- c("Basophils.Abs",
#                             "Basophils..",
#                             "Blood.Urea.Nitrogen",
#                             "Bun.Creat.Ratio",
#                             "Cholesterol",
#                             "Diencephalon",
#                             "GGT",
#                             "Globulin",
#                             "Hippocampus",
#                             "Hypothalamus",
#                             "Maximum.lifespan..y.",
#                             "M_Max_M_NeonateWeight_g",
#                             "MOB",
#                             "Piriform.Lobe",
#                             "Schizocortex",
#                             "Ser.Glut.Oxaloacetic.Transamin.",
#                             "Striatum",
#                             "Striatum..incl..NAcc.",
#                             "Thalamus",
#                             "Uric.Acid",
#                             "XII")
#test_numb <- test_CAAS[which(!test_CAAS$Trait1 %in% traits_with_no_contrast &
#                               !test_CAAS$Trait2 %in% traits_with_no_contrast),] %>%
#  
 # arrange(desc(logpvalue))

#View(test_numb[,c(1,2,3,8,9,12,13,16)])
############################################################################
####IX) Filter by fdr p-value and obtain median of nº of overlapping genes#
############################################################################

FINAL_DATASET_CAAS_OVERLAPS_domain2$fdr <- p.adjust(as.numeric(as.character(FINAL_DATASET_CAAS_OVERLAPS_domain2$pvalue)), 
                                                    method = p.adjust.methods, 
                                                    n = length(as.numeric(as.character(FINAL_DATASET_CAAS_OVERLAPS_domain2$pvalue))))


totals_CAAS <- FINAL_DATASET_CAAS_OVERLAPS_domain2 %>%
  group_by(G1 = pmin(as.character(Secondary_Domain.x), 
                     as.character(Secondary_Domain.y)), 
           G2 = pmax(as.character(Secondary_Domain.x), 
                     as.character(Secondary_Domain.y))) %>%
  count() 

sign_CAAS <-  FINAL_DATASET_CAAS_OVERLAPS_domain2[which(as.numeric(as.character(FINAL_DATASET_CAAS_OVERLAPS_domain2$fdr)) < 0.05),] %>%
  group_by(G1 = pmin(as.character(Secondary_Domain.x), 
                     as.character(Secondary_Domain.y)), 
           G2 = pmax(as.character(Secondary_Domain.x), 
                     as.character(Secondary_Domain.y))) %>%
  count()

prop_CAAS <- left_join(totals_CAAS, sign_CAAS, by=c("G1",
                                                    "G2")) %>%
  replace_na(list(n.y = 0)) %>%
  mutate(Perc = n.y/n.x)

#View(prop_CAAS)

FINAL_DATASET_CAAS_OVERLAPS_domain2 <- FINAL_DATASET_CAAS_OVERLAPS_domain2 %>%
  mutate(G1 = pmin(as.character(Secondary_Domain.x), 
                   as.character(Secondary_Domain.y)), 
         G2 = pmax(as.character(Secondary_Domain.x), 
                   as.character(Secondary_Domain.y)))

FILTERED_CAAS <- left_join(FINAL_DATASET_CAAS_OVERLAPS_domain2, 
                           prop_CAAS,
                           by=c("G1",
                                "G2"))

#View(FILTERED_CAAS)

joined_overlap <- bind_rows(FINAL_DATASET_RERS_OVERLAPS_domain2[which(as.numeric(as.character(FINAL_DATASET_RERS_OVERLAPS_domain2$fdr))<0.05),],
                            FINAL_DATASET_CAAS_OVERLAPS_domain2[which(as.numeric(as.character(FINAL_DATASET_CAAS_OVERLAPS_domain2$fdr))<0.05),],
                            .id = 'Method')
joined_overlap$Method <- ifelse(joined_overlap$Method == 2, "CAAS", "RERS")
#write.table(joined_overlap, "joined_overlap_FEB24.txt", 
#            sep="\t", quote=F, row.names=F,
#            col.names=T)

FILTERED_CAAS_CURATED <- FILTERED_CAAS[, c("Secondary_Domain.x", "Secondary_Domain.y", 
                                           "Perc.y")]
require(data.table)

FILTERED_CAAS_MAT <- FILTERED_CAAS_CURATED %>% drop_na() %>%
  group_by(G1 = pmin(as.character(Secondary_Domain.x), 
                     as.character(Secondary_Domain.y)), 
           G2 = pmax(as.character(Secondary_Domain.x), 
                     as.character(Secondary_Domain.y))) %>%
  summarise(Median_overlap = median(as.numeric(as.character(Perc.y))))

#View(FILTERED_CAAS_MAT)
FILTERED_CAAS_MAT <- as.data.frame(FILTERED_CAAS_MAT)
FILTERED_CAAS_MAT[, 'G1'] <- as.factor(FILTERED_CAAS_MAT[, 'G1'])
FILTERED_CAAS_MAT[, 'G2'] <- as.factor(FILTERED_CAAS_MAT[, 'G2'])

FILTERED_CAAS_MAT$G1 <- factor(FILTERED_CAAS_MAT$G1, levels =levels(trait_domains$Secondary_Domain)) 
FILTERED_CAAS_MAT$G2 <- factor(FILTERED_CAAS_MAT$G2, levels = levels(trait_domains$Secondary_Domain)) 


FILTERED_CAAS_MAT_COMPLETED <- complete(FILTERED_CAAS_MAT, G1, G2, fill=list(CAAS =0))
FILTERED_CAAS_MAT_COMPLETED$Median_overlap[is.na(FILTERED_CAAS_MAT_COMPLETED$Median_overlap)] <- 0

MATRIX_CAAS <- with(FILTERED_CAAS_MAT_COMPLETED, 
                    tapply(as.numeric(as.character(Median_overlap)), 
                           list(G1, G2), median))
MATRIX_CAAS[is.na(MATRIX_CAAS)] <- 0

#X) Melt into dataframe
upper_tri <- get_upper_tri(MATRIX_CAAS)
meltupper <- melt(upper_tri, na.rm = T)


#XI) Merge both CAAS and RERs dataframes into the same df
meltlower_filt <- meltlower[which(as.character(meltlower$Var1) != as.character(meltlower$Var2)),]
meltupper_filt <- meltupper[which(as.character(meltupper$Var1) != as.character(meltupper$Var2)),]

colnames(meltlower_filt) <- c("Var1", "Var2","RERs")
colnames(meltupper_filt) <- c("Var1", "Var2","CAAS")

meltupper_filt$cross <- ifelse(as.numeric(as.character(meltupper_filt$CAAS)) == 0, 
                               TRUE,
                               FALSE)
meltlower_filt$cross <- ifelse(as.numeric(as.character(meltlower_filt$RERs)) == 0, 
                               TRUE,
                               FALSE)

#View(FINAL_DATASET_RERS_OVERLAPS_domain2)
#View(FILTERED_RERS_CURATED)
#View(FILTERED_RERS_MAT)


CAAS <- meltupper_filt[!is.na(colnames(meltupper_filt))]
colnames(CAAS) <- c("Var1", "Var2", "Method", "cross")
RERS <- meltlower_filt[!is.na(colnames(meltlower_filt))] %>%
  subset(., select=c(2,1,3,4))
colnames(RERS) <- c("Var1", "Var2", "Method", "cross")
all_cross <- rbind(CAAS, RERS)
levels(all_cross$Var2)
#View(RERS)

#View(FILTERED_CAAS)

#domains sorting
domains_sorting <- c("Biometabolites", "Body Mass", "Brain", "Climate",
                                      "Coloration-Ornamentation", "Conservation",
                                      "Diet", "Energetics", "Habitat Use", "Longevity",
                                      "Organs", "Reproduction", "Sexual Selection", 
                                      "Social organization", "Vestibular/Visual Systems")

#XI) Create heatmap plot and save object
gg3 <- ggplot() +
  labs(x = NULL, y = NULL) + 
  geom_tile(data = meltupper_filt, 
            mapping = aes(Var1, Var2, 
                          fill = CAAS),
            drop=FALSE) +
  scale_fill_gradient(low = "white",
                      high = "darkgoldenrod3") +
  new_scale_fill() +
  geom_tile(data = meltlower_filt,
            mapping = aes(Var2, Var1,
                          fill = RERs)) +
  scale_fill_gradient(low = "white",
                      high = "steelblue4") +
  scale_y_discrete(
    limits = rev(domains_sorting)
  )  +
  geom_segment(data=transform(subset(all_cross, !!cross), Var1=as.numeric(Var1), 
                              Var2=as.numeric(Var2)), 
               aes(x=Var1-.49, xend=Var1+.49, y=16-Var2-.49, yend=16-Var2+.49), 
               color="black", size=0.1) +
  
  theme(plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.text.x = element_text(size = 9, angle = 30,
                                   hjust =1),
        axis.text.y = element_text(size = 9),
        axis.ticks.x =  element_blank(),
        axis.ticks.y = element_blank(),
        legend.key.size = unit(0.8, 'cm'),
        text = element_text(family = "Helvetica"),
        plot.margin = margin(t = 0,  # Top margin
                             r = 0,  # Right margin
                             b = 0,  # Bottom margin
                             l = 0)) +
  coord_fixed() 

saveRDS(gg3, "gg3_FEB24_contr")


