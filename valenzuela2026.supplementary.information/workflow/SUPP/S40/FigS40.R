##########################################################
##########################################################
##########################################################
##########################################################
##########SORTED HOMOGENEITY PLOT ANALYSIS ###############
##########################################################

library(ggplot2)
library(gridExtra)
library(tibble)
library(data.table)
library(ggrepel)
library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)


#SPECIES POWER
Species_power <- read.csv("SPECIES_COUNTS_TRAITS_PASTED.txt", 
                          header=F, sep = "\t")
#head(Species_power)
colnames(Species_power) <- c("Trait", "N_spc")


#LAMBDA VALIDATION
Validation_pvalues_LAMBDA <- read.csv("../../Figures/Fig2/ALL_LAMBDA_OCT23_BETAS.txt", 
                                      header=F, sep = "\t")

colnames(Validation_pvalues_LAMBDA) <- c("Trait", "Pos", "tailed2pval", "b", "gene")
Validation_pvalues_LAMBDA <- Validation_pvalues_LAMBDA %>%
  mutate(pval = ifelse(b>0, tailed2pval/2, 1-(tailed2pval/2))) %>%
  select(Trait, Pos, pval, gene)

#BM VALIDATION
Validation_pvalues_BM <- read.csv("../../Figures/Fig2/ALL_BM_OCT23_BETAS.txt", 
                                  header=F, sep = "\t")

colnames(Validation_pvalues_BM) <- c("Trait", "Pos", "tailed2pval", "b", "gene")
Validation_pvalues_BM <- Validation_pvalues_BM %>%
  mutate(pval = ifelse(b>0, tailed2pval/2, 1-(tailed2pval/2))) %>%
  select(Trait, Pos, pval, gene)

#PHYLOLM
Validation_pvalues_PHYLOG <- read.csv("../../Figures/Fig2/ALL_PHYLOGISTIC_OCT23_BETAS.txt", 
                                  header=F, sep = "\t")

colnames(Validation_pvalues_PHYLOG) <- c("Trait", "Pos", "tailed2pval", "b", "gene")
Validation_pvalues_PHYLOG <- Validation_pvalues_PHYLOG %>%
  mutate(pval = ifelse(b>0, tailed2pval/2, 1-(tailed2pval/2))) %>%
  select(Trait, Pos, pval, gene)

power <- read.csv("ALL_POWER_TRAITS_MAR23.txt", 
                  header=F, sep = "\t")
#head(power)
colnames(power) <- c("gene", "Trait", "Pos", "top", "bottom")

power_summed <- power %>%
  rowwise() %>%
  mutate(total_studied = sum(top, bottom))


#Table specific TOPs

Merged_power <- left_join(power_summed, Species_power, by=c("Trait" = "Trait"))

SUPPTABLE <- read.csv("Supp2Table3.tsv", sep="\t", header=T)

INFO <- SUPPTABLE %>% 
  distinct(Trait, count_contrasts)
#View(INFO)

Merged_power_contr <- left_join(Merged_power, INFO, by=c("Trait" = "Trait"))

#PLOT WITH % Internal validation


ALL_validations <- left_join(Validation_pvalues_LAMBDA, Validation_pvalues_BM,
                             by=c("Pos"="Pos", "gene"="gene", "Trait"="Trait"))



ALL_DF_quant <- left_join(Validation_pvalues_LAMBDA, Merged_power_contr,
                    by=c("Pos"="Pos", "gene"="gene", "Trait"="Trait")) %>%
  distinct()
ALL_DF_qual <- left_join(Validation_pvalues_PHYLOG, Merged_power_contr,
                          by=c("Pos"="Pos", "gene"="gene", "Trait"="Trait")) %>%
  distinct()


ALL_DF_UNIQUE <- rbind(ALL_DF_quant,ALL_DF_qual) %>%
  distinct()

###############################################################################
#MERGE BOTH QUANTITATIVE AND QUALITATIVE - GET COUNTS VALIDABLE AND %VALIDATED#
###############################################################################

#length(unique(ALL_DF_UNIQUE$Trait))
#length(unique(table_spc$Trait))
table_spc <- ALL_DF_UNIQUE[which(!is.na(ALL_DF_UNIQUE$pval) &
                      ALL_DF_UNIQUE$top>5 &
                      ALL_DF_UNIQUE$bottom>5),] %>%
  group_by(count_contrasts, Trait) %>%
  summarize(
    'Total nº seggregating AA positions' = length(top),
    '%PGLS<0.05' = length(pval[as.numeric(as.character(pval))<0.05])/length(top),
    'PGLS<0.05' = length(pval[as.numeric(as.character(pval))<0.05]))

#saveRDS(table_spc, "table_spc.rds")
#View(table_spc)

#ADD ENCODING INFO AND PLOT
table_spc$TRAIT <- table_spc$Trait
levels(table_spc$TRAIT) <- c(levels(table_spc$TRAIT) , 
                             c("Colour_vision", "Social_organization",
                               "Mating_system", "DielActivity",
                               "TrophicGuild", "Locomotion",
                               "Habitat", "SocialMatingSystem", "Monopolization", "NatalDispersal"))

table_spc$TRAIT[startsWith(as.character(table_spc$TRAIT), "Colour_vision_")] <- "Colour_vision"
table_spc$TRAIT[startsWith(as.character(table_spc$TRAIT), "Social_organization_")] <- "Social_organization"
table_spc$TRAIT[startsWith(as.character(table_spc$TRAIT), "Mating_system_")] <- "Mating_system"
table_spc$TRAIT[startsWith(as.character(table_spc$TRAIT), "DielActivity_")] <- "DielActivity"
table_spc$TRAIT[startsWith(as.character(table_spc$TRAIT), "TrophicGuild_")] <- "TrophicGuild"
table_spc$TRAIT[startsWith(as.character(table_spc$TRAIT), "Locomotion_")] <- "Locomotion"
table_spc$TRAIT[startsWith(as.character(table_spc$TRAIT), "Habitat_")] <- "Habitat"
table_spc$TRAIT[startsWith(as.character(table_spc$TRAIT), "SocialMatingSystem_")] <- "SocialMatingSystem"
table_spc$TRAIT[startsWith(as.character(table_spc$TRAIT), "Monopolization_")] <- "Monopolization"
table_spc$TRAIT[startsWith(as.character(table_spc$TRAIT), "NatalDispersal_")] <- "NatalDispersal"

#PLOT
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
                        "IUCN_ER",
                        "IUCN_Pop_Trend",
                        "Locomotion",
                        "L_Pct_ND_InfMort_lt30d",
                        "Mating_system",
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
                        "StableGroups",
                        "Monopolization")

table_spc$Encoding <- ifelse(table_spc$TRAIT %in% qualitative_traits,
                             "Qualitative",
                             "Quantitative")

traits_palette <- readRDS("../../Figures/Fig1/traits_palette")
table_spc <- left_join(table_spc, Species_power, by=c("Trait"="Trait"))

colnames(table_spc) <- c("Nº contrasts","Trait", "Total nº seggregating AA positions", 
                         "% PGLS<0.05", "PGLS<0.05", "TRAIT",
                         "Encoding", "Nº species measured")


#CON LOG10
p2 <- ggplot(data = table_spc) +
  geom_point(aes(x=log10(`Total nº seggregating AA positions`), 
                 y=`% PGLS<0.05`,
                 size=`Nº species measured`,
                 shape=Encoding,
                 color=as.factor(`Nº contrasts`))) +
  theme_classic() +
  ylab("% Internal p-val<0.05")+
  xlab("log10(Total nº seggregating AA positions)") +
  scale_color_manual(values = traits_palette,
                     name = "Nº family contrasts") 

saveRDS(p2, "FigureS40")



