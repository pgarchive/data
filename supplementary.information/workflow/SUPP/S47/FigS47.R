##########################################################
##########################################################
##########################################################
##########################################################
##########SORTED HOMOGENEITY PLOT ANALYSIS ###############
##########################################################

#.libPaths()

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
Species_power <- read.csv("../S45/SPECIES_COUNTS_TRAITS_PASTED.txt", 
                          header=F, sep = "\t")
#head(Species_power)
colnames(Species_power) <- c("Trait", "N_spc")


#LAMBDA VALIDATION
Validation_pvalues_LAMBDA <- read.csv("../S45/ALL_LAMBDA_APR23.txt", 
                                      header=F, sep = "\t")

colnames(Validation_pvalues_LAMBDA) <- c("Trait", "Pos", "tailed2pval", "b", "gene")
Validation_pvalues_LAMBDA <- Validation_pvalues_LAMBDA %>%
  mutate(pval = ifelse(b>0, tailed2pval/2, 1-(tailed2pval/2))) %>%
  select(Trait, Pos, pval, gene)


#BM VALIDATION
Validation_pvalues_BM <- read.csv("../S45/ALL_BM_APR23.txt", 
                                  header=F, sep = "\t")

colnames(Validation_pvalues_BM) <- c("Trait", "Pos", "tailed2pval", "b", "gene")
Validation_pvalues_BM <- Validation_pvalues_BM %>%
  mutate(pval = ifelse(b>0, tailed2pval/2, 1-(tailed2pval/2))) %>%
  select(Trait, Pos, pval, gene)


power <- read.csv("../S45/ALL_POWER_APR23.txt", 
                  header=F, sep = "\t")
#head(power)
colnames(power) <- c("gene", "Trait", "Pos", "top", "bottom")

power_summed <- power %>%
  rowwise() %>%
  mutate(total_studied = sum(top, bottom))


#Table specific TOPs

Merged_power <- left_join(power_summed, Species_power, by=c("Trait" = "Trait"))

SUPPTABLE <- read.csv("../S40/Supp2Table3.tsv", sep="\t", header=T)

INFO <- SUPPTABLE %>% 
  distinct(Trait, count_contrasts)
#View(INFO)

Merged_power_contr <- left_join(Merged_power, INFO, by=c("Trait" = "Trait"))

#PLOT WITH % Internal validation


ALL_validations <- left_join(Validation_pvalues_LAMBDA, Validation_pvalues_BM,
                             by=c("Pos"="Pos", "gene"="gene", "Trait"="Trait"))


ALL_DF <- left_join(Merged_power_contr, ALL_validations,
                    by=c("Pos"="Pos", "gene"="gene", "Trait"="Trait"))


ALL_DF_UNIQUE <- ALL_DF %>%
  distinct()

#View(ALL_DF_UNIQUE)
#View(ALL_DF_UNIQUE)
ALL_DF_UNIQUE$Prop_spc <- as.numeric(as.character(ALL_DF_UNIQUE$top))/as.numeric(as.character(ALL_DF_UNIQUE$total_studied))
ALL_DF_UNIQUE$decile <- cut(as.numeric(as.character(ALL_DF_UNIQUE$Prop_spc)),
                            breaks = seq(0, 1, by = 0.1), na.rm = TRUE)

ALL_DF_UNIQUE_FDR <- ALL_DF_UNIQUE[which(ALL_DF_UNIQUE$top>5 &
                                           ALL_DF_UNIQUE$bottom>5),] %>%
  ungroup() %>%
  group_by(Trait, count_contrasts) %>%
  mutate(pval_lambda.adj = p.adjust(pval.x, method='fdr', 
                                    n = length(pval.x[!is.na(pval.x)]))) %>%
  mutate(pval_BM.adj = p.adjust(pval.y, method='fdr', 
                                n = length(pval.y[!is.na(pval.y)]))) %>%
  mutate(pval_lambda.bonf = p.adjust(pval.x, method='bonferroni', 
                                    n = length(pval.x[!is.na(pval.x)]))) %>%
  mutate(pval_BM.bonf = p.adjust(pval.y, method='bonferroni', 
                                n = length(pval.y[!is.na(pval.y)]))) 
  
ALL_DF_UNIQUE_FDR[which(ALL_DF_UNIQUE_FDR$Trait != "Maximum.lifespan..y." &
                          ALL_DF_UNIQUE_FDR$pval_lambda.adj <0.05),]

table <-ALL_DF_UNIQUE[which(ALL_DF_UNIQUE$top>5 &
                              ALL_DF_UNIQUE$bottom>5),] %>%
  drop_na() %>%
  group_by(Trait, count_contrasts) %>%
  summarize(
    'Total' = length(top),
    'lambda<0.05' = length(pval.x[as.numeric(as.character(pval.x))<0.05]),
    'BM<0.05' = length(pval.y[as.numeric(as.character(pval.y))<0.05]),
    '%lambda<0.05' = length(pval.x[as.numeric(as.character(pval.x))<0.05])/length(top),
    '%BM<0.05' = length(pval.y[as.numeric(as.character(pval.y))<0.05])/length(top))


table_spc <- left_join(table, Species_power, by=c("Trait"="Trait"))
#MELTED_DF <- melt(table_spc, id.vars = c("Trait", "N_spc", "count_contrasts"))

#ADD QUALITATIVE DISCR
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

#View(table_spc)
traits_palette <- readRDS("traits_palette")

colnames(table_spc) <- c("Trait", "Nº contrasts", "Total", 
                         "PGLS <0.05", "RRPP <0.05",
                         "% external <0.05",  "% RRPP <0.05",
                         "N_spc", "TRAIT",
                         "Encoding")




#PLOT FOR decile significant test vs conserv Top AA
table <- ALL_DF_UNIQUE_FDR %>%
  drop_na() %>%
  group_by(Trait, count_contrasts) %>%
  summarize(
    'Total' = length(top),
    'lambda-FDR<0.05' = length(pval_lambda.adj[as.numeric(as.character(pval_lambda.adj))<0.05])/length(top),
    'BM-FDR<0.05' = length(pval_BM.adj[as.numeric(as.character(pval_BM.adj))<0.05])/length(top))

table_spc <- left_join(table, Species_power, by=c("Trait"="Trait"))
#MELTED_DF <- melt(table_spc, id.vars = c("Trait", "N_spc", "count_contrasts"))

#ADD QUALITATIVE DISCR
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

#View(table_spc)
traits_palette <- readRDS("traits_palette")

colnames(table_spc) <- c("Trait", "Nº contrasts", "Total", 
                         "% PGLS FDR<0.05",  "% RRPP FDR<0.05",
                         "N_spc", "TRAIT",
                         "Encoding")


#Permulaciones

CAAS_PERM <- read.csv("../../Figures/Fig2/CAAS_FINAL_MAR23.tab", sep="\t", header=F)
colnames(CAAS_PERM) <- c("gene", "Trait", "Pos", "AA", "pval", "scenario", 
                         "top_n", "bottom_n")

Merged_perm_contr <- left_join(INFO, CAAS_PERM, by=c("Trait" = "Trait"))



#unique(Merged_power$Trait[which(is.na(Merged_power$N_spc))])
table_perm <- Merged_perm_contr[which(Merged_perm_contr$Trait %in% 
                                        ALL_DF_UNIQUE$Trait),] %>%
  group_by(Trait, count_contrasts) %>%
  mutate(pval_perm.adj = p.adjust(pval, method='fdr', 
                                  n = length(pval[!is.na(pval)]))) %>%
  summarize(
    'Total' = length(Pos),
    'SignPerm' = length(pval_perm.adj[as.numeric(as.character(pval_perm.adj)) < 0.05])/length(Pos))

table_perm <- left_join(table_perm, Species_power, by=c("Trait"="Trait"))
table_perm$TRAIT <- table_perm$Trait
levels(table_perm$TRAIT) <- c(levels(table_perm$TRAIT) , 
                             c("Colour_vision", "Social_organization",
                               "Mating_system", "DielActivity",
                               "TrophicGuild", "Locomotion",
                               "Habitat", "SocialMatingSystem", "Monopolization", "NatalDispersal"))

table_perm$TRAIT[startsWith(as.character(table_perm$TRAIT), "Colour_vision_")] <- "Colour_vision"
table_perm$TRAIT[startsWith(as.character(table_perm$TRAIT), "Social_organization_")] <- "Social_organization"
table_perm$TRAIT[startsWith(as.character(table_perm$TRAIT), "Mating_system_")] <- "Mating_system"
table_perm$TRAIT[startsWith(as.character(table_perm$TRAIT), "DielActivity_")] <- "DielActivity"
table_perm$TRAIT[startsWith(as.character(table_perm$TRAIT), "TrophicGuild_")] <- "TrophicGuild"
table_perm$TRAIT[startsWith(as.character(table_perm$TRAIT), "Locomotion_")] <- "Locomotion"
table_perm$TRAIT[startsWith(as.character(table_perm$TRAIT), "Habitat_")] <- "Habitat"
table_perm$TRAIT[startsWith(as.character(table_perm$TRAIT), "SocialMatingSystem_")] <- "SocialMatingSystem"
table_perm$TRAIT[startsWith(as.character(table_perm$TRAIT), "Monopolization_")] <- "Monopolization"
table_perm$TRAIT[startsWith(as.character(table_perm$TRAIT), "NatalDispersal_")] <- "NatalDispersal"

table_perm$Encoding <- ifelse(table_perm$TRAIT %in% qualitative_traits,
                             "Qualitative",
                             "Quantitative")

colnames(table_perm) <- c("Trait", "Nº contrasts", "Total", 
                          "% PermuL FDR<0.05", "N_spc",
                         "TRAIT",
                         "Encoding")



#PLOT VALIDATION DECILES

ALL_DF_UNIQUE$TRAIT <- ALL_DF_UNIQUE$Trait
levels(ALL_DF_UNIQUE$TRAIT) <- c(levels(ALL_DF_UNIQUE$TRAIT) , 
                                      c("Colour_vision", "Social_organization",
                                        "Mating_system", "DielActivity",
                                        "TrophicGuild", "Locomotion",
                                        "Habitat", "SocialMatingSystem", "Monopolization", "NatalDispersal"))

ALL_DF_UNIQUE$TRAIT[startsWith(as.character(ALL_DF_UNIQUE$TRAIT), "Colour_vision_")] <- "Colour_vision"
ALL_DF_UNIQUE$TRAIT[startsWith(as.character(ALL_DF_UNIQUE$TRAIT), "Social_organization_")] <- "Social_organization"
ALL_DF_UNIQUE$TRAIT[startsWith(as.character(ALL_DF_UNIQUE$TRAIT), "Mating_system_")] <- "Mating_system"
ALL_DF_UNIQUE$TRAIT[startsWith(as.character(ALL_DF_UNIQUE$TRAIT), "DielActivity_")] <- "DielActivity"
ALL_DF_UNIQUE$TRAIT[startsWith(as.character(ALL_DF_UNIQUE$TRAIT), "TrophicGuild_")] <- "TrophicGuild"
ALL_DF_UNIQUE$TRAIT[startsWith(as.character(ALL_DF_UNIQUE$TRAIT), "Locomotion_")] <- "Locomotion"
ALL_DF_UNIQUE$TRAIT[startsWith(as.character(ALL_DF_UNIQUE$TRAIT), "Habitat_")] <- "Habitat"
ALL_DF_UNIQUE$TRAIT[startsWith(as.character(ALL_DF_UNIQUE$TRAIT), "SocialMatingSystem_")] <- "SocialMatingSystem"
ALL_DF_UNIQUE$TRAIT[startsWith(as.character(ALL_DF_UNIQUE$TRAIT), "Monopolization_")] <- "Monopolization"
ALL_DF_UNIQUE$TRAIT[startsWith(as.character(ALL_DF_UNIQUE$TRAIT), "NatalDispersal_")] <- "NatalDispersal"

ALL_DF_UNIQUE$Encoding <- ifelse(ALL_DF_UNIQUE$TRAIT %in% qualitative_traits,
                                      "Qualitative",
                                      "Quantitative")

table_spc <- ALL_DF_UNIQUE[which(ALL_DF_UNIQUE$top>5 &
                                   ALL_DF_UNIQUE$bottom>5),] %>%
  group_by( Trait, count_contrasts, Encoding, decile) %>%
  drop_na() %>%
  summarize(
    'Total' = length(top),
    'lambda<0.05' = length(pval.x[as.numeric(as.character(pval.x))<0.05 &
                                             !is.na(pval.x)])/length(top))

colnames(table_spc) <- c("Trait", "Nº contrasts", 
                         "Encoding", "decile", "Total", "% PGLS <0.05")

MELTED_DF <- melt(table_spc, id.vars = c("Trait", "Nº contrasts",
                                         "Encoding", "decile"))
#View(MELTED_DF)

SUBSET <- MELTED_DF[which(MELTED_DF$variable ==  "% PGLS <0.05" &
                            MELTED_DF$Encoding == "Quantitative" ),] %>%
  drop_na()


#% validation
p5 <- ggplot(data = SUBSET) +
  geom_bar(aes(x=decile, y=value,
               group=variable,
                 fill=as.factor(`Nº contrasts`)),
           stat="identity",
           position = position_dodge()) +
  theme_classic() +
  ylab("Median proportion external <0.05")+
  xlab("Decile Top AA conservation") +
  scale_fill_manual(values = traits_palette,
                     name = "Nº family contrasts") +
  guides(color = guide_legend(order=1),
         size = guide_legend(order=2),
         shape = guide_legend(order=3)) +
  theme(axis.text.x = element_text(angle=45, vjust = 0.4),
        axis.title.x = element_blank())+
  theme(legend.position = "none") +
  facet_wrap(~Trait)

saveRDS(p5, "FigureS47")
