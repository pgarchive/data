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
Species_power <- read.csv("../S40/SPECIES_COUNTS_TRAITS_PASTED.txt", 
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

power <- read.csv("../S40/ALL_POWER_TRAITS_MAR23.txt", 
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

#ALL_validations_NONBM <- left_join(Validation_pvalues_LAMBDA, Validation_pvalues_PHYLOG,
 #                            by=c("Pos"="Pos", "gene"="gene"))


ALL_DF_quant <- left_join(Validation_pvalues_LAMBDA, Merged_power_contr,
                    by=c("Pos"="Pos", "gene"="gene", "Trait"="Trait")) %>%
  distinct()
ALL_DF_qual <- left_join(Validation_pvalues_PHYLOG, Merged_power_contr,
                          by=c("Pos"="Pos", "gene"="gene", "Trait"="Trait")) %>%
  distinct()

#View(ALL_DF_quant)
ALL_DF <- rbind(ALL_DF_quant, ALL_DF_qual)

ALL_DF_UNIQUE <- ALL_DF %>%
  distinct()

ALL_DF_UNIQUE <- ALL_DF_UNIQUE[which(!is.na(ALL_DF_UNIQUE$count_contrasts)),]

ALL_DF_UNIQUE$Prop_spc <- as.numeric(as.character(ALL_DF_UNIQUE$top))/as.numeric(as.character(ALL_DF_UNIQUE$total_studied))
ALL_DF_UNIQUE$decile <- cut(as.numeric(as.character(ALL_DF_UNIQUE$Prop_spc)),
                            breaks = seq(0, 1, by = 0.1), na.rm = TRUE)

ALL_DF_UNIQUE_FDR <- ALL_DF_UNIQUE %>% ungroup() %>% group_by(Trait) %>%
  mutate(pval.adj = p.adjust(pval, method='fdr', 
                                   n = length(pval[!is.na(pval)])))
 
#WITHOUT FDR CORRECTION NUMBER (QUANTITATIVE PGLS)

#head(test)
#nrow(ALL_DF_quant)
table <- ALL_DF_quant[which(!is.na(ALL_DF_quant$pval) &
                              ALL_DF_quant$top>5 &
                              ALL_DF_quant$bottom>5),] %>%
   group_by(count_contrasts) %>%
  summarize(
    'Total' = length(top),
    'lambda<0.05' = length(pval[as.numeric(as.character(pval))<0.05])/length(top),
    'lambda' = length(pval[as.numeric(as.character(pval))<0.05]))

#WITHOUT FDR CORRECTION NUMBER (QUALITATIVE PHYLOLM)
table <- ALL_DF_qual[which(!is.na(ALL_DF_qual$pval) &
                             ALL_DF_qual$top>5 &
                             ALL_DF_qual$bottom>5),] %>%
  group_by(count_contrasts) %>%
  summarize(
    'Total' = length(top),
    'phylolm<0.05' = length(pval[as.numeric(as.character(pval))<0.05])/length(top),
    'phylolm' = length(pval[as.numeric(as.character(pval))<0.05]))

#PLOT FOR decile significant test vs conserv Top AA
table <- ALL_DF_UNIQUE_FDR[which(!is.na(ALL_DF_UNIQUE_FDR$pval) &
                                   ALL_DF_UNIQUE_FDR$top>5 &
                                   ALL_DF_UNIQUE_FDR$bottom>5),] %>%
  group_by(Trait, count_contrasts) %>%
  summarize(
    'Total' = length(top),
    'lambda-FDR<0.05' = length(pval.adj[as.numeric(as.character(pval.adj))<0.05])/length(top))

#View(table)
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

traits_palette <- readRDS("traits_palette")

#View(table_spc)
colnames(table_spc) <- c("Trait", "Nº contrasts", "Total", 
                         "% FDR<0.05",
                         "N_spc", "TRAIT",
                         "Encoding")




CAAS_PERM <- read.csv("../../Figures/Fig2/CAAS_FINAL_MAR23.tab", sep="\t", header=F)
colnames(CAAS_PERM) <- c("gene", "Trait", "Pos", "AA", "pval", "scenario", 
                         "top_n", "bottom_n")

Merged_perm_contr <- left_join(CAAS_PERM, INFO, by=c("Trait" = "Trait"))


unique(setdiff(unique(INFO$Trait), unique(CAAS_PERM$Trait)))

#unique(Merged_power$Trait[which(is.na(Merged_power$N_spc))])
table_perm <- Merged_perm_contr %>%
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

ALL_DF_UNIQUE_FDR$TRAIT <- ALL_DF_UNIQUE_FDR$Trait
levels(ALL_DF_UNIQUE_FDR$TRAIT) <- c(levels(ALL_DF_UNIQUE_FDR$TRAIT) , 
                                      c("Colour_vision", "Social_organization",
                                        "Mating_system", "DielActivity",
                                        "TrophicGuild", "Locomotion",
                                        "Habitat", "SocialMatingSystem", "Monopolization", "NatalDispersal"))

ALL_DF_UNIQUE_FDR$TRAIT[startsWith(as.character(ALL_DF_UNIQUE_FDR$TRAIT), "Colour_vision_")] <- "Colour_vision"
ALL_DF_UNIQUE_FDR$TRAIT[startsWith(as.character(ALL_DF_UNIQUE_FDR$TRAIT), "Social_organization_")] <- "Social_organization"
ALL_DF_UNIQUE_FDR$TRAIT[startsWith(as.character(ALL_DF_UNIQUE_FDR$TRAIT), "Mating_system_")] <- "Mating_system"
ALL_DF_UNIQUE_FDR$TRAIT[startsWith(as.character(ALL_DF_UNIQUE_FDR$TRAIT), "DielActivity_")] <- "DielActivity"
ALL_DF_UNIQUE_FDR$TRAIT[startsWith(as.character(ALL_DF_UNIQUE_FDR$TRAIT), "TrophicGuild_")] <- "TrophicGuild"
ALL_DF_UNIQUE_FDR$TRAIT[startsWith(as.character(ALL_DF_UNIQUE_FDR$TRAIT), "Locomotion_")] <- "Locomotion"
ALL_DF_UNIQUE_FDR$TRAIT[startsWith(as.character(ALL_DF_UNIQUE_FDR$TRAIT), "Habitat_")] <- "Habitat"
ALL_DF_UNIQUE_FDR$TRAIT[startsWith(as.character(ALL_DF_UNIQUE_FDR$TRAIT), "SocialMatingSystem_")] <- "SocialMatingSystem"
ALL_DF_UNIQUE_FDR$TRAIT[startsWith(as.character(ALL_DF_UNIQUE_FDR$TRAIT), "Monopolization_")] <- "Monopolization"
ALL_DF_UNIQUE_FDR$TRAIT[startsWith(as.character(ALL_DF_UNIQUE_FDR$TRAIT), "NatalDispersal_")] <- "NatalDispersal"

ALL_DF_UNIQUE_FDR$Encoding <- ifelse(ALL_DF_UNIQUE_FDR$TRAIT %in% qualitative_traits,
                                      "Qualitative",
                                      "Quantitative")
#View(table_spc)
#View(table_spc)
test <- ALL_DF_UNIQUE_FDR[which(!is.na(ALL_DF_UNIQUE_FDR$pval) &
                                  ALL_DF_UNIQUE_FDR$top>5 &
                                  ALL_DF_UNIQUE_FDR$bottom>5),]
#View(test)

#View(table_spc)
#Extract traits without validable positions
table_spc <- ALL_DF_UNIQUE[which(!is.na(ALL_DF_UNIQUE$pval) &
                                   ALL_DF_UNIQUE$top>5 &
                                   ALL_DF_UNIQUE$bottom>5),] %>%
  group_by(count_contrasts, Trait) %>%
  summarize(
    'Total nº seggregating AA positions' = length(top),
    '%PGLS<0.05' = length(pval[as.numeric(as.character(pval))<0.05])/length(top),
    'PGLS<0.05' = length(pval[as.numeric(as.character(pval))<0.05]))
traits_significant <- table_spc$Trait[table_spc$`%PGLS<0.05` != 0]

#View(table_spc)
table_spc <- ALL_DF_UNIQUE_FDR[which(!is.na(ALL_DF_UNIQUE_FDR$pval) &
                                       ALL_DF_UNIQUE_FDR$top>5 &
                                       ALL_DF_UNIQUE_FDR$bottom>5),]  %>%
                                       #!ALL_DF_UNIQUE_FDR$Trait %in% traits_with_no_contrast),] %>%
  group_by(Trait, decile, count_contrasts, Encoding) %>%
  summarize(
    'Total' = length(top),
    'lambda<0.05' = length(pval.adj[as.numeric(as.character(pval.adj))<0.05 &
                                             !is.na(pval.adj)])/length(top)) %>%
  ungroup() %>%
  group_by(decile, count_contrasts, Encoding) %>%
  summarize(
    'lambda<0.05' = mean(`lambda<0.05`))

#View(table_spc)
colnames(table_spc) <- c("Decile", "Nº contrasts", 
                         "Encoding", "% PGLS FDR<0.05")

MELTED_DF <- melt(table_spc, id.vars = c("Decile", "Nº contrasts",
                                         "Encoding"))

SUBSET <- MELTED_DF[which(MELTED_DF$variable == '% PGLS FDR<0.05'),] %>%
  drop_na()
#View(SUBSET)


#Conservation max
Merged_power_contr$TRAIT <- Merged_power_contr$Trait
levels(Merged_power_contr$TRAIT) <- c(levels(Merged_power_contr$TRAIT) , 
                              c("Colour_vision", "Social_organization",
                                "Mating_system", "DielActivity",
                                "TrophicGuild", "Locomotion",
                                "Habitat", "SocialMatingSystem", "Monopolization", "NatalDispersal"))

Merged_power_contr$TRAIT[startsWith(as.character(Merged_power_contr$TRAIT), "Colour_vision_")] <- "Colour_vision"
Merged_power_contr$TRAIT[startsWith(as.character(Merged_power_contr$TRAIT), "Social_organization_")] <- "Social_organization"
Merged_power_contr$TRAIT[startsWith(as.character(Merged_power_contr$TRAIT), "Mating_system_")] <- "Mating_system"
Merged_power_contr$TRAIT[startsWith(as.character(Merged_power_contr$TRAIT), "DielActivity_")] <- "DielActivity"
Merged_power_contr$TRAIT[startsWith(as.character(Merged_power_contr$TRAIT), "TrophicGuild_")] <- "TrophicGuild"
Merged_power_contr$TRAIT[startsWith(as.character(Merged_power_contr$TRAIT), "Locomotion_")] <- "Locomotion"
Merged_power_contr$TRAIT[startsWith(as.character(Merged_power_contr$TRAIT), "Habitat_")] <- "Habitat"
Merged_power_contr$TRAIT[startsWith(as.character(Merged_power_contr$TRAIT), "SocialMatingSystem_")] <- "SocialMatingSystem"
Merged_power_contr$TRAIT[startsWith(as.character(Merged_power_contr$TRAIT), "Monopolization_")] <- "Monopolization"
Merged_power_contr$TRAIT[startsWith(as.character(Merged_power_contr$TRAIT), "NatalDispersal_")] <- "NatalDispersal"

Merged_power_contr$Encoding <- ifelse(Merged_power_contr$TRAIT %in% qualitative_traits,
                              "Qualitative",
                              "Quantitative")

table_cons<- Merged_power_contr %>%
  group_by(Trait, count_contrasts, Encoding) %>%
  summarize(
    'AbsenceTop' = length(top[as.numeric(as.character(top)) == 0])/length(top),
    'ConservedTop' = length(top[as.numeric(as.character(top)) == as.numeric(as.character(total_studied))])/length(top)) %>%
  ungroup() %>%
  group_by(count_contrasts, Encoding) %>%
  summarize(
    'AbsenceTop' = median(`AbsenceTop`),
    'ConservedTop' = median(`ConservedTop`))


colnames(table_cons) <- c("Nº contrasts", "Encoding",
                          "Absence Top", "Conserved Top")
MELTED_DF <- melt(table_cons, id.vars = c( "Nº contrasts",
                                           "Encoding"))

#colfunc_ANC <- colorRampPalette(c(traits_palette[2], "white"))
p5 <- ggplot(data = MELTED_DF) +
  geom_bar(aes(x=variable, y=as.numeric(as.character(value)),
               fill=as.factor(`Nº contrasts`)),
           stat="identity",
           position = position_dodge()) +
  theme_classic() +
  ylab("median % CAAS positions")+
  xlab("Total conservation direction") +
  scale_fill_manual(values = traits_palette,
                    name = "Nº family contrasts") +
  guides(color = guide_legend(order=1),
         size = guide_legend(order=2),
         shape = guide_legend(order=3)) +
  facet_wrap(~Encoding) +
  theme(legend.position = "none")
  
saveRDS(p5, "FigureS44")
