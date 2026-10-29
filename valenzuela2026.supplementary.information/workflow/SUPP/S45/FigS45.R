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
Species_power <- read.csv("SPECIES_COUNTS_TRAITS_PASTED.txt", 
                          header=F, sep = "\t")
#head(Species_power)
colnames(Species_power) <- c("Trait", "N_spc")


#LAMBDA VALIDATION
Validation_pvalues_LAMBDA <- read.csv("ALL_LAMBDA_APR23.txt", 
                                      header=F, sep = "\t")

colnames(Validation_pvalues_LAMBDA) <- c("Trait", "Pos", "tailed2pval", "b", "gene")
Validation_pvalues_LAMBDA <- Validation_pvalues_LAMBDA %>%
  mutate(pval = ifelse(b>0, tailed2pval/2, 1-(tailed2pval/2))) %>%
  dplyr::select(Trait, Pos, pval, gene, b)


#BM VALIDATION
Validation_pvalues_BM <- read.csv("ALL_BM_APR23.txt", 
                                  header=F, sep = "\t")

colnames(Validation_pvalues_BM) <- c("Trait", "Pos", "tailed2pval", "b", "gene")
Validation_pvalues_BM <- Validation_pvalues_BM %>%
  mutate(pval = ifelse(b>0, tailed2pval/2, 1-(tailed2pval/2))) %>%
  dplyr::select(Trait, Pos, pval, gene, b)


power <- read.csv("ALL_POWER_APR23.txt", 
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


ALL_DF_UNIQUE$Prop_spc <- as.numeric(as.character(ALL_DF_UNIQUE$top))/as.numeric(as.character(ALL_DF_UNIQUE$total_studied))
ALL_DF_UNIQUE$decile <- cut(as.numeric(as.character(ALL_DF_UNIQUE$Prop_spc)),
                            breaks = seq(0, 1, by = 0.1), na.rm = TRUE)

ALL_DF_UNIQUE_FDR <- ALL_DF_UNIQUE %>%
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



table <-ALL_DF_UNIQUE[which(ALL_DF_UNIQUE$top>5 &
                              ALL_DF_UNIQUE$bottom>5),] %>%
  drop_na() %>%
  dplyr::group_by(Trait, count_contrasts) %>%
  summarise(
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
#View(table_spc)
#View(table_spc)
p3 <- ggplot(data = table_spc) +
  geom_point(aes(x=N_spc, y=log10(Total),
                 size=`% external <0.05`,
                 shape=Encoding,
                 color=as.factor(`Nº contrasts`))) +
  theme_classic() +
  ylab("log10(Nº CAAS)")+
  xlab("Nº spc x trait") +
  scale_color_manual(values = traits_palette,
                     name = "Nº family contrasts") 

saveRDS(p3, "FigureS45")

