#SCRIPT
library(ggplot2)
library(gridExtra)
library(tibble)
library(data.table)
library(ggrepel)
library(dplyr)
library(tidyr)
library(tidyverse)
library(reshape2)


#I) READ RERCONVERGE TABLE INFORMATION

RER_perms <- read.csv("../../Figures/Fig2/ALL_RER_CONVERGE_UNIQUE_MAR23.txt", 
                      header=T, sep = "\t")
colnames(RER_perms) <- c("Rho", "N", "P", "p.adj",	"permpval",
                         "permpvaladj",	"gene",	"trait")

RERS <- RER_perms

#II) READ CAAS TABLE WITH DISCOVERY + PERMULATION P-VALUE

CAAS_permul <- read.csv("../../Figures/Fig2/CAAS_FINAL_MAR23.tab", 
                        header=F, sep = "\t")
nrow(CAAS_permul)

colnames(CAAS_permul) <- c("gene", "trait", "Pos", "aa", "pval", "scenario", "top", "bottom")


#III) READ INTERNAL BM TEST P-VALUES

Internal_BM <- read.csv("../../Figures/Fig2/ALL_BM_OCT23_BETAS.txt", 
                        header=F, sep = "\t")
#Internal_BM <- Internal_BM[,c(1,2,3,5)]
colnames(Internal_BM) <- c("trait", "Pos", "tailed2pval", "coefficient", "gene")

Internal_BM <- Internal_BM %>%
  mutate(pval = ifelse(coefficient>0, tailed2pval/2, 1-(tailed2pval/2))) %>%
  select(trait, Pos, pval, gene)


#IV) READ INTERNAL LAMBDA TEST P-VALUES

Internal_lambda <- read.csv("../../Figures/Fig2/ALL_LAMBDA_OCT23_BETAS.txt", 
                            header=F, sep = "\t")
colnames(Internal_lambda) <- c("trait", "Pos", "tailed2pval", "coefficient", "gene")

Internal_lambda <- Internal_lambda %>%
  mutate(pval = ifelse(coefficient>0, tailed2pval/2, 1-(tailed2pval/2))) %>%
  select(trait, Pos, pval, gene)


#V) CREATE DF1 FOR CAAS AND DF2 FOR RERCONVERGE AND MERGE biological
# trait domain information to datasets

CAAS_sub1 <- left_join(CAAS_permul, Internal_BM,
                       by=c("gene" = "gene", "trait" = "trait", "Pos" = "Pos"))

CAAS_sub2 <- left_join(CAAS_sub1, Internal_lambda,
                  by=c("gene" = "gene", "trait" = "trait", "Pos" = "Pos"))

#power
power <- read.csv("ALL_POWER_TRAITS_MAR23.txt", 
                  header=F, sep = "\t")
#head(power)
colnames(power) <- c("gene", "Trait", "Pos", "top", "bottom")

power_summed <- power %>%
  rowwise() %>%
  mutate(total_studied = sum(top, bottom))

CAAS <- left_join(CAAS_sub2, power_summed,
                  by=c("gene" = "gene", "trait" = "Trait", "Pos" = "Pos"))
CAAS$pval[CAAS$top.y<=5 | CAAS$bottom.y<=5] <- NA
CAAS$pval.y[CAAS$top.y<=5 | CAAS$bottom.y<=5] <- NA


Trait_domains <- read.csv("../../Figures/Fig2/Traits_domains_FIG2.txt", 
                          header=F, sep = "\t")
Trait_secondary_domains <- Trait_domains[,c(1,3)]
colnames(Trait_secondary_domains) <- c("trait", "Domain")

#CAAS addition

CAAS$TRAIT <- CAAS$trait
levels(CAAS$TRAIT) <- c(levels(CAAS$TRAIT) , 
                        c("Colour_vision", "Social_organization",
                          "Mating_system", "DielActivity",
                          "TrophicGuild", "Locomotion",
                          "Habitat", "SocialMatingSystem", "Monopolization", "NatalDispersal"))

CAAS$TRAIT[startsWith(as.character(CAAS$trait), "Colour_vision_")] <- "Colour_vision"
CAAS$TRAIT[startsWith(as.character(CAAS$trait), "Social_organization_")] <- "Social_organization"
CAAS$TRAIT[startsWith(as.character(CAAS$trait), "Mating_system_")] <- "Mating_system"
CAAS$TRAIT[startsWith(as.character(CAAS$trait), "DielActivity_")] <- "DielActivity"
CAAS$TRAIT[startsWith(as.character(CAAS$trait), "TrophicGuild_")] <- "TrophicGuild"
CAAS$TRAIT[startsWith(as.character(CAAS$trait), "Locomotion_")] <- "Locomotion"
CAAS$TRAIT[startsWith(as.character(CAAS$trait), "Habitat_")] <- "Habitat"
CAAS$TRAIT[startsWith(as.character(CAAS$trait), "SocialMatingSystem_")] <- "SocialMatingSystem"
CAAS$TRAIT[startsWith(as.character(CAAS$trait), "Monopolization_")] <- "Monopolization"
CAAS$TRAIT[startsWith(as.character(CAAS$trait), "NatalDispersal_")] <- "NatalDispersal"

CAAS_DOMAIN <-  left_join(CAAS, Trait_secondary_domains,
                          by=c("TRAIT" = "trait"))

CAAS_DOMAIN$pval.x[as.numeric(as.character(CAAS_DOMAIN$pval.x)) == 0] <- 0.001

#RERS Addition
RERS$TRAIT <- RERS$trait
levels(RERS$TRAIT) <- c(levels(RERS$TRAIT) , 
                        c("Colour_vision", "Social_organization",
                          "Mating_system", "DielActivity",
                          "TrophicGuild", "Locomotion",
                          "Habitat", "SocialMatingSystem", "Monopolization", "NatalDispersal"))

RERS$TRAIT[startsWith(as.character(RERS$trait), "Colour_vision_")] <- "Colour_vision"
RERS$TRAIT[startsWith(as.character(RERS$trait), "Social_organization_")] <- "Social_organization"
RERS$TRAIT[startsWith(as.character(RERS$trait), "Mating_system_")] <- "Mating_system"
RERS$TRAIT[startsWith(as.character(RERS$trait), "DielActivity_")] <- "DielActivity"
RERS$TRAIT[startsWith(as.character(RERS$trait), "TrophicGuild_")] <- "TrophicGuild"
RERS$TRAIT[startsWith(as.character(RERS$trait), "Locomotion_")] <- "Locomotion"
RERS$TRAIT[startsWith(as.character(RERS$trait), "Habitat_")] <- "Habitat"
RERS$TRAIT[startsWith(as.character(RERS$trait), "SocialMatingSystem_")] <- "SocialMatingSystem"
RERS$TRAIT[startsWith(as.character(RERS$trait), "Monopolization_")] <- "Monopolization"
RERS$TRAIT[startsWith(as.character(RERS$trait), "NatalDispersal_")] <- "NatalDispersal"

RERS_DOMAIN <-  left_join(RERS, Trait_secondary_domains,
                          by=c("TRAIT" = "trait"))


#VI) Split information for QUANTITATIVE AND QUALITATIVE TRAITS AND
#OBTAIN COUNTS FOR EACH TYPE OF TRAIT

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

#READ CRITERIA NUMBERS (CONTRASTS)
SUPPTABLE <- read.csv("Supp2Table3.tsv", sep="\t", header=T)

INFO <- SUPPTABLE %>% 
  distinct(Trait, count_contrasts)

traits_twocontr <- INFO$Trait[which(INFO$count_contrasts==2)]
traits_onecontr <- INFO$Trait[which(INFO$count_contrasts==1)]

#FILTER NON-TESTABLE POSITION


#####################
#QUALITATIVE TRAITS##
#####################

#PHYLOLM
Validation_pvalues_PHYLOG <- read.csv("../../Figures/Fig2/ALL_PHYLOGISTIC_OCT23_BETAS.txt", 
                                      header=F, sep = "\t")

colnames(Validation_pvalues_PHYLOG) <- c("Trait", "Pos", "tailed2pval", "b", "gene")
Validation_pvalues_PHYLOG <- Validation_pvalues_PHYLOG %>%
  mutate(pval = ifelse(b>0, tailed2pval/2, 1-(tailed2pval/2))) %>%
  select(Trait, Pos, pval, gene)

CAAS_DOMAIN_QUAL <- left_join(CAAS_DOMAIN, Validation_pvalues_PHYLOG,
                       by=c("gene" = "gene", "trait" = "Trait", "Pos" = "Pos"))

CAAS_DOMAIN_QUAL$pval.y.y[CAAS_DOMAIN_QUAL$top.y<=5 | CAAS_DOMAIN_QUAL$bottom.y<=5] <- NA

table_CAAS_QUAL <- CAAS_DOMAIN_QUAL[which(CAAS_DOMAIN_QUAL$TRAIT %in% qualitative_traits),] %>%
  group_by(trait, Domain) %>%
  mutate(pval_perm.adj = p.adjust(pval.x, method='fdr', 
                                  n = length(pval.x[!is.na(pval.x)]))) %>%
  mutate(pval_PHYLO.adj = p.adjust(pval.y.y, method='fdr',
                                n = length(pval.y.y[!is.na(pval.y.y)]))) %>%
  group_by(trait, Domain, gene) %>%
  summarize(
    "Permulations bootstrap" = length(pval_perm.adj[pval_perm.adj<0.05 &
                                                      !is.na(pval_perm.adj)]),
    "phyloLM" = length(pval_PHYLO.adj[pval_PHYLO.adj<0.05&
                                                  !is.na(pval_PHYLO.adj)]))


MELTED_CAAS_GENES_QUAL_ALL <- melt(table_CAAS_QUAL, id.vars=c("trait", "Domain", "gene")) %>%
  group_by(variable, Domain, trait, gene) %>%
  mutate(significant_gene = ifelse(value == 0, "Non-significant", "Significant")) %>%
  filter(significant_gene == "Significant") %>%
  group_by(variable, Domain, trait) %>%
  summarize(N_genes = n_distinct(gene))


table_RERS_QUAL_ALL <- RERS_DOMAIN[which(RERS_DOMAIN$TRAIT %in% qualitative_traits),] %>%
  group_by(Domain, trait) %>%
  summarize( "Parametric test" = length(p.adj[as.numeric(as.character(p.adj))<0.05]),
             "Permulations bootstrap" = length(permpvaladj[as.numeric(as.character(permpvaladj))<0.05]))

unique(table_RERS_QUANT_ALL$trait)

MELTED_RERS_GENES_QUAL_ALL <- melt(table_RERS_QUAL_ALL)
#merge both caas and rers

ALL_INFO_TRAITS <- full_join(MELTED_CAAS_GENES_QUAL_ALL, MELTED_RERS_GENES_QUAL_ALL, by=c("trait" = "trait"))
ALL_INFO_TRAITS$N_genes[is.na(ALL_INFO_TRAITS$N_genes)] <- 0
ALL_INFO_TRAITS$value[is.na(ALL_INFO_TRAITS$value)] <- 0
ALL_INFO_TRAITS<- ALL_INFO_TRAITS[which(ALL_INFO_TRAITS$variable.x != "Brownian Motion-RRPP" |
                                          is.na(ALL_INFO_TRAITS$variable.x)),]

ALL_INFO_TRAITS$trait <- droplevels(ALL_INFO_TRAITS$trait)

ALL_INFO_TRAITS$trait <- factor(ALL_INFO_TRAITS$trait ,
                                levels =sort(levels(ALL_INFO_TRAITS$trait)))


#View(ALL_INFO)
ALL_INFO_TRAITS$variable.x <- factor(ALL_INFO_TRAITS$variable.x, levels=c("phyloLM",
                                                                          "Permulations bootstrap"))

ALL_INFO_TRAITS$Domain.x <- as.factor(gsub("_", " ", ALL_INFO_TRAITS$Domain.x))
ALL_INFO_TRAITS$Domain.x <- factor(ALL_INFO_TRAITS$Domain.x, levels=c("Vestibular/Visual Systems",rev(
  levels(ALL_INFO_TRAITS$Domain.x))))
ALL_INFO_TRAITS$Domain.x[which(ALL_INFO_TRAITS$Domain.x == "Vestibular and Visual Systems")] <-
  "Vestibular/Visual Systems"
#View(ALL_INFO)

#RERS
RERS_df <- ALL_INFO_TRAITS[,c(3,5,6,7)]
RERS_df <- RERS_df[!duplicated(RERS_df),]

#View(ALL_INFO)

gg1 <- ggplot(data=RERS_df) +
  geom_bar(aes(x=trait,
               fill=variable.y, y=log10(as.numeric(as.character(value))+1)),
           stat="identity",
           position = position_dodge(width=1)) +
  scale_fill_manual(values = rev(RERs_palette), na.translate = F) +
  theme_minimal() + coord_flip() + xlab("Biological Domain") +
  scale_y_continuous( trans = 'reverse') +
  ylab("log10(nº median Genes)")  +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=10),
        legend.text = element_text(size=8),
        legend.position="top",
        legend.title = element_text(size=8)) +
  xlab("") + labs(fill = "RERs") + guides(fill = guide_legend(reverse = TRUE))

#unique(ALL_INFO_TRAITS$trait)
#CAAS
#View(ALL_INFO_TRAITS)
CAAS_df <- ALL_INFO_TRAITS[,c(1,2,3,4)]
CAAS_df <- CAAS_df[!duplicated(CAAS_df),]
gg2 <- ggplot(data=CAAS_df) +
  geom_bar(aes(x=trait, 
               fill=variable.x, y=log10(as.numeric(as.character(N_genes)))), 
           stat="identity",
           position = position_dodge(width=1)) +
  scale_fill_manual(values = rev(CAAS_palette), na.translate = F) +
  theme_minimal() + coord_flip() + xlab("Biological domain") +
  theme(axis.text.y = element_text(hjust = 0.5, size=9),
        legend.position="top",
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=8),
        axis.title.x = element_text(size=10),
        axis.title.y = element_blank(),
        legend.text = element_text(size=8),
        legend.title = element_text(size=8)) +
  ylab("log10(nº median Genes)") + labs(fill ="CAAS")



#cowplot::plot_grid(gg2, gg1, ncol = 2, align = "h")
#nrow(MELTED_RERS_GENES_QUAL)

p5 <- grid.arrange(gg1,
             gg2,
             widths=c(0.2,0.3),
             ncol=2)
saveRDS(p5, "FigureS31")
#ggarrange(p5, ncol = 2)
