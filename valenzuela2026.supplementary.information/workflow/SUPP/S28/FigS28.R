#################

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

Internal_phylolm <- read.csv("../../Figures/Fig2/ALL_PHYLOGISTIC_OCT23_BETAS.txt", 
                             header=F, sep = "\t")
colnames(Internal_phylolm) <- c("trait", "Pos", "tailed2pval", "coefficient", "gene")


Internal_phylolm <- Internal_phylolm %>%
  mutate(pval = ifelse(coefficient>0, tailed2pval/2, 1-(tailed2pval/2))) %>%
  select(trait, Pos, pval, gene)


#V) CREATE DF1 FOR CAAS AND DF2 FOR RERCONVERGE AND MERGE biological
# trait domain information to datasets

CAAS_sub1 <- left_join(CAAS_permul, Internal_BM,
                       by=c("gene" = "gene", "trait" = "trait", "Pos" = "Pos"))

CAAS_sub2 <- left_join(CAAS_sub1, Internal_lambda,
                       by=c("gene" = "gene", "trait" = "trait", "Pos" = "Pos"))

CAAS <- left_join(CAAS_sub2, Internal_phylolm,
                  by=c("gene" = "gene", "trait" = "trait", "Pos" = "Pos"))

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

table_CAAS_QUAL <- CAAS_DOMAIN[which(CAAS_DOMAIN$TRAIT %in% qualitative_traits),] %>%
  #drop_na() %>%
  group_by(trait, Domain) %>%
  mutate(pval_perm.adj = p.adjust(pval.x, method='fdr', 
                                  n = length(pval.x[!is.na(pval.x)]))) %>%
  mutate(pval_phylolm.adj = p.adjust(pval.y.y, method='fdr',
                                    n = length(pval.y.y[!is.na(pval.y.y)]))) %>%
  group_by(trait, Domain, gene) %>%
  summarize(
    "Permulations bootstrap (fdr<0.05)" = length(pval_perm.adj[pval_perm.adj<0.05 &
                                                                 !is.na(pval_perm.adj)]),
    "phyloGLM (fdr<0.05)" = length(pval_phylolm.adj[pval_phylolm.adj<0.05&
                                                     !is.na(pval_phylolm.adj)]))

MELTED_CAAS_GENES_QUAL <- melt(table_CAAS_QUAL, id.vars=c("trait", "Domain", "gene")) %>%
  group_by(variable, Domain, trait, gene) %>%
  mutate(significant_gene = ifelse(value == 0, "Non-significant", "Significant")) %>%
  group_by(variable, Domain, trait, significant_gene) %>%
  summarize(N_genes = n_distinct(gene, significant_gene)) %>%
  ungroup() %>%
  complete(nesting(variable, Domain, trait), significant_gene, fill=list(N_genes=0)) %>%
  filter(significant_gene == "Significant") %>% 
  group_by(Domain, variable) %>%
  summarize(Median_gene_count = median(N_genes))

Ns_qual <- table_CAAS_QUAL %>%
  group_by(Domain) %>%
  summarise(count = n_distinct(trait))

MELTED_CAAS_GENES_QUAL <- left_join(MELTED_CAAS_GENES_QUAL, Ns_qual)

table_RERS_QUAL <- RERS_DOMAIN[which(RERS_DOMAIN$TRAIT %in% qualitative_traits),] %>%
  group_by(Domain, trait) %>%
  summarize( "Parametric test (fdr<0.05)" = length(p.adj[as.numeric(as.character(p.adj))<0.05 &
                                              !is.na(as.numeric(as.character(p.adj)))]),
             "Permulations bootstrap (fdr<0.05)" = length(permpvaladj[as.numeric(as.character(permpvaladj))<0.05 &
                                                             !is.na(as.numeric(as.character(permpvaladj))) ]))

#View(table_RERS_QUAL)
MELTED_RERS_GENES_QUAL <- melt(table_RERS_QUAL, id.vars=c("trait", "Domain")) %>%
  group_by(variable, Domain) %>%
  summarize(Median_gene_count = median(value))

Ns_qual <- table_RERS_QUAL %>%
  group_by(Domain) %>%
  summarise(count = n_distinct(trait))

MELTED_RERS_GENES_QUAL <- left_join(MELTED_RERS_GENES_QUAL, Ns_qual)


#View(MELTED_RERS_GENES)

#############
#QUALITATIVE#
#############

ALL_INFO <- full_join(MELTED_CAAS_GENES_QUAL, MELTED_RERS_GENES_QUAL, by=c("Domain" = "Domain"))
ALL_INFO$Median_gene_count.x[is.na(ALL_INFO$Median_gene_count.x)] <- 0
ALL_INFO$count.x[is.na(ALL_INFO$count.x)] <- 0
ALL_INFO$Median_gene_count.y[is.na(ALL_INFO$Median_gene_count.y)] <- 0
ALL_INFO$count.y[is.na(ALL_INFO$count.y)] <- 0
ALL_INFO <- ALL_INFO[which((ALL_INFO$variable.x== "Permulations bootstrap (fdr<0.05)" &
                              ALL_INFO$variable.y== "Permulations bootstrap (fdr<0.05)" ) |
                             (is.na(ALL_INFO$variable.x) &
                                ALL_INFO$variable.y== "Permulations bootstrap (fdr<0.05)")),]

#ALL_INFO <-  ALL_INFO[which(!ALL_INFO$variable == "Parametric test"),]

#View(ALL_INFO)
ALL_INFO$variable.y <- factor(ALL_INFO$variable.y, levels=c("Permulations bootstrap (fdr<0.05)"
))
ALL_INFO$variable.x <- factor(ALL_INFO$variable.x, levels=c("Permulations bootstrap (fdr<0.05)"
))

#View(ALL_INFO)
#View(ALL_INFO)
colfunc_CAAS <- colorRampPalette(c("darkgoldenrod3", "white"))
CAAS_palette <- colfunc_CAAS(10)[c(7,1)]
colfunc_RERs <- colorRampPalette(c("steelblue4", "white"))
RERs_palette <- colfunc_RERs(10)[c(7,1)]


ALL_INFO$Domain <- as.factor(gsub("_", " ", ALL_INFO$Domain))
ALL_INFO$Domain <- factor(ALL_INFO$Domain, levels=c("Vestibular/Visual Systems",rev(
  levels(ALL_INFO$Domain))))
ALL_INFO$Domain[which(ALL_INFO$Domain == "Vestibular and Visual Systems")] <-
  "Vestibular/Visual Systems"
#View(ALL_INFO)

#RERS
RERS_df <- ALL_INFO[,c(1,5,6,7)]
RERS_df <- RERS_df[!duplicated(RERS_df),]

#View(ALL_INFO)
gg1 <- ggplot(data=RERS_df) +
  geom_bar(aes(x=Domain,
               fill=variable.y, y=as.numeric(as.character(Median_gene_count.y))),
           stat="identity",
           width=0.8,
           position = position_dodge(width=0)) +
  scale_fill_manual(values = rev(RERs_palette), na.translate = F) +
  theme_minimal() + coord_flip() + xlab("Biological Domain") +
  ylab("nº median Genes")  +
  geom_text(aes(x=Domain, y = -2, fontface="italic",
                label=paste("n=", count.y, sep="")),
            color = "grey42", fontface="italic") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size=12),
        legend.position="top") +
  scale_y_continuous( trans = 'reverse') +
  xlab("") + labs(fill = "RERs") + guides(fill = guide_legend(reverse = TRUE))

#CAAS
CAAS_df <- ALL_INFO[,c(1,2,3,4)]
CAAS_df <- CAAS_df[!duplicated(CAAS_df),]
gg2 <- ggplot(data=CAAS_df) +
  geom_bar(aes(x=Domain, 
               fill=variable.x, y=as.numeric(as.character(Median_gene_count.x))), 
           stat="identity",
           width=0.8,
           position = position_dodge(width=0)) +
  scale_fill_manual(values = rev(CAAS_palette), na.translate = F) +
  theme_minimal() + coord_flip() + xlab("Biological domain") +
  geom_text(aes(x=Domain, y = -20, label=paste("n=", count.x, sep="")),
            fontface="italic", color = "grey42") +
  theme(axis.text.y = element_text(hjust = 0.5, size=15),
        legend.position="top",
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        legend.text = element_text(size=12)) +
  ylab("nº median Genes") + labs(fill ="CAAS")



#cowplot::plot_grid(gg2, gg1, ncol = 2, align = "h")
#nrow(MELTED_RERS_GENES_QUAL)
fig2 <- grid.arrange(gg1,
             gg2,
             widths=c(0.25,0.4),
              ncol=2)

saveRDS(fig2, "FigureS28")
    #,
   # 'informative>10spc' = length(top[top>10 & bottom>10]))


