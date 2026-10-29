##########################################################
##########################################################
##########################################################
##########################################################
##########FIGURE 2A) FDR count for our set################
#### of quantitative trait measurements with CAAS#########
############# and RERs methodologies######################
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


#I) READ RERCONVERGE TABLE INFORMATION

RERs <- read.csv("ALL_RER_CONVERGE_UNIQUE_MAR23.txt", 
         header=T, sep = "\t")

colnames(RERs) <- c("Rho", "N", "P", "p.adj",	"permpval",
                         "permpvaladj",	"gene",	"trait")


#II) READ CAAS TABLE WITH DISCOVERY + PERMULATION P-VALUE

CAAS_permul <- read.csv("CAAS_FINAL_MAR23.tab", 
                              header=F, sep = "\t")

colnames(CAAS_permul) <- c("gene", "trait", "Pos", "aa", "pval", "scenario", "top", "bottom")

#CAAS_permul_adjusted1$pval[as.numeric(as.character(CAAS_permul_adjusted1$pval)) == 0] <- 0.001
#CAAS_permul_adjusted2 <- CAAS_permul_adjusted1 %>% group_by(trait) %>%
 # mutate(pval_perm.adj2 = p.adjust(pval, method='fdr', 
  #                                n = length(pval[!is.na(pval)])))

#write.table(CAAS_permul_adjusted2, "CAAS_FINAL_MAR23_ADJUSTED2.tab",
#            sep="\t", col.names=F,row.names=F,quote=F)


#III) READ INTERNAL BM TEST P-VALUES

Internal_BM <- read.csv("ALL_BM_OCT23_BETAS.txt", 
                        header=F, sep = "\t")
#Internal_BM <- Internal_BM[,c(1,2,3,5)]
colnames(Internal_BM) <- c("trait", "Pos", "tailed2pval", "coefficient", "gene")


Internal_BM <- Internal_BM %>%
  mutate(pval = ifelse(coefficient>0, tailed2pval/2, 1-(tailed2pval/2))) %>%
  dplyr::select(trait, Pos, pval, gene)


#IV) READ INTERNAL LAMBDA TEST P-VALUES

Internal_lambda <- read.csv("ALL_LAMBDA_OCT23_BETAS.txt", 
                        header=F, sep = "\t")
colnames(Internal_lambda) <- c("trait", "Pos", "tailed2pval", "coefficient", "gene")

Internal_lambda <- Internal_lambda %>%
  mutate(pval = ifelse(coefficient>0, tailed2pval/2, 1-(tailed2pval/2))) %>%
  dplyr::select(trait, Pos, pval, gene, coefficient)

#Internal_lambda_adjusted <- Internal_lambda %>% group_by(trait) %>%
#  mutate(pval.adj = p.adjust(pval, method='fdr', 
#                                  n = length(pval[!is.na(pval)])))


#V) CREATE DF1 FOR CAAS AND DF2 FOR RERCONVERGE AND MERGE biological
# trait domain information to datasets

CAAS_sub1 <- left_join(CAAS_permul, Internal_BM,
                       by=c("gene" = "gene", "trait" = "trait", "Pos" = "Pos"))

CAAS <- left_join(CAAS_sub1, Internal_lambda,
                  by=c("gene" = "gene", "trait" = "trait", "Pos" = "Pos"))

Trait_domains <- read.csv("Traits_domains_FIG2.txt", 
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
RERS <- RERs
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


table_CAAS_QUANT <- CAAS_DOMAIN[which(!CAAS_DOMAIN$TRAIT %in% qualitative_traits &
                                        !CAAS_DOMAIN$trait %in% traits_with_no_contrast),] %>%
  group_by(trait, Domain) %>%
  mutate(pval_perm.adj = p.adjust(pval.x, method='fdr', 
                                  n = length(pval.x[!is.na(pval.x)]))) %>%
  mutate(pval_BM.adj = p.adjust(pval.y, method='fdr',
                                n = length(pval.y[!is.na(pval.y)]))) %>%
  mutate(pval_lambda.adj = p.adjust(pval, method='fdr',
                                    n = length(pval[!is.na(pval)]))) %>%
  group_by(trait, Domain, gene) %>%
  summarize(
    "Permulations bootstrap (fdr<0.05)" = length(pval_perm.adj[pval_perm.adj<0.05 &
                                                      !is.na(pval_perm.adj)]),
    "Brownian Motion-RRPP (fdr<0.05)" = length(pval_BM.adj[pval_BM.adj<0.05&
                                                  !is.na(pval_BM.adj)]),
    "Lambda-PGLS (fdr<0.05)" = length(pval_lambda.adj[pval_lambda.adj<0.05&
                                             !is.na(pval_lambda.adj)]),
    "Permulations bootstrap (p-val<0.05)" = length(pval.x[pval.x<0.05 &
                                                                 !is.na(pval.x)]),
    "Brownian Motion-RRPP (p-val<0.05)" = length(pval.y[pval.y<0.05&
                                                             !is.na(pval.y)]),
    "Lambda-PGLS (p-val<0.05)" = length(pval[pval<0.05&
                                                        !is.na(pval)]))

  #  mutate("Summary_Counts" = case_when(`Permulations bootstrap` >= 1 &
   #                                     `Lambda-PGLS`>= 1 ~ "BOTH",
    #                              `Permulations bootstrap` >= 1 ~ "Permulations bootstrap",
     #                             `Lambda-PGLS` >= 1 ~ "Lambda-PGLS",
      #                         TRUE ~ "Non-significant"))
  
#sort(as.numeric(as.character(CAAS_DOMAIN$pval.x[which(CAAS_DOMAIN$TRAIT == "Hypothalamus")])))
#test <- p.adjust(sort(as.numeric(as.character(CAAS_DOMAIN$pval.x[which(CAAS_DOMAIN$TRAIT == "Hypothalamus")])))
#, method="fdr")
#View(test)

#test <- table_CAAS_QUANT[which(table_CAAS_QUANT$Summary_Counts == "BOTH"),]
#unique(test$trait)
#View(MELTED_CAAS_GENES_QUANT)
#MELTED_CAAS_GENES_QUANT <- table_CAAS_QUANT %>%
 # select("trait", "Domain", "gene", "Summary_Counts") %>%
  #group_by(Domain, trait, Summary_Counts) %>%
  #summarize(N_genes = n_distinct(gene)) %>%
  #ungroup() %>%
  #complete(nesting(Domain, trait),  Summary_Counts, fill=list(N_genes=0)) %>%
  #filter(Summary_Counts != "Non-significant") %>%
  #group_by(Domain, Summary_Counts) %>%
  #summarize(Median_gene_count = median(N_genes))
  

#test1 <- table_CAAS_QUANT$pval.x[which(table_CAAS_QUANT$trait == "Alkaline.Phosphatase")]
#test2 <- table_CAAS_QUANT$pval_perm.adj[which(table_CAAS_QUANT$trait == "Alkaline.Phosphatase")]

#length(test1)
#length(test2)
#length(sort(as.numeric(as.character(test1[test1<0.05]))))
#length(sort(as.numeric(as.character(test2[test2<0.05]))))
#test3 <- p.adjust(test1, method="fdr")
#length(test3[test3<0.05])
#df <- as.data.frame(test1)

#length(test1[test1<0.02])

#ggplot(as.data.frame(test1)) +
 # geom_histogram(aes(x=test1)) +
  #geom_vline(aes(xintercept=0.05))+
  #geom_vline(aes(xintercept=0.02, col="red"))
  

#View(MELTED_CAAS_GENES_QUANT)

MELTED_CAAS_GENES_QUANT <-  melt(table_CAAS_QUANT, id.vars=c("trait", "Domain", "gene")) %>%
 group_by(variable, Domain, trait, gene) %>%
  mutate(significant_gene = ifelse(value == 0, "Non-significant", "Significant")) %>%
  group_by(variable, Domain, trait, significant_gene) %>%
  summarize(N_genes = n_distinct(gene, significant_gene)) %>%
 ungroup() %>%
  complete(nesting(variable, Domain, trait), significant_gene, fill=list(N_genes=0)) %>%
 filter(significant_gene == "Significant") %>% 
  group_by(Domain, variable) %>%
  summarize(Median_gene_count = median(N_genes))

Ns_quant <- table_CAAS_QUANT %>%
  group_by(Domain) %>%
  summarise(count = n_distinct(trait))

MELTED_CAAS_GENES_QUANT <- left_join(MELTED_CAAS_GENES_QUANT, Ns_quant)


table_RERS_QUANT <- RERS_DOMAIN[which(!RERS_DOMAIN$TRAIT %in% qualitative_traits),] %>%
  group_by(Domain, trait) %>%
  summarize( "Parametric test (fdr<0.05)" = length(p.adj[as.numeric(as.character(p.adj))<0.05 &
                                               !is.na(as.numeric(as.character(p.adj))) ]),
             "Permulations bootstrap (fdr<0.05)" = length(permpvaladj[as.numeric(as.character(permpvaladj))<0.05 &
                                                             !is.na(as.numeric(as.character(permpvaladj)))]),
             "Parametric test (p-val<0.05)" = length(P[as.numeric(as.character(P))<0.05 &
                                                           !is.na(as.numeric(as.character(P)))]),
             "Permulations bootstrap (p-val<0.05)" = length(permpval[as.numeric(as.character(permpval))<0.05 &
                                                                        !is.na(as.numeric(as.character(permpval)))]))


MELTED_RERS_GENES_QUANT <- melt(table_RERS_QUANT, id.vars=c("trait", "Domain")) %>%
  group_by(variable, Domain) %>%
  summarize(Median_gene_count = median(value))
#View(MELTED_RERS_GENES_QUANT)

Ns_quant <- table_RERS_QUANT %>%
  group_by(Domain) %>%
  summarise(count = n_distinct(trait))


MELTED_RERS_GENES_QUANT <- left_join(MELTED_RERS_GENES_QUANT, Ns_quant)

#############################################
##VII) Create figure for quantitative trait##
#############################################

colfunc_CAAS <- colorRampPalette(c("darkgoldenrod3", "white"))
CAAS_palette <- colfunc_CAAS(10)[c(1,5)]
colfunc_RERs <- colorRampPalette(c("steelblue4", "white"))
RERs_palette <- colfunc_RERs(10)[c(1,5)]

#test <- table_RERS_QUANT %>%
 # distinct(Domain, trait)
#View(test)

#test2 <- RERS_DOMAIN %>%
#  distinct(Domain, trait)
#View(test2)

#View(ALL_INFO)
ALL_INFO <- full_join(MELTED_CAAS_GENES_QUANT, MELTED_RERS_GENES_QUANT, by=c("Domain" = "Domain"))
ALL_INFO$Median_gene_count.x[is.na(ALL_INFO$Median_gene_count.x)] <- 0
ALL_INFO$count.x[is.na(ALL_INFO$count.x)] <- 0
ALL_INFO$Median_gene_count.y[is.na(ALL_INFO$Median_gene_count.y)] <- 0
ALL_INFO$count.y[is.na(ALL_INFO$count.y)] <- 0
ALL_INFO_FILT <- ALL_INFO %>%
  filter(variable.y %in% c("Permulations bootstrap (fdr<0.05)",
                         "Permulations bootstrap (p-val<0.05)")) %>%
  filter(variable.x %in% c("Permulations bootstrap (fdr<0.05)",
                           "Permulations bootstrap (p-val<0.05)") |
                           is.na(variable.x))
#View(ALL_INFO_FILT)
#ALL_INFO <-  ALL_INFO[which(!ALL_INFO$variable == "Parametric test"),]
ALL_INFO_FILT$variable.y <- factor(ALL_INFO_FILT$variable.y, levels=c(levels(ALL_INFO_FILT$variable.y),
                                                                      "p-val<0.05",
                                                                      "FDR<0.05"))
ALL_INFO_FILT$variable.x <- factor(ALL_INFO_FILT$variable.x, levels=c(levels(ALL_INFO_FILT$variable.x),
                                                                      "p-val<0.05",
                                                                      "FDR<0.05"))


ALL_INFO_FILT$variable.y[ALL_INFO_FILT$variable.y == "Permulations bootstrap (fdr<0.05)"] <- "FDR<0.05"
ALL_INFO_FILT$variable.x[ALL_INFO_FILT$variable.x == "Permulations bootstrap (p-val<0.05)"] <- "p-val<0.05"
ALL_INFO_FILT$variable.x[ALL_INFO_FILT$variable.x == "Permulations bootstrap (fdr<0.05)"] <- "FDR<0.05"
ALL_INFO_FILT$variable.y[ALL_INFO_FILT$variable.y == "Permulations bootstrap (p-val<0.05)"] <- "p-val<0.05"


ALL_INFO_FILT$Domain <- as.factor(gsub("_", " ", ALL_INFO_FILT$Domain))
ALL_INFO_FILT$Domain <- factor(ALL_INFO_FILT$Domain, levels=c("Vestibular/Visual Systems",rev(
                                                 levels(ALL_INFO_FILT$Domain))))
ALL_INFO_FILT$Domain[which(ALL_INFO_FILT$Domain == "Vestibular and Visual Systems")] <-
  "Vestibular/Visual Systems"

#RERS
RERS_df <- ALL_INFO_FILT[,c(1,5,6,7)]
RERS_df <- RERS_df[!duplicated(RERS_df),]
RERS_df <- RERS_df %>%
  arrange(variable.y)
gg1 <- ggplot(data=RERS_df) +
  geom_bar(aes(x=Domain,
               fill=variable.y, y=log10(as.numeric(as.character(Median_gene_count.y)))),
               stat="identity",
           width=1.5,
           position = position_dodge(width = 0)) +
  #scale_fill_manual(values = rev(RERs_palette), na.translate = F) +
  scale_fill_manual(
    values = rev(RERs_palette), na.translate = F,
    labels = c(
      "p-val<0.05" = expression(italic("p")["permul"] < 0.05),
      "FDR<0.05" = expression(italic("p")["permul-FDR"] < 0.05)
    )
  ) +
  theme_minimal() + coord_flip() + xlab("Biological Domain") +
  ylab("log10(nº median Genes)")  +
    geom_text(aes(x=Domain, y = -0.5, fontface="italic",
                  label=paste("n=", count.y, sep="")),
              color = "grey42", fontface="italic") +
  theme(axis.text.y = element_blank(),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=12),
        legend.text = element_text(size=10),
        legend.position="top") +
  scale_y_continuous( trans = 'reverse') +
  xlab("") + labs(fill = "RERs") + guides(fill = guide_legend(reverse = TRUE))

#CAAS
CAAS_df <- ALL_INFO_FILT[,c(1,2,3,4)]
CAAS_df <- CAAS_df[!duplicated(CAAS_df),]
CAAS_df$Median_gene_count.x[CAAS_df$Median_gene_count.x == 0 &
                              CAAS_df$count.x !=0] <- NA

#View(CAAS_df)
CAAS_df <- CAAS_df %>%
  arrange(variable.x)
gg2 <- ggplot(data=CAAS_df) +
  geom_bar(aes(x=Domain, 
               fill=variable.x, y=log10(as.numeric(as.character(Median_gene_count.x)))), 
               stat="identity",
           width=1.5,
           position = position_dodge(width = 0)) +
  #scale_fill_manual(values = rev(CAAS_palette), na.translate = F) +
  scale_fill_manual(
    values = rev(CAAS_palette), na.translate = F,
    labels = c(
      "p-val<0.05" = expression(italic("p")["permul"] < 0.05),
      "FDR<0.05" = expression(italic("p")["permul-FDR"] < 0.05)
    )
  ) +
  theme_minimal() + coord_flip() + xlab("Biological domain") +
  geom_text(aes(x=Domain, y = -0.5, label=paste("n=", count.x, sep="")),
            fontface="italic", color = "grey42") +
  theme(legend.position="top",
        axis.ticks.y = element_blank(),
        axis.text.y = element_text(size=12, hjust = 0.5),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        legend.text = element_text(size=10)) +
  ylab("log10(nº median Genes)") + labs(fill ="CAAS")



p <- grid.arrange(gg1,
             gg2,
             widths=c(0.25,0.3),
             ncol=2)

#saveRDS(p, "CAAS_RERS_GENE_COUNT_MAR23_with_nom_pval")
saveRDS(p, "CAAS_RERS_GENE_COUNT_FEB24_with_nom_pval_contrasts")



#pdf("Fig2A.pdf")
#grid.arrange(gg1,
#             gg2,
#             widths=c(0.2,0.3),
#             ncol=2)
#dev.off()
