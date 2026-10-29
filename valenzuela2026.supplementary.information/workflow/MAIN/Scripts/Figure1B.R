####################################################################
##########(FIGURE1B) CONTRASTS FIGURE##################
#######Observe total number of extreme contrasts in our set ########
######of Top and Bottom configurations for each primate family #####
####################################################################
####################################################################

#I) READ LIBRARIES AND TRAITS DATASET
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
library(corrplot)


all_traits <- read.csv("PRIMATES_PHENOMIC_DATASET.tsv", 
                       header=T, sep = "\t")


#II) SELECT DESIRED SPECIES

seq_species <- read.csv("names_FINAL.txt",
                        sep = "\t", header = TRUE)
colnames(seq_species)<- c("Tree_Species", "Genome_species", "Trait_species")
covered_species_all <- as.character(seq_species$Trait_species[!is.na(seq_species$Tree_Species)])
covered_species <- as.character(covered_species_all[!is.na(covered_species_all)])
sequenced_traits <- all_traits[which(all_traits$SpeciesBROAD %in% covered_species),]

##III) ADD FAMILY TAXONOMICAL INFORMATION
groups_list <- strsplit(as.character(sequenced_traits$GroupName), "_")
groups <- c()
for(i in 1:length(groups_list)){
  group <- paste(groups_list[[i]][1], "_", groups_list[[i]][2], sep="")
  groups[i] <- group}
groups

#IV) BIND Human quantitative information into the dataset

library(janitor)

human_measurements <- read.csv("human_direct_measurements.tsv", header=TRUE, sep = "\t")
human_transposed <- t(human_measurements) %>% row_to_names(row_number = 1) %>% head(1)
#convert to numeric columns
human_transposed <- as.data.frame(human_transposed) %>% 
  mutate_at(vars(BodyMass_kg,
                 Body.mass..g., 
                 Body_mass,
                 Brain_mass,
                 Neonatal.mass..g.,
                 LitterSz,
                 Interbirth.interval..days.,
                 Gestation..days.,
                 Lactation..days.,
                 Age.of.female.at.first.reproduction..y.,
                 Maximum.lifespan..y.,
                 Maximum.possible.reproductive.success.under.optimal.conditions..rmax.,
                 Albumin,
                 Calcium,
                 Cholesterol,
                 Phosphorus,
                 Potassium,
                 Ser.Glut.Pyruvic.Transaminase,
                 Ser.Glut.Oxaloacetic.Transamin.,
                 Total.Protein,
                 BodyMassMale_kg,
                 BodyMassFemale_kg,
                 MassDimorphism,
                 Offspring.per.litter,
                 Weaning_prop,
                 Alkaline.Phosphatase,
                 Basophils..,
                 Chloride,
                 Hematocrit,
                 Hemoglobin,
                 Lactic.Dehydrogenase,
                 Lymphocytes..,
                 Magnesium,
                 Monocytes..,
                 Neutrophil.Seg..,
                 Red.Blood.Cells,
                 Sodium,
                 Triglycerides,
                 Uric.Acid,
                 White.Blood.Cells,
                 High.Density.Lipoprotein,
                 Iron,
                 BasalMetabolicRate_ml_O2.h,
                 RestingMetabolicRate_kcal.day), ~as.numeric(as.character(.))) 
#convert to integer
human_transposed <- as.data.frame(human_transposed) %>% 
  mutate_at(vars(Relative_tail_length), ~as.integer(as.character(.))) 

sequenced_traits <- dplyr::bind_rows(sequenced_traits, human_transposed)
levels(sequenced_traits$SpeciesBROAD) <- c(levels(sequenced_traits$SpeciesBROAD),"Homo_sapiens")
sequenced_traits[nrow(sequenced_traits), "SpeciesBROAD"] <- "Homo_sapiens"
sequenced_traits$groups <- c(groups, "APE_Hominidae")
#sequenced_traits <- sequenced_traits[which(sequenced_traits$SpeciesBROAD != "Papio_ursinus"),]


#V) Read tables with computed extremes for qualitative and quantitative
#traits run

QUALI_EXTREMES <- read.csv("Qualitative_extremes",
                           header=T, sep = "\t")


QUANTI_EXTREMES <- read.csv("Quantitative_extremes",
                           header=T, sep = "\t")



#VI)  Count contrasts (Families with both Top and Bottom species)
#for each type of trait

#colnames(QUALI_EXTREMES)
contrasts_counts_qual <- QUALI_EXTREMES %>%
  group_by(label, Trait) %>%
  filter(!Trait %in% c("R_MostCommon_LitterSize_1",
                       "R_MostCommon_LitterSize_2",
                       "R_MostCommon_LitterSize_3")) %>%
  filter(!is.na(value)) %>%
  mutate(outliers = case_when(
    value_FG == 1~ "1", 
    value_BG == 0 ~ "0")) %>%
  ungroup() %>%
  select(Family.x, Trait, outliers) %>%
  distinct() %>%
  group_by(Trait, Family.x) %>%
  mutate(contrast = any(outliers == 0)) %>%
  dplyr::filter(outliers == 1 & contrast == TRUE) %>%
  ungroup() 


contrasts_count_quant <- QUANTI_EXTREMES %>% filter(!is.na(value_FG) | !is.na(value_BG)) %>%
  mutate(outliers = case_when(
    !is.na(value_FG)~ "1", 
    !is.na(value_BG) ~ "0")) %>%
  ungroup() %>%
  select(Family.x, variable, outliers) %>%
  distinct() %>%
  group_by(variable, Family.x) %>%
  mutate(contrast = any(outliers == 0)) %>%
  dplyr::filter(outliers == 1 & contrast == TRUE) %>%
  ungroup() 
colnames(contrasts_count_quant) <- c("Family.x", "Trait", "outliers", "contrast")

#VII) Label corresponding biological trait for binarized qualitative categories
levels(contrasts_counts_qual$Trait) <- c(levels(contrasts_counts_qual$Trait) , 
                                         c("Colour_vision", "Social_organization",
                                           "Mating_system", "DielActivity",
                                           "TrophicGuild", "Locomotion",
                                           "Habitat", "SocialMatingSystem", "Monopolization", "NatalDispersal"))
contrasts_counts_qual$Trait[startsWith(as.character(contrasts_counts_qual$Trait), "Colour_vision_")] <- "Colour_vision"
contrasts_counts_qual$Trait[startsWith(as.character(contrasts_counts_qual$Trait), "Social_organization_")] <- "Social_organization"
contrasts_counts_qual$Trait[startsWith(as.character(contrasts_counts_qual$Trait), "Mating_system_")] <- "Mating_system"
contrasts_counts_qual$Trait[startsWith(as.character(contrasts_counts_qual$Trait), "DielActivity_")] <- "DielActivity"
contrasts_counts_qual$Trait[startsWith(as.character(contrasts_counts_qual$Trait), "TrophicGuild_")] <- "TrophicGuild"
contrasts_counts_qual$Trait[startsWith(as.character(contrasts_counts_qual$Trait), "Locomotion_")] <- "Locomotion"
contrasts_counts_qual$Trait[startsWith(as.character(contrasts_counts_qual$Trait), "Habitat_")] <- "Habitat"
contrasts_counts_qual$Trait[startsWith(as.character(contrasts_counts_qual$Trait), "SocialMatingSystem_")] <- "SocialMatingSystem"
contrasts_counts_qual$Trait[startsWith(as.character(contrasts_counts_qual$Trait), "Monopolization_")] <- "Monopolization"
contrasts_counts_qual$Trait[startsWith(as.character(contrasts_counts_qual$Trait), "NatalDispersal_")] <- "NatalDispersal"


#VIII) Merge all contrasts
all_contrasts <- as.data.frame(full_join(contrasts_count_quant, contrasts_counts_qual))


#IX) Read primate family information palette

list_Fams <- strsplit(as.character(all_contrasts$Family.x),
                      "_", fixed = TRUE)
for (i in 1:length(list_Fams)){
  levels( all_contrasts$Family.x) <- c(levels(all_contrasts$Family.x), list_Fams[[i]][2])
  all_contrasts$Family.x[i] <- list_Fams[[i]][2]
}

primates_palette <- readRDS("primate_families_palette.rds")

names(primates_palette) <- levels(factor(groups))
list_Fams_names <- strsplit(as.character(names(primates_palette)),
                            "_", fixed = TRUE)

for (i in 1:length(list_Fams_names)){
  names(primates_palette)[i] <- list_Fams_names[[i]][2]
}

#X) Read biological trait information palette

trait_domains <- read.csv("Traits_domains_sorted_FIG1.txt", sep = "\t")

traits_palette <- readRDS("traits_palette")
list_Primary <- strsplit(as.character(trait_domains$Domain_Secondary_Domain),
                         "_", fixed = TRUE)
for (i in 1:length(list_Primary)){
  trait_domains$Primary_domain[i] <- list_Primary[[i]][1]
}

levels(all_contrasts$Trait) <- c(levels(all_contrasts$Trait), "IUCN_ER")
all_contrasts$Trait[all_contrasts$Trait == "sequenced_traits$IUCN_ER"] <- "IUCN_ER"

domains <- c()
for (i in 1:length(all_contrasts$Trait)){
  if (all_contrasts$Trait[i] %in% trait_domains$Measurement){
  domain_curr <- as.character(trait_domains$Primary_domain[trait_domains$Measurement == 
                                                                             as.character(all_contrasts$Trait[i])])
domains <- c(domains, domain_curr)}
else{
  print(all_contrasts$Trait[i])
}}

#XI) Replace name for life-history traits
all_contrasts$domains <- domains
all_contrasts$domains[all_contrasts$domains == "Life"] <- "Life-history"


#XII) Finally, make corresponding barplot with contrasts

all_contrasts$Family.x <- factor(all_contrasts$Family.x, levels=rev(c("Cercopithecidae", "Hylobatidae", "Hominidae" ,
                                                                      "Callitrichidae", "Aotidae", "Cebidae", "Atelidae",
                                                                  "Pitheciidae", "Tarsiidae", "Lemuridae", 
                                                                  "Indriidae",   "Lepilemuridae", "Cheirogaleidae", 
                                                                  "Daubentoniidae", "Lorisidae", "Galagidae")))

#View(all_contrasts)
p2 <- ggplot(data=all_contrasts) + geom_bar(aes(x=Family.x, fill=domains), stat="count") +
  coord_flip() + scale_fill_manual(values = traits_palette[c(1,2,3,4,6)]) + theme_classic() + theme(legend.title = element_blank()) +
  ylab("Nº binarized traits with family contrasts") + xlab("Family") +
  theme(legend.position = "none",
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12))

saveRDS(p2, "gg2")

