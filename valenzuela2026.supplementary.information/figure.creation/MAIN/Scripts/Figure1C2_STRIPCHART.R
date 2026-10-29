#####################################
#####FIG1C) BODY_MASS extreme########
##selection for our primate phenomic#
################dataset##############
#####################################

#I) READ LIBRARIES AND TRAITS DATASET
library(ape)
library(ggtree)
library(phytools)
library(dplyr)
library(tidyverse)
#install.packages("~/Descargas/phangorn_2.6.3.tar.gz", repos = NULL, type="source")
#I) Create Heatmap of traits and clusters consistency to show possible and comparable analyses on this
#research topic


all_traits <- read.csv("PRIMATES_PHENOMIC_DATASET.tsv", header=T, sep = "\t")


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

-log10(0.05)

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

sequenced_traits$SpeciesBROAD
##V) READ PHYLOGENY INFORMATION AND MATCH TIP NAMES NOMENCLATURE for BODY MASS
#levels(sequenced_traits$SpeciesBROAD) <- c(levels(sequenced_traits$SpeciesBROAD), "Tarsius_bancanus", "Leontocebus_fuscicolis",
#                                           "Plecturocebus_hoffmansi", "Lagothrix_lagothricha") 
#sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Cephalopachus_bancanus"] <- "Tarsius_bancanus"
#sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Leontocebus_fuscicollis"] <- "Leontocebus_fuscicolis"
#sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Plecturocebus_hoffmannsi"] <- "Plecturocebus_hoffmansi"
#sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Lagothrix_lagotricha"] <- "Lagothrix_lagothricha"

#setdiff(sequenced_traits$SpeciesBROAD[!is.na(sequenced_traits$BodyMass_kg)], tree$tip.label)
tree <- read.tree("FINAL_TREE_BL_WITHOUTNODELAB.newick")
species_match <- match(sequenced_traits$SpeciesBROAD[!is.na(sequenced_traits$BodyMass_kg)],
                       tree$tip.label)  

pruned_primates.tree<-drop.tip(tree,tree$tip.label[-species_match])
pruned_primates.tree <- force.ultrametric(pruned_primates.tree)


##VI) PREPARE PHYLOGENETIC PLOT WITH y axis corresponding to sorted species in
##phylogeny

tree_y <-  function(ggtree, data){
  if(!inherits(ggtree, "ggtree"))
    stop("not a ggtree object")
  left_join(dplyr::select(data, label), dplyr::select(ggtree$data, label, y)) %>%
    pull(y)
}

library(randomcoloR)
library(ggtree)
p<- ggtree(pruned_primates.tree)

families <- c()
for(i in 1:length(p$data$label)){
  if (as.character(p$data$label[i]) %in% as.character(sequenced_traits$SpeciesBROAD)){
    family <- as.character(sequenced_traits$groups[sequenced_traits$SpeciesBROAD == as.character(p$data$label[i])])
    families <- c(families, family)}
  else{
    family <- NA
    families <- c(families, family)
  }}

p$data$Family <- families


##VII) PLOT SPREADNESS SPECIES FOR BODY MASS

primate_families_palette<- readRDS("primate_families_palette.rds")
t1 <-  p + geom_tippoint(aes(color=Family), size=1.5,
                         show.legend = FALSE) +
  scale_color_manual(values = primate_families_palette) 



QUANTI_EXTREMES <- read.csv("Quantitative_extremes",
                            header=T, sep = "\t")

Max_lifespan <- QUANTI_EXTREMES[which(QUANTI_EXTREMES$variable == "Maximum.lifespan..y."),]



#ADD FAMILY INFO

list_Fams <- strsplit(as.character(Max_lifespan$Family.x),
                      "_", fixed = TRUE)
for (i in 1:length(list_Fams)){
  levels(Max_lifespan$Family.x) <- c(levels(Max_lifespan$Family.x), 
                                   list_Fams[[i]][2])
  Max_lifespan$Family.x[i] <- list_Fams[[i]][2]
}


Max_lifespan$Family.x <- factor(Max_lifespan$Family.x, levels=rev(c("Cercopithecidae", "Hylobatidae", "Hominidae" ,
                                                              "Callitrichidae", "Aotidae", "Cebidae", "Atelidae",
                                                              "Pitheciidae", "Tarsiidae", "Lemuridae", 
                                                              "Indriidae",   "Lepilemuridae", "Cheirogaleidae", 
                                                              "Daubentoniidae", "Lorisidae", "Galagidae")))


Max_lifespan_summary <-aggregate(value ~ Family.x, median, 
                              data=Max_lifespan)
total_median_LONG <- median(as.numeric(as.character(Max_lifespan$value)))

View(Max_lifespan)

p5 <- ggplot(data =Max_lifespan,
             aes(y = as.numeric(as.character(value)),
                 x=Family.x)) +
  #geom_text(aes(label=label), position = position_jitter(width = 0.001, height=0.001 ))  +
  geom_jitter(
    aes(fill = Family.x,
        # color = Family.x,
        shape = dplyr::case_when(!is.na(condition_FG) ~ "17", 
                                 !is.na(condition_BG) ~ "15",
                                 TRUE ~ "1"),
        size = dplyr::case_when(!is.na(condition_FG) ~ 4, 
                                !is.na(condition_BG) ~ 4,
                                TRUE ~ 3)),
    position = position_jitter(0.4),
    show.legend = F
  ) +
  geom_crossbar(data=Max_lifespan_summary,
                aes(x = Family.x,
                    col=Family.x,
                    ymin = value, 
                    ymax = value),
                size=0.3, width = 2,
                show.legend = F) +
  scale_fill_manual(values = rev(primate_families_palette[c(3,2,1,
                                                         6,4,7,5,8,
                                                         16,13,12,9,
                                                         10,15,11)])) +
  scale_color_manual(values = rev(primate_families_palette[c(3,2,1,
                                                         6,4,7,5,8,
                                                         16,13,12,9,
                                                         10,15,11)])) +
  scale_size_continuous(range = c(2, 3)) +
  scale_shape_manual(values = c(21,22,24)) +
  geom_hline(aes(yintercept=total_median_LONG),
             linetype=5, colour="black") +
  theme_minimal() +
  xlab("Family") +
  ylab(expression("Residual from Maximum lifespan in years (10"^"x"*")")) +
  theme_classic() +
  theme(strip.text = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() 


log10(QUANTI_EXTREMES$value[which(QUANTI_EXTREMES$variable == "BodyMass_kg")])
QUANTI_EXTREMES$value[which(QUANTI_EXTREMES$variable == "BodyMass_kg")]

#indices <- tree_y(p1, NONALLO_EXTREMES[which(NONALLO_EXTREMES$variable == "BodyMass_kg"),])
#test <- NONALLO_EXTREMES[which(NONALLO_EXTREMES$variable == "BodyMass_kg"),]

#order_test <- test[order(indices),]
#order_families <- rev(unique(order_test$Family.x))

saveRDS(p5, "gg5_STRIPCHART")

#QUALITATIVE EXAMPLE

QUALI_EXTREMES <- read.csv("Qualitative_extremes",
                           header=T, sep = "\t")

Gummivore <- QUALI_EXTREMES[which(QUALI_EXTREMES$Trait == "TrophicGuild_Gummivore"),]



#ADD FAMILY INFO

list_Fams <- strsplit(as.character(Gummivore$Family.x),
                      "_", fixed = TRUE)
for (i in 1:length(list_Fams)){
  levels(Gummivore$Family.x) <- c(levels(Gummivore$Family.x), 
                                     list_Fams[[i]][2])
  Gummivore$Family.x[i] <- list_Fams[[i]][2]
}


Gummivore$Family.x <- factor(Gummivore$Family.x, levels=rev(c("Cercopithecidae", "Hylobatidae", "Hominidae" ,
                                                                    "Callitrichidae", "Aotidae", "Cebidae", "Atelidae",
                                                                    "Pitheciidae", "Tarsiidae", "Lemuridae", 
                                                                    "Indriidae",   "Lepilemuridae", "Cheirogaleidae", 
                                                                    "Daubentoniidae", "Lorisidae", "Galagidae")))

p6 <- ggplot(data =Gummivore,
             aes(y = as.numeric(as.character(value)),
                 x=Family.x)) +
  geom_point(
    aes(fill = Family.x,
        # color = Family.x,
        shape = dplyr::case_when(!is.na(value_FG) ~ "17", 
                                 !is.na(value_BG) ~ "15",
                                 TRUE ~ "1"),
        size = dplyr::case_when(!is.na(value_FG) ~ 4, 
                                !is.na(value_BG) ~ 4,
                                TRUE ~ 3)),
    position=position_jitter(width = 0.3, height = 0.1),
    show.legend = F
  ) +
  scale_fill_manual(values = rev(primate_families_palette[c(3,2,1,
                                                            6,4,7,5,8,
                                                            16,13,12,14,9,
                                                            10,15,11)])) +
  scale_color_manual(values = rev(primate_families_palette[c(3,2,1,
                                                             6,4,7,5,8,
                                                             16,13,12,14,9,
                                                             10,15,11)])) +
  scale_size_continuous(range = c(2, 3)) +
  scale_shape_manual(values = c(21, 22,24)) +
  theme_minimal() +
  xlab("Family") +
  ylab("Non-gummivore diet=0/Gummivore diet=1") +
  theme_classic() +
  theme(strip.text = element_blank(),
        axis.ticks.y = element_blank()) +
  coord_flip() 


log10(QUANTI_EXTREMES$value[which(QUANTI_EXTREMES$variable == "BodyMass_kg")])
QUANTI_EXTREMES$value[which(QUANTI_EXTREMES$variable == "BodyMass_kg")]

#indices <- tree_y(p1, NONALLO_EXTREMES[which(NONALLO_EXTREMES$variable == "BodyMass_kg"),])
#test <- NONALLO_EXTREMES[which(NONALLO_EXTREMES$variable == "BodyMass_kg"),]

#order_test <- test[order(indices),]
#order_families <- rev(unique(order_test$Family.x))

saveRDS(p6, "gg6_STRIPCHART")
