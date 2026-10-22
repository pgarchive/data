library(ape)
library(ggtree)
library(phytools)
library(dplyr)
#install.packages("~/Descargas/phangorn_2.6.3.tar.gz", repos = NULL, type="source")
#I) Create Heatmap of traits and clusters consistency to show possible and comparable analyses on this
#research topic


all_traits <- read.csv("../../Figures/Fig1/PRIMATES_PHENOMIC_DATASET.tsv", 
                       header=T, sep = "\t")


#II) Select desired species

seq_species <- read.csv("../../Figures/Fig1/names_FINAL.txt",
                        sep = "\t", header = TRUE)
colnames(seq_species)<- c("Tree_Species", "Genome_species", "Trait_species")
covered_species_all <- as.character(seq_species$Trait_species[!is.na(seq_species$Tree_Species)])
covered_species <- as.character(covered_species_all[!is.na(covered_species_all)])
sequenced_traits <- all_traits[which(all_traits$SpeciesBROAD %in% covered_species),]

groups_list <- strsplit(as.character(sequenced_traits$GroupName), "_")
groups <- c()
for(i in 1:length(groups_list)){
  group <- paste(groups_list[[i]][1], "_", groups_list[[i]][2], sep="")
  groups[i] <- group}
groups

#BIND HUMAN INFO

library(janitor)

human_measurements <- read.csv("../../Figures/Fig1/human_direct_measurements.tsv", header=TRUE, sep = "\t")
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
sequenced_traits <- sequenced_traits[which(sequenced_traits$SpeciesBROAD != "Papio_ursinus"),]

#READ TREE
levels(sequenced_traits$SpeciesBROAD) <- c(levels(sequenced_traits$SpeciesBROAD), "Tarsius_bancanus", "Leontocebus_fuscicolis",
                                           "Plecturocebus_hoffmansi", "Lagothrix_lagothricha") 
sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Cephalopachus_bancanus"] <- "Tarsius_bancanus"
sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Leontocebus_fuscicollis"] <- "Leontocebus_fuscicolis"
sequenced_traits$SpeciesBROAD[sequenced_traits$SpeciesBROAD == "Plecturocebus_hoffmannsi"] <- "Plecturocebus_hoffmansi"


tree <- read.tree("../../Figures/Fig1/FINAL_TREE_BL_WITHOUTNODELAB.newick")

############
#FIT MODELS#
############

#######################
#READ NON-ALLOMETRIC CONTINOUS TRAITS#
#######################
allometric_traits <- read.csv("Allometric_traits.txt", header=T, sep = "\t")

results_path <- "../../../results/Metadata/"


domains <- read.csv("../../Figures/Fig1/Traits_domains_sorted_FIG1.txt", sep="\t")
selected_domains <- domains$Measurement[!domains$Domain_Secondary_Domain %in% c("Parasite-Host_Health", "Morphology_Geometric Morphometrics")]
continous_metadata <- read.csv("Metadata_continous_traits_FINAL.tsv",
                               sep="\t")
redundant <- read.csv("Redundant_measurements.tsv", header=FALSE, sep = "\t")

continous_metadata <- continous_metadata[which(!continous_metadata$trait %in% redundant$V1 &
                                                 continous_metadata$trait %in% selected_domains &
                                                 !continous_metadata$trait %in% colnames(allometric_traits) &
                                                 !continous_metadata$trait %in% c("Total.T4",
"Bilirubin.Direct",
"Iron",
"Fibrinogen",
"Total.Solids",
"Rectal.Body.Temperature",
"Bicarbonate",
"High.Density.Lipoprotein",
"Axillary.Body.Temperature",
"Bone.Thickness")
),]
continous_traits <- as.data.frame(cbind(sequenced_traits$SpeciesBROAD,
                                        sequenced_traits$groups,
                                        sequenced_traits[,which(colnames(sequenced_traits) %in% continous_metadata$trait &
                                                                  !colnames(sequenced_traits) %in% c("InterbirthInterval_days") )
                                        ]))
colnames(continous_traits)[1:2] <- c("label", "Family")



library(mvMORPH)
#fit several models
list <- c()
for (i in 3:ncol(continous_traits)){
  trait <- continous_traits[,i]
  species <- continous_traits[,1]
  names(trait) <- species
  trait <- trait[!is.na(trait)]
  species_match <- match(names(trait),
                         tree$tip.label)  
  pruned_primates.tree<-drop.tip(tree,tree$tip.label[-species_match])
  if (class(try(mvOU(pruned_primates.tree,data=trait,model="OU1"))) == "try-error"){
    next
  }
  OUmodel <-mvOU(pruned_primates.tree,data=trait,model="OU1")
  BMModel <-mvBM(pruned_primates.tree,data=trait,model="BM1")
  EBModel <-mvEB(pruned_primates.tree,data=trait)
  results <- list(OUmodel,BMModel,EBModel)
  table_results <- aicw(results, aicc=TRUE)
  table <- as.data.frame(cbind(as.character(table_results$models),
                               as.numeric(as.character(table_results$AIC)),
                               as.numeric(as.character(table_results$diff)),
                               as.numeric(as.character(table_results$wi)),
                               as.numeric(as.character(table_results$aicweights)), 
                               colnames(continous_traits)[i]))
  colnames(table) <- c(names(table_results),   "trait")
  table <- table %>%
    arrange(desc(as.numeric(as.character(aicweights))))
  list[[i-2]] <- table
}

#############################
#Allometric corrected traits#
#############################

#Create Allometric set of traits
allometric_traits$spc[allometric_traits$spc == "Cephalopachus_bancanus"] <- "Tarsius_bancanus"
allometric_traits$spc[allometric_traits$spc == "Leontocebus_fuscicollis"] <- "Leontocebus_fuscicolis"
allometric_traits$spc[allometric_traits$spc == "Plecturocebus_hoffmannsi"] <- "Plecturocebus_hoffmansi"
allometric_traits <- allometric_traits[which(allometric_traits$spc != "Papio_ursinus"),]
sequenced_allometric_traits <- allometric_traits[which(allometric_traits$spc %in% sequenced_traits$SpeciesBROAD),]
#sort(sequenced_traits$spc)

for (i in 1:nrow(sequenced_allometric_traits)){
  if (as.character(sequenced_allometric_traits$spc[i]) == "Homo_sapiens"){
    sequenced_allometric_traits$GroupName[i] <- "APE_Hominidae_Homo_sapiens"
    print("here")
  }
  else{
    sequenced_allometric_traits$GroupName[i] <- as.character(sequenced_traits$GroupName[sequenced_traits$SpeciesBROAD == as.character(sequenced_allometric_traits$spc[i])])
  }}

groups_list <- strsplit(as.character(sequenced_allometric_traits$GroupName), "_")
groups <- c()
for(i in 1:length(groups_list)){
  group <- paste(groups_list[[i]][1], "_", groups_list[[i]][2], sep="")
  groups[i] <- group}
groups

sequenced_allometric_traits<- cbind(sequenced_allometric_traits, groups)
sequenced_allometric_traits <- sequenced_allometric_traits[,which(!colnames(sequenced_allometric_traits) %in% c("SpeciesBROAD"))]
colnames(sequenced_allometric_traits)[length(sequenced_allometric_traits)] <- "Family"


#make list
list_allometric <- c()
for (i in 2:ncol(sequenced_allometric_traits)){
  trait <- sequenced_allometric_traits[,i]
  species <- sequenced_allometric_traits[,1]
  names(trait) <- species
  trait <- trait[!is.na(trait)]
  species_match <- match(names(trait),
                         tree$tip.label)  
  pruned_primates.tree<-drop.tip(tree,tree$tip.label[-species_match])
  if (class(try(mvOU(pruned_primates.tree,data=trait,model="OU1"))) == "try-error"){
    next
  }
  OUmodel <-mvOU(pruned_primates.tree,data=trait,model="OU1")
  OUmodel$alpha
  BMModel <-mvBM(pruned_primates.tree,data=trait,model="BM1")
  EBModel <-mvEB(pruned_primates.tree,data=trait)
  results <- list(OUmodel,BMModel,EBModel)
  table_results <- aicw(results, aicc=TRUE)
  table <- as.data.frame(cbind(as.character(table_results$models),
                               as.numeric(as.character(table_results$AIC)),
                               as.numeric(as.character(table_results$diff)),
                               as.numeric(as.character(table_results$wi)),
                               as.numeric(as.character(table_results$aicweights)), 
                               colnames(sequenced_allometric_traits)[i]))
  colnames(table) <- c(names(table_results),   "trait")
  table <- table %>%
    arrange(desc(as.numeric(as.character(aicweights))))
  list_allometric[[i-1]] <- table
}


#COMBINE LISTS

ALL_LIST <- c(list, list_allometric)

#ALL_TABLE
All_dataframe <- do.call(rbind, ALL_LIST)
All_dataframe[which(All_dataframe$trait == "BodyMass_kg"),]

#table with most probable model
TOP_dataframe <- do.call(rbind, lapply(ALL_LIST, head, 1))

#palette
#AAs_palette <- randomcoloR::distinctColorPalette(k = 21)
AAs_palette <- readRDS("AAs_colors")

library(ggplot2)
p1 <- ggplot(data=TOP_dataframe) +
  geom_histogram(aes(x=models, fill=models),
                 stat="count") +
  scale_fill_manual(values = AAs_palette)

saveRDS(p1, "FigureS22")

#head(All_dataframe)
