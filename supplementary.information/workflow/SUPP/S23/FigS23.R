##########################################
#####DESCRIPTIVE PLOTS FROM ORTHOLOG######
########CLUSTERS + RESEQUENCED############
##########################################

library(phytools)
library(ggplot2)
library(ggnewscale)
library(ggtree)
packageVersion("phytools")
library(ggstance)
library(dplyr)
library(tidyr)
library(gridExtra)
require(reshape2)


#READ ALL SPECIES FILE
#Here your path for files
#path <- "../../../3rd_year/GENOME_PHENOME_FINAL/src/Quality_alignments_decisions_CAAS/"

#get species names
#list.directories <- list.dirs(path)
#fields <- strsplit(list.directories, "/", fixed = TRUE)
#genes <- c()

#ldf <- list()
#count <- 0
#for (k in 2:length(list.directories)){
 # gene_name <- fields[[3]][14]
 # files <- list.files(list.directories[k])
  #for (j in 1:length(files)){
  # count <- count + 1
  #filename <- list.files(list.directories[k], pattern = "fasta$")[j]
  #gene_name <- strsplit(filename, "\\.")[[1]][1]
  #print(gene_name)
  #genes[count] <- gene_name
  #filepath <- paste(list.directories[k], filename, sep="/")
  #if (class(try(read.csv(filepath, sep ="\t", na.strings = c('NaN', '')))) == "try-error"){
  #  print("repeat gene")
  #}
  #else{
  #spc_df  <- read.csv(filepath, sep ="\t", na.strings = c('NaN', ''), header = FALSE)
  #gene_names <- rep(gene_name, nrow(spc_df))
  #ldf[[count]] <- as.data.frame(cbind(gene_names,
  #                    spc_df))
#}}}

#length(ldf)

#First plot with all genes counts and all species (without family distinction)
#n_species_df <- do.call(rbind, ldf)
#saveRDS(n_species_df, "n_species_df")

n_species_df <- readRDS("n_species_df")

count_species_global<- n_species_df %>%
  group_by(gene_names) %>%
  count(gene_names) %>%
  filter(n<=230)

#View(count_species_global)

#length(unique(count_species_global$gene_names))

Ratio_of_species <- count_species_global %>%
  ungroup() %>%
  count(n)

#View(Ratio_of_species)

#sort(count_species_global$n, decreasing = TRUE)
#cumsum
count_species_global$n <- sort(count_species_global$n)
count_species_global$cumulative <- cumsum(count_species_global$n/sum(count_species_global$n))
count_species_global$n_decreasing <- sort(count_species_global$n, decreasing = TRUE)
count_species_global$inv_cumulative <- cumsum(count_species_global$n_decreasing/sum(count_species_global$n_decreasing))



genes <- count_species_global %>% ungroup() %>%
  count(n) %>% arrange(n)
genes$cumulative <- cumsum(genes$nn/sum(genes$nn))
#genes$inv_cumulative <- cumsum(count_species_global$n_decreasing/sum(count_species_global$n_decreasing))

#View(count_species_global)
#ggplot(data=count_species_global,aes(x=n)) + 
#  geom_histogram(binwidth = 0.01, col="gold")

ymax=max(genes$nn)
allgenes=sum(genes$nn)

P1 <- ggplot(data = genes) + 
  geom_bar(aes(x = n, y = nn),
           fill = "gold", color = "black",
           stat = "identity", size = 0.1) +
  geom_line(aes(x = n, y = cumulative * ymax), 
            col = "darkred", lwd = 1) +
  scale_y_continuous(
    name = "Cumulative percentage of gene alignments [%]",
    breaks = seq(0, 1, 0.05) * ymax,   # mismo sitio que antes
    labels = scales::percent(seq(0, 1, 0.05)), 
    limits = c(0, ymax),
    sec.axis = sec_axis(~ .,
                        name = "Nº of gene alignments")
  ) +
  scale_x_continuous(
    limits = c(0, 238), expand = c(0, 0),
    name = "Total nº of species"
  ) +
  theme(
    axis.title.x = element_text(size = 20),
    axis.title.y = element_text(size = 20),
    axis.text.y  = element_text(size = 10),
    axis.text.x  = element_text(size = 10)
  )
  
  
saveRDS(P1, "FigS23")
