###Analysis of Vandeputte et al data
###Dada 2 pipleline

## Loading libraries##################################################
library(dada2)
library(tidyverse)

## We have 66 healthy controls and 29 Crohn's disease samples.

## Reading in the metadata############################################

## Metadata has disease status as well as avg cell count per data
metadata = read.csv(file.path("Vandputte", "data","41586_2017_BFnature24460_MOESM10_ESM.csv"))

## Renaming the cell counts
metadata = metadata %>%
  dplyr::rename("CellCount" = "Average.cell.count..per.gram.of.frozen.feces.")

## Dada2 pipeline#####################################################


fastq.dir = file.path("Vandputte", "data","fastq","submitted")

## Forward and reverse file names
fnFs = sort(list.files(fastq.dir, pattern = "R1.fastq.gz",full.names=TRUE,recursive = TRUE))
fnRs = sort(list.files(fastq.dir, pattern = "R2.fastq.gz",full.names=TRUE,recursive = TRUE))

## Creating the sample names sample names
sample.names = c(paste0("SC0",seq(1,9)), paste0("SC",seq(10,40)),paste0("DC0",seq(1,9)),paste0("DC",seq(10,95)))

###Filter and trim
filtFs = file.path(fastq.dir,"filtered",paste0(sample.names,"_F_filt.fastq.gz"))
filtRs = file.path(fastq.dir,"filtered",paste0(sample.names,"_R_filt.fastq.gz"))

names(filtFs) = sample.names
names(filtRs) = sample.names

###Filter and trim function
filtered = filterAndTrim(fnFs, filtFs, fnRs, filtRs, maxN=0, maxEE = 4, truncQ = 2, rm.phix = TRUE, compress = TRUE, multithread = FALSE)
head(filtered)

##Learning error rates
errF = learnErrors(filtFs)

errR = learnErrors(filtRs)

###Core algorithm
dadaFs = dada(filtFs, err = errF)
dadaRs = dada(filtRs, err = errR)

##Merging sequences
mergers = mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab = makeSequenceTable(mergers)

seqtab.nochim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

write.csv(seqtab.nochim, file.path("Vandputte", "data", "OTU_nochim.csv"))

taxa = assignTaxonomy(seqtab.nochim, file.path("Vandputte", "taxonomy", "rdp_train_set_18.fa.gz"), tryRC = TRUE)
taxa = addSpecies(taxa, file.path("Vandputte", "taxonomy", "rdp_species_assignment_18.fa.gz"))
write.csv(taxa, file.path("Vandputte", "data", "taxa_assignments.csv"))


