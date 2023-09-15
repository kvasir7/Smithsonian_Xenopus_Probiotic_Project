### By: Carly Muletz Wolz

### DADA2 pipeline and creating feature table, taxonomy table and sequence file for later analyses


#############  INITIAL processing of files in terminal to get this into dada2 format #############

## From basespace, the files are downloaded with each sample having a folder 
# and within that folder are the forward and reverse reads

## To make it compatible with dada2 and R we need to copy all the fastq files to one folder 
# Then you can delete the old folder that used to hold the reverse and forward reads. 

## You need to navigate to the project folder (most likely called FASTQ_Generation...) in terminal 

#cd /Users/lkgentry/Library/CloudStorage/OneDrive-SmithsonianInstitution/GrayferLabExperiments/XenoProbioticsBothRunsFastQ/XenoProbiotics-371840469/629625305
# Move all the files within each folders up one folder to the FASTQ files, which will then hold all
# of the fastq files in one main folder as opposed to per sample folders
# then remove the folders with nothing in them now, and unzip

# mv -v *ds*/* .     ### move all files within those folders that have L001 in them up one folder (forward and reverse reads) 

# rm -r *ds*/        ### remove all the empty folders

### YOU MUST UNZIP all of the files

# gunzip *_L001*

###############  INSTALL DADA2 if you don't already have it.  ################
### Follow tutorial on how to install (follow 1. and 2.) https://benjjneb.github.io/dada2/dada-installation.html

library(dada2)
packageVersion("dada2")

library(devtools)
devtools::install_github("benjjneb/dada2")

## if updates required, i always update all and say no to from source


########## DADA2 tutorial is very helpful: https://benjjneb.github.io/dada2/tutorial.html ##########


## Run data is here. On home desktop. 

setwd("/Users/lkgentry/OneDrive - Smithsonian Institution/GrayferLabExperiments/XenoProbioticsBothRunsFastQ/XenoProbiotics-371840469/636614131/")


path <- "/Users/lkgentry/OneDrive - Smithsonian Institution/GrayferLabExperiments/XenoProbioticsBothRunsFastQ/XenoProbiotics-371840469/636614131/"
list.files(path)



#FILTER AND TRIM 

# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


#INSPECTION OF QUALITY PROFILES

plotQualityProfile(fnFs[8:16])


plotQualityProfile(fnRs[8:16])


## Quality looks ok, starts dropping on reverse around 150
## Since need to trim forward primer off 19 bp, will do 275. We want 450 bp merged

#FILTERING AND TRIMMING 

filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))


## parameters for filtering data

## Not uncommon to lose 6 to 10k sequences between reads.in and reads.out
##maxee to 2,5 with low quality samples. This loosens the parameters with the reverse reads

#out3 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,190),
                   #   maxN=0, maxEE=c(2,2), trimLeft = 19, trimRight = 23, 
                    #  truncQ=2, rm.phix=TRUE,
                    #  compress=TRUE, multithread=TRUE) 
out2 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(260,190),
                      maxN=0, maxEE=c(2,5), trimLeft = 19, trimRight = 23, 
                      truncQ=2, rm.phix=TRUE,
                      compress=TRUE, multithread=TRUE)
##out1 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(250,190),
                     # maxN=0, maxEE=c(2,5), trimLeft = 19, trimRight = 23, 
                      #truncQ=2, rm.phix=TRUE,
                      #compress=TRUE, multithread=TRUE)

#out4 <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,190),
                     # maxN=0, maxEE=c(2,5), trimLeft = 19, trimRight = 23, 
                     # truncQ=2, rm.phix=TRUE,
                     # compress=TRUE, multithread=TRUE)


head(out2)
str(out2)
mean(out2[,2])
mean(out2[,1])

## an average pf 80% with this is good. First run 16%. Second run 63%
mean(out2[,2])/mean(out2[,1])


## If comes back as zero, you have to remove the main fastq sequencing file to learn errors and start again
table(file.exists(filtFs))
table(file.exists(filtRs))

## Run 1- Manually had to remove PXC100522 sample and NPC1006222

## Run 2- Manually had to remove NXC92722

## NEED to run whichever parameter last so that filtFs and filtRs are from this


## GOING with out2

errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)


plotErrors(errF, nominalQ=TRUE)

#errors look good for this run. Should follow black line closely

derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names

dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

dadaFs[[1]]


mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers)

setwd("/Users/lkgentry/OneDrive - Smithsonian Institution/GrayferLabExperiments/XenoProbioticAnalysisNov2022/CombinedRunFiles/")
      
## Used out 1 parameters
#saveRDS(seqtab, "XengPro2022_run1_seqtab.rds")

#will merge this with file I create with run 2 data. Re-ran all previous code again with run2 files. Skip like 155 
# second run and save with code below

saveRDS(seqtab, "XengPro2022_run2_seqtab.rds")

##Had to install packages and updates below. If already installed, skip to line 164

#if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

#if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
#BiocManager::install("dada2", version = "3.16")

#BiocManager::install("DECIPHER")

library(BiocManager)  
library(DECIPHER)
  
#Do this code later and see what it looks like
  

  
  
# Merge multiple runs
st1 <- readRDS("XengPro2022_run1_seqtab.rds")
st2 <- readRDS("XengPro2022_run2_seqtab.rds")
  
## You get an error message "Duplicated sample names detected in rownames", but this is ok
## dada2 is just letting you know this
st.all <- mergeSequenceTables(st1, st2, repeats = "sum")
  
dim(st.all)


## Distribution of amplicon sizes in bp
table(nchar(getSequences(st.all)))

seqtab.nochim <- removeBimeraDenovo(st.all, method="consensus", multithread=TRUE, verbose=TRUE)

dim(seqtab.nochim)

## Still retain 94.5% of sequences after chimera removal
sum(seqtab.nochim)/sum(st.all)

rowSums(seqtab.nochim)


## This won't work if you start at line 167, which is fine
getN <- function(x) sum(getUniques(x))

## Make sure you change 'out1' to whatever output you end up selecting
track <- cbind(out2, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))

# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

rownames(track) <- sample.names

head(track)


write.csv(track, "dada2_output_XengPro_combined.csv")


## https://benjjneb.github.io/dada2/training.html
## Assign taxonomy is up three directories so that I can use these files for multiple projects
taxa <- assignTaxonomy(seqtab.nochim, "/Users/lkgentry/Documents/Microbiome Analysis/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

taxa <- addSpecies(taxa, "/Users/lkgentry/Documents/Microbiome Analysis/silva_species_assignment_v138.1 (1).fa.gz")



#  inspect the taxonomic assignments:
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)


library(phyloseq); packageVersion("phyloseq")

## combine feature table and taxonomy table in same order
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               tax_table(taxa))
ps

## rename ASVs to numbers
new.names <- paste0("ASV", seq(ntaxa(ps))) # Define new names ASV1, ASV2, ...
seqs <- taxa_names(ps) # Store sequences
names(seqs) <- new.names # Make map from ASV1 to full sequence
taxa_names(ps) <- new.names # Rename to human-friendly format


## convert feature table to matrix
site_species <-as(otu_table(ps), "matrix")

## need to change this to match mapping file later
rownames(site_species)

# not sure what this is doing. needed it for another project. Might be useful for underscores?
samples.out <- rownames(site_species)
rownames(site_species) <- sapply(strsplit(samples.out, "f"), `[`, 1)

rownames(site_species)

## transpose to make a species by site matrix

species_site <- t(site_species)

# taxon table 
tax <- as(tax_table(ps), "matrix")



getwd()
setwd("/Users/lkgentry/OneDrive - Smithsonian Institution/GrayferLabExperiments/XenoProbioticAnalysisNov2022/CombinedRunFiles/")


##Did not know what to do here. Need to clarify with Carly. Decontam is run later and no big contamination
#in this set of samples. Will just use decontam.  
## Write this file out and look at it. Determine what to do with ASVs in Neg Controls
## select whole worksheet and filter by reads in neg controls, 
## MAKE SURE you have selected all columns/rows before you filter

#this is done with decontam in process file. I did not do any of this for this set of data
## Remove ASVs if they are in all negative controls or in one negative control, but almost all samples
## pretty subjective, no guidelines really on what to do here
## Once done making neg control decisions, remove negative controls from file 
## and rename _final.csv
## Can also use 'decontam' package. That is what we are doing now in most cases.




write.csv(species_site, "XENGPRO2022_feature_table_comborun.csv")

write.csv(tax, "XENGPRO2022_taxonomy_comborun.csv")
write.csv(seqs, 'XENGPRO2022_feature_DNAsequences_comborun.csv')

#this replaces seqRFLP function that is no longer available
writeFasta<-function(data, filename){
  fastaLines = c()
  for (rowNum in 1:nrow(data)){
    fastaLines = c(fastaLines, as.character(paste(">", data[rowNum,"name"], sep = "")))
    fastaLines = c(fastaLines,as.character(data[rowNum,"seq"]))  }
  
fileConn<-file(filename)
writeLines(fastaLines, fileConn)
close(fileConn)
}

seqsDf <- as.data.frame(seqs)
seqsDf$ASV <- row.names(seqsDf)

seqsDf <- seqsDf[ ,c(2,1)]

colnames(seqsDf) <- c('name','seq')

writeFasta(seqsDf, "XENGPRO2022_DNAsequences_combo.fasta")



##library(seqRFLP) No longer available, do not run anything below this line. Will filter using decontam 
#in preprocess file
#done with this file. Load preprocess for next steps








