if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")
a
##16s rRNA-seq analysis of CRC metagenomes
#Quality check: R, xml2, httpuv, latticeExtra, dada2
# China Data
fastq_path <- "C:/Users/Admin/Downloads/SraAccessionList/"
list.files(fastq_path)

#Assign Sample Names:
library(dada2)
library(phyloseq)
library(ggplot2)

fnFs1 <- sort(list.files(fastq_path, pattern="_1.fastq", full.names = TRUE))
fnFs1
fnRs1<- sort(list.files(fastq_path, pattern="_2.fastq", full.names = TRUE))
fnRs1
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs1), "_"), `[`, 1)

#sample.names <- sapply(strsplit(basename(fnFs), "\\.") , `[`, 1)
sample.names

#Plot Quality Control charts:
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:4])

#Assign names of the Filtered sequences:
filtFs1 <- file.path(fastq_path, "filter1", paste0(sample.names, "_1_filt.fastq.gz"))
filtRs1 <- file.path(fastq_path, "filter1", paste0(sample.names, "_2_filt.fastq.gz"))
filtFs1
names(filtFs1) <- sample.names
#names(filtRs) <- sample.names
#filtRs
#Trimming of the sequences:
#out <- filterAndTrim(fnFs, filtFs,fnRs,filtRs,truncLen =(240,160)  ,maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE, compress=TRUE, multithread=TRUE, minQ=20, verbose=TRUE)
#out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs,minQ=20, verbose=TRUE)
out # This will spit out the number of the input and output reads.
out3 <- filterAndTrim(fnFs1, filtFs1, fnRs1, filtRs1, truncLen=c(240,160), maxEE=c(2,2), truncQ=2, minQ=20,compress=TRUE)
out3
##Generate the Filtered Profiles
#Perform Quality check and generate new profiles for filtered sequences:
plotQualityProfile(filtFs[1:4])
#plotQualityProfile(filtRs[1:4])

##Error Rate Profiles
#Generate the Estimated Error Rates for Forward and Reverse Sequences
errF1 <- learnErrors(filtFs1, multithread=TRUE)
errR1 <- learnErrors(filtRs1, multithread=TRUE)

#To visualize the estimated error rates:
#plotErrors(errF, nominalQ=TRUE)

#Third step is the dereplication:
derepFs1 <- derepFastq(filtFs1, verbose=TRUE)
derepFs1
derepRs1 <- derepFastq(filtRs1, verbose=TRUE)
derepFs1

##Sample Inference:
#Fourth step is the sample inference:
dadaFs1 <- dada(derepFs1, err=errF1, multithread=TRUE)
dadaRs1 <- dada(derepRs1, err=errR1, multithread=TRUE)

#5th step: Merging paired reads:
#mergers <- mergePairs(dadaFs, derepFs, verbose=TRUE)
#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#mergers
#6th: Last step is sequence table construction:
seqtab1 <- makeSequenceTable(dadaFs1)
seqtab1
# View No. of samples and No. of variants
dim(seqtab1)

# View the distribution of the sequences
table(nchar(getSequences(seqtab1)))

# Then, removing the chimeras: # This step could take long if verbose=TRUE
data1 <- removeBimeraDenovo(seqtab1, method="consensus", multithread=TRUE, verbose=FALSE)
write.table(data1,"C:/Users/Admin/Downloads/SraAccessionList/metadata_PRJNA818796.txt",sep="\t",append = TRUE)
View(data1)
# Check the no. of samples and no. of variants when chimeras are removed
dim(data1)
#Finally, Assigning the taxonomy:


# To see the first lines of taxonomy:
#unname(head(taxa))

# We will construct the phyloseq object :

##Phyloseq step for 16s rRNA analysis:
#cases <- c(rep("SRR15182562", 2), rep("SRR15182563", 2), rep("SRR15182564", 2), rep("SRR15182565", 2), rep("SRR15182566", 2), rep("SRR15182567", 2), rep("SRR15182568", 2), rep("SRR15182569", 2), rep("SRR15182570", 2), rep("SRR15182571", 2), rep("SRR15182572", 2), rep("SRR15182573", 2), rep("SRR15182574", 2), rep("SRR15182575", 2), rep("SRR15182576", 2), rep("SRR15182577", 2), rep("SRR15182578", 2), rep("SRR15182579", 2), rep("SRR15182580", 2), rep("SRR15182581", 2), rep("SRR15182582", 2), rep("SRR15182583", 2), rep("SRR15182584", 2), rep("SRR15182585", 2), rep("SRR15182586", 2), rep("SRR15182587", 2), rep("SRR15182588", 2), rep("SRR15182589", 2), rep("SRR15182590", 2), rep("SRR15182591", 2), rep("SRR15182592", 2), rep("SRR15182593", 2), rep("SRR15182594", 2), rep("SRR15182595", 2), rep("SRR15182596", 2), rep("SRR15182597", 2), rep("SRR15182598", 2), rep("SRR15182599", 2), rep("SRR15182600", 2), rep("SRR15182601", 2), rep("SRR15182602", 2), rep("SRR15182603", 2), rep("SRR15182604", 2), rep("SRR15182605", 2), rep("SRR15182606", 2), rep("SRR15182607", 2), rep("SRR15182608", 2), rep("SRR15182609", 2), rep("SRR15182610", 2), rep("SRR15182611", 2), rep("SRR15182612", 2), rep("SRR15182613", 2), rep("SRR15182614", 2), rep("SRR15182615", 2), rep("SRR15182616", 2), rep("SRR15182617", 2), rep("SRR15182618", 2), rep("SRR15182619", 2), rep("SRR15182620", 2), rep("SRR15182621", 2), rep("SRR15182622", 2), rep("SRR15182623", 2), rep("SRR15182624", 2), rep("SRR15182625", 2), rep("SRR15182626", 2), rep("SRR15182627", 2), rep("SRR15182628", 2), rep("SRR15182629", 2), rep("SRR15182630", 2), rep("SRR15182631", 2))

#df <- data.frame(Samples= sample.names, Status= cases)
#ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), sample_names(df), tax_table(taxa))
#ps <- prune_samples(sample_names(ps) != "Mock", ps)

## To check the phyloseq table
#ps

## Alpha Diversity Commands
## Alpha Diversity of 16S rRNA data
#theme_set(theme_bw())
#plot_richness(ps, measures=c("Shannon", "Simpson"))

##Genus Level Profile
##Generate the plots: Genus, Phylum level
#plot_bar(ps, fill="Genus")
#plot_bar(ps, fill="Phylum")
setwd("~/my_shared_data_folder/bactovid/bactovid/data")
# Read metadata table
metadata1<- read.table("~/my_shared_data_folder/bactovid/bactovid/accessions/metadata_Poland.txt", sep="\t", header = TRUE)
# fonvert to dataframe
View(metadata1)
str(metadata1)

#convert data to dataframe
data1<-data.frame(data1)
View(data1)
row.names(data1)

# convert run accession index column to colum name
df1 <- cbind(run_accession = rownames(data1), data1)
View(df1)
# convert index column to intager index
rownames(df1) <- 1:nrow(df1)
View(df1)
View(metadata1)
install.packages("dplyr")
library("dplyr")
library("rlang")
#install.packages("rlang")
#concatenate metadata and ASV table based on run accesion columncomb<-full_join(USA,df,by="run_accession")
comb1<- full_join(metadata1,df1,by="run_accession",copy=TRUE)
View(comb1)
write.csv(comb1,"~/my_shared_data_folder/bactovid/bactovid/data/PolandMetadataJoinedASV2.csv")
