# Before installing DADA2, install some required dependencies
install.packages("BiocManager")
BiocManager::install("ShortRead")
BiocManager::install("Biostrings", force = TRUE)

# install dada2
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2", force = TRUE)

#load package
library(dada2)
packageVersion("dada2")
library(stats)

# 16s rRNA-seq analysis of CYD 
# If the package successfully loaded and your listed files match those here, you are ready to go through the DADA2 pipeline. 
fastq_path <- "C:/Users/Admin/Documents/DADA2_pipeline/qiime2/demux"
list.files(fastq_path)

# Single-end fastq filenames have format: SAMPLENAME.fastq.gz
fnFastq <- sort(list.files(fastq_path, pattern=".fastq.gz", full.names = TRUE))
fnFastq

# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq.gz
sample.names <- sapply(strsplit(basename(fnFastq), "_"), `[`, 1)
sample.names

##Plot Quality Control charts:visualizing the quality profiles of the forward reads
plotQualityProfile(fnFastq[1:4])

# Inspect quality profiles
quality_plots <- plotQualityProfile(fnFastq[1:4])

#Save quality plots pdf
library(ggplot2)
ggsave(plot = quality_plots, filename = "quality_reads.pdf",
       width = 10, height = 10, dpi = 300)

#Assign names of the Filtered sequences:
##Place filtered files in filtered/subdirectory
filtFiles <- file.path(fastq_path, "filtered", paste0(sample.names, "_filt.fastq.gz"))
names(filtFiles) <- sample.names

#Trimming of the sequences:
filtered_out <- filterAndTrim(fnFastq, filtFiles,
                              truncLen=c(240), trimLeft=20, maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                              compress=TRUE, verbose = TRUE)
head(filtered_out)

#Generate the Filtered Profiles and Perform Quality check and generate new profiles for filtered sequences:
quality_filt_plots = plotQualityProfile(filtFiles[1:4])

#Save quality plots pdf
library(ggplot2)
ggsave(plot = quality_filt_plots, filename = "quality_filt_reads.pdf",
       width = 10, height = 10, dpi = 300)

#Error Rate Profiles and Generate the Estimated Error Rates for Forward and Reverse Sequences
errF <- learnErrors(filtFiles)
plotErrors(errF, nominalQ=TRUE)
#To visualize the estimated error rates:Plot the error rates
error_rates <- plotErrors(errF, nominalQ=TRUE)

#Save error plots pdf
pdf("estimated_error_rates.pdf")
ggsave(plot = error_rates, filename = "estimated_error_rates.pdf",
       width = 10, height = 10, dpi = 300)
# Close the PDF device
dev.off()

#dereplication:combines all identical sequencing reads into into “unique sequences” 
#with a corresponding “abundance” equal to the number of reads with that unique sequence
derep_forward <- derepFastq(filtFiles, verbose=TRUE)
derep_forward
names(derep_forward) <- sample.names

##Sample Inference:
dadaFs <- dada(derep_forward, err=errF, multithread=TRUE, pool=TRUE)
#Inspecting the returned dada-class object:
dadaFs[[1]]

#5th step: Merging paired reads:but i have SE
#mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#mergers

#sequence table construction:
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)

# View No. of samples and No. of variants
dim(seqtab1)

# Inspect distribution of sequence lengths
table(nchar(getSequences(seqtab)))

# Then, removing the chimeras: # This step could take long if verbose=TRUE
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
View(seqtab.nochim)
# Check the no. of samples and no. of variants when chimeras are removed
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

# Track reads through the pipeline so far
getN <- function(x) sum(getUniques(x))
track <- cbind(filtered_out, getN(dadaFs), getN(derep_forward), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "dereplicated_reads", "nonchim")
rownames(track) <- sample.names

getN <- function(x) sum(getUniques(x))
track <- cbind(filtered_out, sapply(dadaFs, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
print(track)

#saving the track info in tsv file
write.table(track, "read-count-tracking.tsv", quote=FALSE, sep="\t", col.names=NA)

#Assign taxonomy
# used this at first install.packages("DECIPHER") but it didnt work so installed using bioconductor
BiocManager::install("DECIPHER")
library(DECIPHER)

IHAVE REACHED HERE!!!!!!
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
