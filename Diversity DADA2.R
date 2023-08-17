# https://astrobiomike.github.io/amplicon/dada2_workflow_ex#assigning-taxonomy
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
library(ggplot2)
library(phyloseq)


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
ggsave(plot = quality_plots, filename = "quality_reads.pdf",
       width = 10, height = 10, dpi = 300)

#Assign names of the Filtered sequences:
##Place filtered files in filtered/subdirectory
filtFiles <- file.path(fastq_path, "filtered", paste0(sample.names, "_filt.fastq.gz"))
names(filtFiles) <- sample.names

#Trimming of the sequences:
filtered_out <- filterAndTrim(fnFastq, filtFiles,
                              truncLen=c(240), maxEE=2, truncQ=2, maxN=0, rm.phix=TRUE,
                              compress=TRUE, verbose = TRUE)
head(filtered_out)
print(filtered_out)
#Generate the Filtered Profiles and Perform Quality check and generate new profiles for filtered sequences:
quality_filt_plots <- plotQualityProfile(filtFiles[1:4])

conflicts()
detach("package:dplyr", unload = TRUE)

#Save quality plots pdf
library(ggplot2)
ggsave(plot = quality_filt_plots, filename = "quality_filt_reads.pdf",
       width = 10, height = 10, dpi = 300)

#Error Rate Profiles and Generate the Estimated Error Rates for Forward and Reverse Sequences
errF <- learnErrors(filtFiles)
plotErrors(errF, nominalQ=TRUE)
#To visualize the estimated error rates:Plot the error rates
error_rates <- plotErrors(errF, nominalQ=TRUE)

#Save error plots png
ggsave(plot = error_rates, filename = "estimated_error_rates.png",
       width = 10, height = 10, dpi = 300)

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

#Create a DNAStringSet from the ASVs
dna <- DNAStringSet(getSequences(seqtab.nochim))
url="http://www2.decipher.codes/Classification/TrainingSets/SILVA_SSU_r138_2019.binder.RData"
destfile="SILVA_SSU_r138_2019.binder"
# download.file(url, destfile)
load("SILVA_SSU_r138_2019.binder.RData")
id_info <- IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE)

#Giving the seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

#Make and write out a fasta of our final ASV seqs:
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "ASVs.fa")

#Count table: genrates asv tables with count and the sequences
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts.tsv", sep="\t", quote=F, col.names=NA)

#Convert the output object of class "Taxa" to a matrix analogous to the output from assignTaxonomy
#Tax table:
#Create table of taxonomy and setting any that are unclassified as "NA"
ranks <- c("domain", "phylum", "class", "order", "family", "genus", "species")
asv_tax <- t(sapply(id_info, function(x) {
  m <- match(ranks, x$rank)
  taxa <- x$taxon[m]
  taxa[startsWith(taxa, "unclassified_")] <- NA
  taxa
}))

colnames(asv_tax) <- ranks
rownames(asv_tax) <- gsub(pattern=">", replacement="", x=asv_headers)

write.table(asv_tax, "ASVs_taxonomy.tsv", sep = "\t", quote=F, col.names=NA)

# remove contaminants using decontam
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

BiocManager::install("decontam")
library(decontam)
packageVersion("decontam")
vector_for_decontam <- c(rep(TRUE), rep(FALSE))
contam_df <- isContaminant(t(asv_tab), neg=vector_for_decontam)

#Getting vector holding the identified contaminant IDs
contam_asvs <- row.names(contam_df[contam_df$contaminant == TRUE, ])

#Make new fasta file
contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
dont_want <- sort(c(contam_indices, contam_indices + 1))
asv_fasta_no_contam <- asv_fasta[- dont_want]

#Make new count table
asv_tab_no_contam <- asv_tab[!row.names(asv_tab) %in% contam_asvs, ]

#Make new taxonomy table
asv_tax_no_contam <- asv_tax[!row.names(asv_tax) %in% contam_asvs, ]

#Write them out to files
write(asv_fasta_no_contam, "ASVs-no-contam.fa")
write.table(asv_tab_no_contam, "ASVs_counts-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)
write.table(asv_tax_no_contam, "ASVs_taxonomy-no-contam.tsv",
            sep="\t", quote=F, col.names=NA)

#Reading the data in R
count_tab <- read.table("ASVs_counts-no-contam.tsv", header=T, row.names=1,
                        check.names=F, sep="\t")
tax_tab <- as.matrix(read.table("ASVs_taxonomy-no-contam.tsv", header=T,
                                row.names=1, check.names=F, sep="\t"))

sample_info_tab <- read.table("C:/Users/Admin/Documents/DADA2_pipeline/qiime2/sample-metadata.tsv", header=F, row.names=1,
                              check.names=F, sep="\t")
# Provide explicit column names
col_names <- c("Counties", "BarcodeSequence", "LinkerPrimerSequence", "BarcodeName", "ReversePrimer","ProjectName", "Description")  # Replace with actual column names
colnames(sample_info_tab) <- col_names
View(sample_info_tab)

BiocManager::install("phyloseq")
library(phyloseq)

#Make phyloseq object with transformed table
count_tab_phy <- otu_table(count_tab, taxa_are_rows=T)
tax_tab_phy <- tax_table(tax_tab)
sample_info_tab_phy <- sample_data(sample_info_tab)
sample_info_tab$Counties <- as.character(sample_info_tab$Counties)
View(sample_info_tab)

# first we need to create a phyloseq object using our un-transformed count table
ASV_physeq <- phyloseq(count_tab_phy, tax_tab_phy, sample_info_tab_phy)

#Plot richness using the plot_richness() function on the phyloseq object
alpha_plot <- plot_richness(ASV_physeq, measures=c("Shannon", "Simpson"), color = "Counties")
alpha_plot

ggsave(plot = alpha_plot, filename = "Alpha_diversity_plot.png",
       width = 10, height = 10, dpi = 300)
#plot2
alpha_plot2 <- plot_richness(ASV_physeq,x="Counties", color="Description", measures=c("Shannon", "Simpson"))
alpha_plot2

ggsave(plot = alpha_plot, filename = "Alpha_diversity_plot2.png",
       width = 10, height = 10, dpi = 300)


#Box_plot1
alpha_box_plot <- plot_richness(ASV_physeq, color="Counties", measures=c("Observed", "Shannon")) + geom_boxplot()
alpha_box_plot

ggsave(plot = alpha_box_plot, filename = "Alpha_diversity_boxplot.png",
       width = 10, height = 10, dpi = 300)

#Box_plot2
alpha_box_plot2 <- plot_richness(ASV_physeq,  x="Counties", color="Description", measures=c("Observed", "Shannon")) + geom_boxplot()
alpha_box_plot2

ggsave(plot = alpha_box_plot2, filename = "Alpha_diversity_boxplot2.png",
       width = 10, height = 10, dpi = 300)
# i did not like the box plot here, i prefer the dot plots, better to explain

#Beta diversity
#PCOA plot using Bray-Curtis as distance
bray_curtis_dist = phyloseq::distance(ASV_physeq, method="bray", weighted=F)
ord = ordinate(ASV_physeq, method="PCoA", distance=bray_curtis_dist)
Bray_plt <- plot_ordination(ASV_physeq, ord, color="Counties") + theme(aspect.ratio=1)
Bray_plt

BiocManager::install("DESeq2")
library(DESeq2)

# first we need to make a DESeq2 object
deseq_counts <- DESeqDataSetFromMatrix(count_tab, colData = sample_info_tab, design = ~Counties) 
# we have to include the "colData" and "design" arguments because they are 
# required, as they are needed for further downstream processing by DESeq2, 
# but for our purposes of simply transforming the data right now, they don't 
# matter

deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)
# NOTE: If you get this error here with your dataset: "Error in
# estimateSizeFactorsForMatrix(counts(object), locfunc =locfunc, : every
# gene contains at least one zero, cannot compute log geometric means", that
# can be because the count table is sparse with many zeroes, which is common
# with marker-gene surveys. In that case you'd need to use a specific
# function first that is equipped to deal with that. You could run:
# deseq_counts <- estimateSizeFactors(deseq_counts, type = "poscounts")
# now followed by the transformation function:
# deseq_counts_vst <- varianceStabilizingTransformation(deseq_counts)

# and here is pulling out our transformed table
vst_trans_count_tab <- assay(deseq_counts_vst)

# and calculating our Euclidean distance matrix
euc_dist <- dist(t(vst_trans_count_tab))

# Now that we have our Euclidean distance matrix, let’s make and plot 
# a hierarchical clustering of our samples
euc_clust <- hclust(euc_dist, method="ward.D2")

# hclust objects like this can be plotted with the generic plot() function
hclust = plot(euc_clust) 
# but i like to change them to dendrograms for two reasons:
# 1) it's easier to color the dendrogram plot by groups
# 2) if wanted you can rotate clusters with the rotate() 
#    function of the dendextend package
pdf("Hierachical_clustering1.pdf")
ggsave(plot = hclust, filename = "Hierachical_clustering1.pdf",
       width = 10, height = 10, dpi = "retina")

euc_dend <- as.dendrogram(euc_clust, hang=0.1)

# making our phyloseq object with transformed table
vst_count_phy <- otu_table(vst_trans_count_tab, taxa_are_rows=T)
sample_info_tab_phy <- sample_data(sample_info_tab)
vst_physeq <- phyloseq(vst_count_phy, sample_info_tab_phy)

# generating and visualizing the PCoA with phyloseq
vst_pcoa <- ordinate(vst_physeq, method="MDS", distance="euclidean")
eigen_vals <- vst_pcoa$values$Eigenvalues # allows us to scale 

# Define a color palette for counties
county_colors <- c("Kwale" = "red", "Kilifi" = "blue", "Lamu" = "green")

# Plot the ordination
Beta_PCOA <- plot_ordination(vst_physeq, vst_pcoa, color = "Counties") + 
  geom_point(size = 1) + labs(col = "Counties") + 
  geom_text(aes(label = rownames(sample_info_tab), hjust = 0.3, vjust = -0.4)) + 
  coord_fixed(sqrt(eigen_vals[2] / eigen_vals[1])) + 
  ggtitle("PCoA") + 
  scale_color_manual(values = county_colors) + 
  theme_bw() + theme(legend.position = "none")

ggsave(plot = Beta_PCOA, filename = "Beta_PCOA_plot.png",
       width = 10, height = 10, dpi = 300)

#Abundance bar plot for phylum
sorting <- names(sort(taxa_sums(ASV_physeq), decreasing=TRUE))
ASV_physeq.sorting <- transform_sample_counts(ASV_physeq, function(OTU) OTU/sum(OTU))
ASV_physeq.sorting <- prune_taxa(sorting, ASV_physeq.sorting)
abundance2 <- plot_bar(ASV_physeq.sorting, fill="phylum") + facet_wrap(~Counties, scales="free_x")

ggsave(plot = abundance2, filename = "phylum_bar_plot.png",
       width = 15, height = 10, dpi = 300)

#Abundance bar plot for genus
sorting <- names(sort(taxa_sums(ASV_physeq), decreasing=TRUE))
ASV_physeq.sorting <- transform_sample_counts(ASV_physeq, function(OTU) OTU/sum(OTU))
ASV_physeq.sorting <- prune_taxa(sorting, ASV_physeq.sorting)
abundance <- plot_bar(ASV_physeq.sorting, fill="genus") + facet_wrap(~Counties, scales="free_x")

ggsave(plot = abundance, filename = "Classification_bar_plot.png",
       width = 15, height = 10, dpi = 300)

# if i want to plot using ggplot2
library(ggplot2)

# Convert the phyloseq object to a data frame
abundance_data <- as.data.frame(otu_table(ASV_physeq.sorting))
print(abundance_data)

# Plot the taxonomic bar graph using ggplot2
# need to look into this command still not working
ggplot(abundance_data, aes(x = Sample, y = Abundance, fill = "Phylum")) +
  geom_bar(stat = "identity") +
  facet_wrap(~ Counties, scales = "free_x") +
  labs(x = "Samples", y = "Abundance") +
  theme_minimal()

ggsave(plot = abundance_data, filename = "Tax_Phylum_bar_plot.png",
       width = 15, height = 10, dpi = 300)

#In order to train the machine learning models, we need to have labels of patients
#and controls for each sample (counties). This is why we need the metadata including three 
#importing columns: same as the one on the ASV table and 
#the column indicating counties. After you have that 
#metadata, this  part of the script does this. 
#The output of this script will be the input of the notebook named 
#preprocessing_machine_learning.ipy

install.packages("dplyr")
install.packages("rlang")
library("dplyr")
library("rlang")

# Concatenate metadata and ASV table based on column "Sample"
comb <- full_join(sample_info_tab, abundance_data, by = "Sample", copy = TRUE)

# Pivot the data into long format
data_long <- comb %>%
  pivot_longer(cols = starts_with("ASV_"), names_to = "ASV", values_to = "Abundance")

# View the combined and pivoted data
View(data_long)

# Write the combined and pivoted data to a CSV file
write.csv(data_long, "path/to/your/output/file.csv")

