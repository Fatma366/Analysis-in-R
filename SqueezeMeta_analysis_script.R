title: "dataset2"
author: "Stephen Ogada"
date: "7/26/2023"

# load package
library('SQMtools')

# load dataset
test = loadSQM('/home/administrator/Desktop/DERLEEN/SQUEEZEMETA/test')

# Now, there is a faster option to load your data if your project is large:
test1 = loadSQM('/home/administrator/Desktop/DERLEEN/SQUEEZEMETA/test', engine = 'data.table')

# to get abundance information for super kingdom
test$taxa$superkingdom$abund

# Matrix genera vs. abundances
genus_tax=test$taxa$genus$abund
genus_tax[1:5, 1:5]

# Matrix with COGs vs. abundances
COG_table=test$functions$COG$abund
COG_table[1:5, 1:5]

# Raw table from SQM: 13.<SQM project name>.orftable
orf_table=test$orfs$table
orf_table[1:5, 1:5]

# Then, you can always make a summary of the project to discard occasional problems.
sum_test = summary(test)
sum_test$project_name

sum_test$samples

sum_test$reads_in_orfs

# Plot taxonomy at the phylum and family levels
plotTaxonomy(test, rank='phylum', count='percent')
plotTaxonomy(test, rank='family', count='percent')

# save your plots using R functions, please check the documentation!
# Option A
#Save a svg/png file
svg('myfigure.svg') # To save a png, it would be png('myfigure.png')
plotTaxonomy(test, rank='family', count='percent') #Plot
dev.off()

# Option B
#Use ggplot2
library(ggplot2)
library(svglite)
image = plotTaxonomy(test, rank='family', count='percent')
ggsave(file='myfigure.svg', plot = image, width=10, height=8) #svg
ggsave(file='myfigure.png', plot = image, width=10, height=8) #png

# Picking up taxa without rescaling and without 'Other'
plotTaxonomy(test, rank='phylum', count='percent', tax = c('Proteobacteria','Actinobacteria','Spirochaetes'), rescale = F, color = c('yellowgreen','tan3','cornflowerblue'), others = F)

# Choosing taxa rescaling and without 'Other'
# When rescaling is turned on, the abundance values are normalized to percentages, allowing for easier comparison between different taxa. When rescaling is turned off, the abundance values are not adjusted and are displayed as they are.

plotTaxonomy(test, rank='phylum', count='percent', tax=c('Proteobacteria','Actinobacteria','Spirochaetes'), rescale = T, color = c('yellowgreen','tan3','cornflowerblue'), others = F)

# check the functional profile(KEGG, COG and PFAM annotations) of the samples
plotFunctions(test, fun_level = 'KEGG', count = 'tpm', N = 15)

# Extract the orf annotated with that KEGG
orf = test$orfs$table[Hadza$orfs$table$`KEGG ID` == 'K07133', 1:15] # see COGs, PFAM...
# Extract the sequence
sequence=test$orfs$seqs[rownames(orf)]

# We decided to pick up all contigs from 'Proteobacteria' and plot their taxonomy at the genus level and their most abundant functions
proteo=subsetTax(test, 'phylum',tax='Proteobacteria', rescale_copy_number = F)
plotTaxonomy(proteo, 'genus','percent', N=10, rescale = T, others = T)
plotFunctions(proteo, fun_level = 'KEGG',count = 'copy_number')

# subset all the genes whose functional annotations contain the word 'antibiotic
antibiotic = subsetFun(test, fun = 'antibiotic', rescale_copy_number = F) 
plotTaxonomy(antibiotic, 'genus','percent', N = 10, rescale = T, others = T)
plotFunctions(antibiotic, fun_level = 'KEGG',count = 'tpm')

# you can subset all the genes related to a particular KEGG pathway
aromatic_aa = subsetFun(test, fun = 'Phenylalanine, tyrosine and tryptophan biosynthesis', rescale_tpm = F,rescale_copy_number = F)

# useful to analyze which taxa are involved in a specific KEGG pathway
# Plot taxonomy
plotTaxonomy(aromatic_aa, 'genus', 'percent', N = 10, rescale = F, others = T)
plotFunctions(aromatic_aa, fun_level = 'KEGG', count = 'tpm')

# The function subsetFun also accepts regular expressions if the parameter fixed = F
# Select genes annotated with K06001 OR K01696
trp_synthesis = subsetFun(test, fun = 'K06001|K01696', fixed = F)
plotTaxonomy(trp_synthesis, rank = 'genus')

# if you are interested in bins, you can subset them from the SQM object
maxbin005 = subsetBins(test, 'maxbin.005.fasta.contigs' )
# Plot functions
plotFunctions(maxbin005, fun_level = 'KEGG', count = 'tpm')

# Combine different subsets
bacteroidetes = subsetTax(test, 'phylum','Bacteroidetes')
bacteroidetes_carbohydrates = subsetFun(bacteroidetes, 'Carbohydrate metabolism')
proteobacterias = subsetTax(test, 'phylum','Proteobacteria')
proteobacterias_aminoacids = subsetFun(proteobacterias, 'Amino acid metabolism')
bact_carbo_proteo_amino = combineSQM(bacteroidetes_carbohydrates, proteobacterias_aminoacids, rescale_copy_number = F, tax_source = "contigs")
bact_carbo_proteo_amino = combineSQM(bacteroidetes_carbohydrates, proteobacterias_aminoacids, rescale_copy_number = F)

plotTaxonomy(bact_carbo_proteo_amino, 'phylum','percent', N = 10, rescale = T, others = T)

plotFunctions(bact_carbo_proteo_amino, fun_level = 'KEGG', count = 'tpm')

# Export table
exportTable(bact_carbo_proteo_amino$functions$COG$tpm, 'test_COG.tsv')  
exportKrona(test)

# Export Pathways
# To know if a KEGG is present or not in one sample,
# you can provide a vector of colors.
# Colors must be in the same order that samples
# To distinguish between Hadza and Italians:
# Repeat one color 8 times because the first 8 samples are from Hadza.
# The next 8 samples are from Italians, repeat their color.
colors = c(rep('#006682', 8), rep('#c26e00', 8))
exportPathway(test, '00400', output_suffix = 'aromatic_aa', sample_colors = colors, max_scale_value = 3, count= 'copy_number')

# To calculate the log-fold change of a KEGG, you can provide a list of vectors.
# One vector with the samples names by condition:
# test samples
H.samples = test$misc$samples[grep('H', test$misc$samples)]
IT.samples = test$misc$samples[grep('IT', test$misc$samples)]
condition = list(H.samples, IT.samples)
# Choose one color per condition
colors = c('#006682', '#c26e00')
# Plot log2 fold changes using copy number abundances in the selected KEGG pathway:
exportPathway(test, '00400', output_suffix = 'aromatic_aa.log2FC', fold_change_colors = colors, fold_change_groups = condition, count ='copy_number', max_scale_value = 1.5)

# Differencial abundances with DESeq2
# Tutorial based on http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# Load the package
library('DESeq2')
# Required: matrix with the raw abundances, colData (metadata file)
metadata=as.data.frame(c(rep('H',8), rep('I',8)))
rownames(metadata)=colnames(test$functions$KEGG$abund)
colnames(metadata)='condition'
# Verify sample order
all(rownames(metadata)==colnames(test$functions$KEGG$abund))
# Convert your data to the format required by DESeq2
dds = DESeqDataSetFromMatrix(countData = test$functions$KEGG$abund, colData = metadata, design = ~ condition)
# Remove low abundant KEGGs:
keep = rowSums(counts(dds)) >= 10
dds = dds[keep,]
# Choose factor levels
dds$condition=factor(dds$condition, levels=c('H','I'))
# Run DESeq2
dds2=DESeq(dds)
results=results(dds2, name='condition_I_vs_H')
plotMA(results, ylim=c(-2,2))

#  Use the vegan package with the SQM Object
library('vegan')
metadata = as.data.frame(c(rep('test', 8), rep('Italy', 8)))
rownames(metadata) = colnames(test$functions$KEGG$tpm)
colnames(metadata) = 'condition'
# Tranpose the matrix to have samples in rows.
kegg_tpm = t(test$functions$KEGG$tpm)
MDS = metaMDS(kegg_tpm)
colvec = c('#006682','#c26e00') # colors
plot(MDS, display = 'sites')
with(metadata, points(MDS, display = 'sites', col = colvec[condition], pch = 19))
with(metadata, legend('topright', legend = levels(condition), col = colvec, pch = 21, pt.bg = colvec))

# Additionally, you can run different statistical tests (eg: PERMANOVA test) to know if some variable (like the location) is affecting the clustering of samples in the ordination
library('vegan')
metadata = as.data.frame(c(rep('test', 8), rep('Italy', 8)))
rownames(metadata) = colnames(test$functions$KEGG$tpm)
colnames(metadata) = 'condition'
# Tranpose the matrix to have samples in rows.
kegg_tpm = t(test$functions$KEGG$tpm)
country = c(rep('test', 8), rep('Italy', 8))
adonis(kegg_tpm~country)
Call :
  adonis (formula = kegg_tpm ~ country)
Permutation : free

# Save your R workspace
save.image('mySession.RData') #Save  R workspace
load('mySession.RData') #Load a previous  R workspace

