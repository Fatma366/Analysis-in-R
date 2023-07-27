library('SQMtools')
Hadza = loadSQM('/media/disk6/natalia/Hadza16')
# Now, there is a faster option to load your data if your project is large:
Hadza = loadSQM('/media/disk6/natalia/Hadza16', engine = 'data.table')

# Matrix genera vs. abundances
genus_tax=Hadza$taxa$genus$abund
genus_tax[1:5, 1:5]

# Matrix with COGs vs. abundances
COG_table=Hadza$functions$COG$abund
COG_table[1:5, 1:5]

# Raw table from SQM: 13.<SQM project name>.orftable
orf_table=Hadza$orfs$table
orf_table[1:5, 1:5]

# Then, you can always make a summary of the project to discard occasional problems.
sum_Hadza = summary(Hadza)
sum_Hadza$project_name
[1] "Hadza16"
sum_Hadza$samples
[1] "H1"   "H10"  "H11"  "H12"  "H13"  "H15"  "H14"  "H16"  "IT1"  "IT11"
[11] "IT13" "IT14" "IT2"  "IT3"  "IT4"  "IT5"
sum_Hadza$reads_in_orfs
H1      H10      H11      H12      H13      H15      H14      H16      IT1     IT11     IT13     IT14      IT2      IT3      IT4      IT5 
23847271 52692441 13718056  7833832 23807939 16412877 61547031 25641150  3191223  7090403  3415890  2472831 15259976 12563001  4362936 29855475

# Plot taxonomy at the phylum and family levels
plotTaxonomy(Hadza, rank='phylum', count='percent')
plotTaxonomy(Hadza, rank='family', count='percent')

# save your plots using R functions, please check the documentation!
# Option A
#Save a svg/png file
svg('myfigure.svg') # To save a png, it would be png('myfigure.png')
plotTaxonomy(Hadza, rank='family', count='percent') #Plot
dev.off()

# Option B
#Use ggplot2
library(ggplot2)
image = plotTaxonomy(Hadza, rank='family', count='percent')
ggsave(file='myfigure.svg', plot = image, width=10, height=8) #svg
ggsave(file='myfigure.png', plot = image, width=10, height=8) #png

# Picking up taxa without rescaling and without 'Other'
plotTaxonomy(Hadza, rank='phylum', count='percent', tax = c('Proteobacteria','Actinobacteria','Spirochaetes'), rescale = F, color = c('yellowgreen','tan3','cornflowerblue'), others = F)

# Choosing taxa rescaling and without 'Other'
plotTaxonomy(Hadza, rank='phylum', count='percent', tax=c('Proteobacteria','Actinobacteria','Spirochaetes'), rescale = T, color = c('yellowgreen','tan3','cornflowerblue'), others = F)

# check the functional profile(KEGG, COG and PFAM annotations) of the samples
plotFunctions(Hadza, fun_level = 'KEGG', count = 'tpm', N = 15)

# Extract the orf annotated with that KEGG
orf = Hadza$orfs$table[Hadza$orfs$table$`KEGG ID` == 'K07133', 1:15] # see COGs, PFAM...
# Extract the sequence
sequence=Hadza$orfs$seqs[rownames(orf)]

megahit_1163508_371-1111 
"MVEKSNIDMLNRKIYSYLRDFFETEKKALLVSGARQVGKTFAIRKVGAECFADVVEFNFLNNPKYREAFKSPSDAKEILLRLSALAEKKLIPGTTLVFFDEVQECPEMVTAIKFLVEEGSYRYVMSGSLLGVELKDIRSVPVGYMAERDENYKVEGYEKGMQRFYAYLSLFTMSMLGLVVATNIFQMYLFWELVGVCSYLLIGFYYPKHAAVHASKKAFKKLELDPPSHLATRCAMGTSHSYRHHR*" 

# We decided to pick up all contigs from 'Proteobacteria' and plot their taxonomy at the genus level and their most abundant functions
proteo=subsetTax(Hadza, 'phylum',tax='Proteobacteria', rescale_copy_number = F)
plotTaxonomy(proteo, 'genus','percent', N=10, rescale = T, others = T)
plotFunctions(proteo, fun_level = 'KEGG',count = 'copy_number')

# subset all the genes whose functional annotations contain the word 'antibiotic
antibiotic = subsetFun(Hadza, fun = 'antibiotic', rescale_copy_number = F) 
plotTaxonomy(antibiotic, 'genus','percent', N = 10, rescale = T, others = T)
plotFunctions(antibiotic, fun_level = 'KEGG',count = 'tpm')

# you can subset all the genes related to a particular KEGG pathway
aromatic_aa = subsetFun(Hadza, fun = 'Phenylalanine, tyrosine and tryptophan biosynthesis', rescale_tpm = F,rescale_copy_number = F)

# useful to analyze which taxa are involved in a specific KEGG pathway
# Plot taxonomy
plotTaxonomy(aromatic_aa, 'genus', 'percent', N = 10, rescale = F, others = T)
plotFunctions(aromatic_aa, fun_level = 'KEGG', count = 'tpm')

# The function subsetFun also accepts regular expressions if the parameter fixed = F
# Select genes annotated with K06001 OR K01696
trp_synthesis = subsetFun(Hadza, fun = 'K06001|K01696', fixed = F)
plotTaxonomy(trp_synthesis, rank = 'genus')

# if you are interested in bins, you can subset them from the SQM object
maxbin005 = subsetBins(Hadza, 'maxbin.005.fasta.contigs' )
# Plot functions
plotFunctions(maxbin005, fun_level = 'KEGG', count = 'tpm')

# Combine different subsets
bacteroidetes = subsetTax(Hadza, 'phylum','Bacteroidetes')
bacteroidetes_carbohydrates = subsetFun(bacteroidetes, 'Carbohydrate metabolism')
proteobacterias = subsetTax(Hadza, 'phylum','Proteobacteria')
proteobacterias_aminoacids = subsetFun(proteobacterias, 'Amino acid metabolism')
bact_carbo_proteo_amino = combineSQM(bacteroidetes_carbohydrates, proteobacterias_aminoacids, rescale_copy_number = F, tax_source = "contigs")
bact_carbo_proteo_amino = combineSQM(bacteroidetes_carbohydrates, proteobacterias_aminoacids, rescale_copy_number = F)

plotTaxonomy(bact_carbo_proteo_amino, 'phylum','percent', N = 10, rescale = T, others = T)

plotFunctions(bact_carbo_proteo_amino, fun_level = 'KEGG', count = 'tpm')

# Export table
exportTable(bact_carbo_proteo_amino$functions$COG$tpm, 'Hadza_COG.tsv')  
exportKrona(Hadza)

# Export Pathways
# To know if a KEGG is present or not in one sample,
# you can provide a vector of colors.
# Colors must be in the same order that samples
# To distinguish between Hadza and Italians:
# Repeat one color 8 times because the first 8 samples are from Hadza.
# The next 8 samples are from Italians, repeat their color.
colors = c(rep('#006682', 8), rep('#c26e00', 8))
exportPathway(Hadza, '00400', output_suffix = 'aromatic_aa', sample_colors = colors, max_scale_value = 3, count= 'copy_number')

# To calculate the log-fold change of a KEGG, you can provide a list of vectors.
# One vector with the samples names by condition:
#Hadza samples
H.samples = Hadza$misc$samples[grep('H', Hadza$misc$samples)]
IT.samples = Hadza$misc$samples[grep('IT', Hadza$misc$samples)]
condition = list(H.samples, IT.samples)
# Choose one color per condition
colors = c('#006682', '#c26e00')
# Plot log2 fold changes using copy number abundances in the selected KEGG pathway:
exportPathway(Hadza, '00400', output_suffix = 'aromatic_aa.log2FC', fold_change_colors = colors, fold_change_groups = condition, count ='copy_number', max_scale_value = 1.5)

# Differencial abundances with DESeq2
# Tutorial based on http://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
# Load the package
library('DESeq2')
# Required: matrix with the raw abundances, colData (metadata file)
metadata=as.data.frame(c(rep('H',8), rep('I',8)))
rownames(metadata)=colnames(Hadza$functions$KEGG$abund)
colnames(metadata)='condition'
# Verify sample order
all(rownames(metadata)==colnames(Hadza$functions$KEGG$abund))
# Convert your data to the format required by DESeq2
dds = DESeqDataSetFromMatrix(countData = Hadza$functions$KEGG$abund, colData = metadata, design = ~ condition)
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
metadata = as.data.frame(c(rep('Hadza', 8), rep('Italy', 8)))
rownames(metadata) = colnames(Hadza$functions$KEGG$tpm)
colnames(metadata) = 'condition'
# Tranpose the matrix to have samples in rows.
kegg_tpm = t(Hadza$functions$KEGG$tpm)
MDS = metaMDS(kegg_tpm)
colvec = c('#006682','#c26e00') # colors
plot(MDS, display = 'sites')
with(metadata, points(MDS, display = 'sites', col = colvec[condition], pch = 19))
with(metadata, legend('topright', legend = levels(condition), col = colvec, pch = 21, pt.bg = colvec))

# Additionally, you can run different statistical tests (eg: PERMANOVA test) to know if some variable (like the location) is affecting the clustering of samples in the ordination
library('vegan')
metadata = as.data.frame(c(rep('Hadza', 8), rep('Italy', 8)))
rownames(metadata) = colnames(Hadza$functions$KEGG$tpm)
colnames(metadata) = 'condition'
# Tranpose the matrix to have samples in rows.
kegg_tpm = t(Hadza$functions$KEGG$tpm)
country = c(rep('Hadza', 8), rep('Italy', 8))
adonis(kegg_tpm~country)
Call :
  adonis (formula = kegg_tpm ~ country)
Permutation : free
Number of permutation: 999
Terms added sequentially (first to last)
Df SumsOfSqs  MeanSqs F.Model  R2      Pr(>F)
country       1  0.034942   0.034942 7.9894 0.36333 0.001 ***
  Residuals     14 0.061231   0.004374        0.63667
Total         15 0.096173                   1.00000
---
  Signif. codes:
  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

# Save your R workspace
save.image('mySession.RData') #Save  R workspace
load('mySession.RData') #Load a previous  R workspace

