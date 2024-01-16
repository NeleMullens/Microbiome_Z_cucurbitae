# Start data analysis ----------------------------

--------------------------------------------------------------------------------
##                                                                            ##
##                          Processing of data                                ##
##                                                                            ##
--------------------------------------------------------------------------------
# First step is quality control of data using FastQC (https://github.com/s-andrews/FastQC)
# This is executed on linux terminal (Bash)using code:
# fastqc filename --outdir output_directory
  
#phyloseq and dada2
if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")}
BiocManager::install("phyloseq")
BiocManager::install("dada2", version = "3.12")

#Load packages and dependencies
library(dada2)
library(phyloseq)
library(ggplot2)


## Begin DADA2 pipeline -------------------------------------------
#DADA2 For amplicon sequencing of V3-V4 16S rRNA genes (Bacteria) Instead of
#working with Uparse for OTU clustering, we will make use of DADA2 for
#identifying sequence variants. DADA2 will infer exact sequence variants from
#the data.


#First we read in the names of the fastq files, and perform some string
#manipulation to get lists of the forward and reverse fastq files in matched
#order

#Sort the reads (Forward and Reverse) so they are in the same order
fnFs <- sort(list.files("data", pattern="_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files("data", pattern="_2.fastq.gz", full.names = TRUE))

#check quality of reads with `rplotQualityProfile`
plotQualityProfile(fnFs[1:6])
plotQualityProfile(fnRs[1:6])


### Filtering and trimming --------------------------------- 
#Extract sample names (take subset ('[') from basename, when the basename string
#is split at "_")
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
#Alternatively, use purr: 
##unlist(purrr::map(strsplit(basename(fnFs),"_"),purrr::pluck,1))

#Prepare filenames/path for filtered data
filtFs <- file.path("data/Filtered_data", paste0(sample.names, "_F_filtered.fq.gz"))
filtRs <- file.path("data/Filtered_data", paste0(sample.names, "_R_filtered.fq.gz"))
#necessary for the dereplication and merging step 
names(filtFs) <- sample.names 
names(filtRs) <- sample.names

#"For the filtering, we’ll use standard filtering parameters: maxN=0 (DADA2
#requires no Ns), truncQ=2 and rm.phix=TRUE. maxEE will be set at 3, whereas
#standard this is set on 2. Here we will use a value of 3 to better compare with
#the Uparse pipeline (where a value of 3 is chosen as a standard). The maxEE
#parameter sets the maximum number of “expected errors” allowed in a read, which
#is a better filter than simply averaging quality scores.

#trimLeft is used to remove primers. Don't use external trimming programmes, as
#they remove primer mismatches, but due to the high number of unknown bacteria,
#you will have these in your dataset and you would thus remove actual biological
#variation

#r filtering_trimming #runtime of 2-3 hours?
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=c(17, 21), truncLen=c(260,240),
                     maxN=0, maxEE=c(1,1), truncQ=2, rm.phix=TRUE,
                     compress=TRUE) # On non-Windows set multithread=TRUE
head(out)
dim(out)
plotQualityProfile(filtFs[1:6])
plotQualityProfile(filtRs[1:6])

### Estimate the error rate --------------------------------- 

#The DADA2 algorithm depends on a parametric error model (err) and every
#amplicon dataset has a different set of error rates. The learnErrors method
#learns the error model from the data, by alternating estimation of the error
#rates and inference of sample composition until they converge on a jointly
#consistent solution. As in many optimization problems, the algorithm must begin
#with an initial guess, for which the maximum possible error rates in this data
#are used (the error rates if only the most abundant sequence is correct and all
#the rest are errors).

#runtime around 45 minutes for each one
#>**Note** Parameter learning is computationally intensive, so by default the
#*learnErrors* function uses only a subset of the data (the first 1M reads). If
#the plotted error model does not look like a good fit, try increasing the
#nreads parameter to see if the fit improves.
errF <- learnErrors(filtFs, nread=1e6, multithread=TRUE)
errR <- learnErrors(filtRs, nread=1e6, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

### Dereplication  ---------------------------------
#Dereplication optinal?
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names
names(derepRs) <- sample.names


### Sample Inference  --------------------------------- 
#With dereplication step done
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#If dereplication is skipped #without dereplication compute time still okay
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

#Inspect the returned dada-class object
dadaFs[[1]]


### Merge paired reads ---------------------------------
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])


### Construct sequence table --------------------------------------------------
#make a sequence table of the merged reads with the counts per sample
seqtab <- makeSequenceTable(mergers)

#Check table
dim(seqtab)
sum(seqtab)

#Inspect distribution of sequence length
table(nchar(getSequences(seqtab)))

#save the sequence table as an RDS file, change the name according to the run
saveRDS(seqtab, "outputs/seqtab_Microbioom_all.rds")


### Remove chimeras -----------------------------------------------------------

#this step takes around 1 hour to run
#seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
seqtab.nochim2 <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE, minFoldParentOverAbundance=2)
dim(seqtab.nochim2)
sum(seqtab.nochim2)

#table(nchar(getSequences(seqtab.nochim)))
table(nchar(getSequences(seqtab.nochim2)))

sum(seqtab.nochim2)/sum(seqtab)

saveRDS(seqtab.nochim2, "outputs/seqtab.nochim_All.rds")


### Taxonomy assignment ---------------------------------------------------------

#Download training dataset. Check if download succesful:
tools::md5sum("data/silva_nr99_v138.1_train_set.fa.gz")

#Make a training dataset. Runs for around 5 minutes
taxa <- assignTaxonomy(seqtab.nochim2, "data/silva_nr99_v138.1_train_set.fa.gz", multithread=TRUE)

#Show taxa in dataset
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
head(taxa.print)

table<-as.data.frame(t(seqtab.nochim2))
write.table(table,file="outputs/Table_with_reads_and_taxonomic_data.txt",sep="\t",col.names=TRUE,row.names=TRUE)

table<-as.data.frame(t(taxa))
write.table(table,file="outputs/Table_with_bacteria.txt",sep="\t",col.names=TRUE,row.names=TRUE)

## End of DADA2 pipeline --------------------------------------------------

## Combine taxonomy and data in one file -------------------------------
dim(seqtab.nochim2)
dim(taxa)

#transpose table
table<-as.data.frame(t(seqtab.nochim2))

#remove rownames
rownames(table)<-NULL

#add extra columns containing taxa, for this, first transform the taxa object to a data frame
taxadf<-as.data.frame(taxa)
table$Kingdom<-taxadf$Kingdom
table$Phylum<-taxadf$Phylum
table$Class<-taxadf$Class
table$Order<-taxadf$Order
table$Family<-taxadf$Family
table$Genus<-taxadf$Genus

#add extra column containing sequences
table$seq<-as.vector(colnames(seqtab.nochim2))

#write table to output file
write.table(table,"outputs/Full_table_with_taxonmy.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)


### Remove contaminations from wet lab procedure --------------------------------

# Go to https://github.com/donaldtmcknight/microDecon
install.packages("devtools")
install.packages("cli")
devtools::install_github("donaldtmcknight/microDecon") #Installs microDecon
library(microDecon)

# transform Full_table_with_taxonmy.txt: add a first column with OTUs, and the negative control(s) should be next. Remove taxonomic data (or reduce to one column, at the end)
clean_table <- read.csv("outputs/Clean_blanks_table.csv")
clean_tab <- as.data.frame(clean_table)
View(clean_tab)
# if error, check if there are empty columns and remove these through clean_tab$X <- NULL
contamdone<-decon(data=clean_tab, numb.blanks=1, numb.ind=c(4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4), taxa=T) #all 96 samples in groups of 4; based on location, fruit, and treatment
contamdone2<-contamdone$decon.table
contamdone2  #with this step we can see which ASV have been removed due to too low counts in samples compared to the negative control and how many counts have been subtracted for the ASVs that have not been removed.

cleaned<-as.data.frame(contamdone2)
write.table(cleaned,file="outputs/Filtered_OTUs.txt",sep="\t",col.names=TRUE,row.names=TRUE)

#Rows and columns are swithed, so switch them around:
library(tibble)
switched <- as_tibble(t(cleaned[,-1]))
write.table(switched,file="outputs/Filtered_OTUs.txt",sep="\t",col.names=TRUE,row.names=TRUE)

#### manually filter out Mitochondria and Chloroplast from file 
#File: outputs/Complete_filtered_data.csv
#File without negative control: Complete_filtered_data_neg_removed.csv


## Create a frequency table --------------------------------
table2 <- read.delim("outputs/Complete_filtered_data.txt")

View(table2)
#The last seven columns are for taxonomic identification, we want to select all but that
samplenumber<-ncol(table2)-7
#The samples are all columbs starting from 1 to the number calculated above
samples<-colnames(table2)[1:samplenumber]

table2 <- read.csv("outputs/Complete_freq.csv", row.names = 1)
samplenumber<-ncol(table2)-7

#Place the samples in a dataframe
freq_table <- data.frame(matrix(nrow = nrow(table2), ncol=samplenumber, 
                                dimnames = list(c(rownames(table2)), 
                                                c(colnames(table2[1:samplenumber])))))

#calculate frequencies
for (i in 1:samplenumber){
  freq_table[,i] <- table2[,i]/colSums(table2[1:samplenumber])[i]
}
# [,i] => select all elements from column i

write.table(freq_table, "outputs/freq_table.txt", sep="\t")

#save(freq_table, "outputs/freq_table.RData" )

#add species information to the frequency table again
freq_table$Kingdom<-taxadf$Kingdom
freq_table$Phylum<-taxadf$Phylum
freq_table$Class<-taxadf$Class
freq_table$Order<-taxadf$Order
freq_table$Family<-taxadf$Family
freq_table$Genus<-taxadf$Genus

#write to new output table
write.table(freq_table,"outputs/frequency_table_complete.txt",sep="\t",col.names=TRUE,row.names=FALSE,quote=FALSE)



  
## Calculate phylogenetic tree ------------------------------------
install.packages("GUniFrac")
install.packages("phangorn")
install.packages("magrittr")
install.packages("ade4")


library(GUniFrac)
library(phangorn)
library(magrittr)
library(ade4)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("DECIPHER")
BiocManager::install("rPlant")

# Write all sequences to a fasta file with the ID's corresponding to the rownumber
#uniquesToFasta(Tab_UniFrac_test2,"seqtab_nochim_reunion2019.fasta",ids = (1:ncol(Tab_UniFrac_test2)))
uniquesToFasta(seqtab.nochim2,"Microbiome_phyl_tree_input.fasta",ids = (1:ncol(seqtab.nochim2)))

#Use this tree with CIPRES for RAxML analysis, Geneious for consensus.

# Read the fasta file and put it in a variable
Zcuc <-Biostrings::readAAStringSet("outputs/Microbiome_phyl_tree_input.fasta")

#An alignment is performed and processors are set to five to spead up the proces but not overload the server.
# If you have enough processors available you can elevate the number
Aligned_seq_rarified<-DECIPHER::AlignSeqs(Zcuc, processors = 5)

# The aligned sequences are written to a file to be used by FASTREE
Biostrings::writeXStringSet(Aligned_seq_rarified,"outputs/Aligned_microbiome.fasta")


--------------------------------------------------------------------------------
##                                                                            ##
##                             Statistical pipeline                           ##
##                                                                            ##
--------------------------------------------------------------------------------

# Start statistical analaysis----------------------
  
## Phyloseq and visualisation -----
# 1) OTU file and put in correct format
input_otu <- read.csv("outputs/Frequency_table.csv", row.names=1)
input_otu <- t(input_otu)
View(input_otu)

# 2) metadata file/sample file
otu_metadata <- read.csv("outputs/Input_metadata_unifrac.csv", row.names=1) 
View(otu_metadata)

# 3) taxonomy file
tax_tab <- read.csv("Outputs/Taxonomy_table.csv", row.names=1)
View(tax_tab)
tax_mat <- as.matrix(tax_tab)

#create phyloseq
OTUs = otu_table(input_otu, taxa_are_rows = TRUE)
TAX = tax_table(tax_mat)
samples = sample_data(otu_metadata)


phyl_tree <- ape::read.nexus("Outputs/Tree/Microbiomics_64_consensus_tree.nex")

# check if tree is rooted
is.rooted(phyl_tree)

phyl_tree_rooted <- root(phyl_tree, outgroup = "OTU733", resolve.root = TRUE)
is.rooted(phyl_tree_rooted)
write.tree(phyl_tree_rooted,"outputs/Microbiome_tree_rooted.tre")

phytree = phy_tree(phyl_tree_rooted)

micro_gut <- phyloseq(OTU, samples, TAX, phytree)
micro_gut



## Alpha diversity -------------------

# Calculate Shannon and Inverse simpson
data_shannon <- diversity(input_otu, index = "shannon")
head(data_shannon)
data_invsimpson <- diversity(input_otu, index = "invsimpson")
head(data_invsimpson)

write.table(data_shannon, "outputs/diversity_shannon.txt", sep="\t")
write.table(data_invsimpson, "outputs/diversity_invsimpson.txt", sep="\t")

# Calculate ACE 
calAce <- read.csv("outputs/Complete_filtered_data_64.csv", row.names=1)
AceIndex <- estimateR(t(calAce))
write.table(AceIndex,file="outputs/ACE_index.txt",sep="\t",col.names=TRUE,row.names=TRUE)

# calculate Faith's phylogenetic diversity (PD)
Faith <- pd(input_otu, phyl_tree_rooted, include.root=TRUE)
write.table(Faith, "outputs/diversity_Faith.txt", sep="\t")

# Combine diversity indices in "outputs/Diversity_indices.csv"

### GAD analysis --------------------------------------
#divdat <- read.csv("outputs/Diversity_indices.csv")
divdat <- read.csv("outputs/Diversity_indices.csv")


CON <- as.fixed(divdat$Conditions)
ALT <- as.fixed(divdat$Altitude)
TRT <- as.fixed(divdat$Treatment)
FRT <- as.fixed(divdat$Fruit)
LOC <- as.random(divdat$LOCtest1)


gad(my_model)

#### Shannon -------------------------- 

H_model <-
  lm(
    Shannon ~ TRT + ALT + FRT + TRT * ALT + ALT * FRT + TRT * ALT + TRT * FRT *
      ALT + FRT * LOC %in%
      (TRT * ALT) + LOC %in%
      (TRT * ALT),
    data = divdat
  )


# Cochran's test -> check homogeneity of variances, if significant transform data
C.test(H_model)
# Shannon not significant
# Inverse_Simpson significant
# ACE significant
# PD

# ANOVA Shannon
gad(H_model)

# snk.test(AltitudeShannon_lm, term='Alt', among = NULL, within = NULL)
snk.test(H_model, term='FRT', among = NULL, within = NULL)
# doesn't work for interactions?
snk.test(H_model, term="TRT:ALT", among = "TRT", within = "ALT")

summary( lsmeans(H_model, pairwise ~ FRT*ALT), infer=TRUE)
# only significant differences between fruits in low altitude

# no significant interactions between Treatment and Fruit for Shannon, but fruit is significant

#### Inverse simpson -------------------------- 
#lm = lm(Inverse_Simpson ~ TRT + ALT + FRT + TRT*ALT + ALT*FRT + TRT*FRT + TRT*FRT*ALT, data=divdat) #LOC%in%ALT*TRT
my_model <-
  lm(
    Inverse_Simpson ~ TRT + ALT + FRT + TRT * ALT + ALT * FRT + TRT * ALT + TRT * FRT *
      ALT + FRT * LOC %in%
      (TRT * ALT) + LOC %in%
      (TRT * ALT),
    data = divdat
  )

C.test(my_model) # significant, data needs to be transformed

#### ACE -------------------------- 
#lm = lm(S.ACE ~ TRT + ALT + FRT + TRT*ALT + ALT*FRT + TRT*FRT + TRT*FRT*ALT, data=divdat) #LOC%in%ALT*TRT
my_model <- lm(
  S.ACE ~ TRT + ALT + FRT + TRT * ALT + ALT * FRT + TRT * ALT + TRT * FRT *
    ALT + FRT * LOC %in%
    (TRT * ALT) + LOC %in%
    (TRT * ALT),
  data = divdat
)
C.test(my_model) # significant, data needs to be transformed

#### Faith PD -------------------------- 
#lm = lm(Faith_Phylogenetic_Diversity ~ TRT + ALT + FRT + TRT*ALT + ALT*FRT + TRT*FRT + TRT*FRT*ALT, data=divdat) #LOC%in%ALT*TRT
PD_model <- lm(
  Faith_Phylogenetic_Diversity ~ TRT + ALT + FRT + TRT * ALT + ALT * FRT + TRT * ALT + TRT * FRT *
    ALT + FRT * LOC %in%
    (TRT * ALT) + LOC %in%
    (TRT * ALT),
  data = divdat
)
C.test(PD_model) # Not significant

gad(PD_model)

# snk.test(AltitudeShannon_lm, term='Alt', among = NULL, within = NULL)
snk.test(PD_model, term='TRT', among = NULL, within = NULL)
# doesn't work for interactions?
snk.test(PD_model, term="TRT:ALT", among = "TRT", within = "ALT")

summary(lsmeans(PD_model, pairwise ~ TRT*ALT), infer=TRUE)
summary(lsmeans(PD_model, pairwise ~ FRT * LOC %in% (TRT * ALT)), infer=TRUE)

#treatment and fruit is significant, no significant interactions

### Transform data ----------------------

#Inverse SimpsOn
divdat <- cbind(divdat, sqrt(divdat$Inverse_Simpson))
colnames(divdat)[16] <- "SQRTtransf_InvS"


divdat <- cbind(divdat, (divdat$Inverse_Simpson)^0.25)
colnames(divdat)[17] <- "FourRTtransf_InvS"

divdat <- cbind(divdat, log(divdat$Inverse_Simpson))
colnames(divdat)[18] <- "LNtransf_InvS"

IS_model = lm(SQRTtransf_InvS ~ TRT + ALT + FRT + TRT * ALT + ALT * FRT + TRT * ALT + TRT * FRT *
                ALT + FRT * LOC %in%
                (TRT * ALT) + LOC %in%
                (TRT * ALT), data=divdat)
C.test(lm) # not significant

IS_model = lm(FourRTtransf_InvS ~ TRT + ALT + FRT + TRT * ALT + ALT * FRT + TRT * ALT + TRT * FRT *
                ALT + FRT * LOC %in%
                (TRT * ALT) + LOC %in%
                (TRT * ALT), data=divdat)
C.test(lm) # not significant

IS_model = lm(LNtransf_InvS ~ TRT + ALT + FRT + TRT * ALT + ALT * FRT + TRT * ALT + TRT * FRT *
                ALT + FRT * LOC %in%
                (TRT * ALT) + LOC %in%
                (TRT * ALT), data=divdat)

C.test(IS_model) # significant

gad(IS_model)
snk.test(lm, term='FRT', among = NULL, within = NULL)
snk.test(IS_model, term="TRT:ALT", among = "TRT", within = "ALT")
summary( lsmeans( lm, pairwise ~ FRT*ALT), infer=TRUE)
# only significant differences between cucumber in low altitude
summary( lsmeans( lm, pairwise ~  TRT*ALT), infer=TRUE)
# no significant values

## ACE

divdat <- cbind(divdat, sqrt(divdat$S.ACE))
colnames(divdat)[19] <- "SQRTtransf_ACE"

divdat <- cbind(divdat, (divdat$S.ACE)^0.25)
colnames(divdat)[20] <- "FourRTtransf_ACE"

divdat <- cbind(divdat, log(divdat$S.ACE))
colnames(divdat)[21] <- "LNtransf_ACE"


ACE <- divdat$S.ACE
b <- boxcox(lm(ACE ~ 1))
lambda <- b$x[which.max(b$y)]
lambda
# -0.6666667
new_x_exact <- (ACE ^ lambda - 1) / lambda

divdat <- cbind(divdat, new_x_exact)
colnames(divdat)[22] <- "box_cox_ACE"

ace_model = lm(new_x_exact ~ TRT + ALT + FRT + TRT * ALT + ALT * FRT + TRT * ALT + TRT * FRT *
                 ALT + FRT * LOC %in%
                 (TRT * ALT) + LOC %in%
                 (TRT * ALT), data=divdat)
C.test(ace_model) # not significant


my_model = lm(SQRTtransf_ACE ~ TRT + ALT + FRT + TRT * ALT + ALT * FRT + TRT * ALT + TRT * FRT *
                ALT + FRT * LOC %in%
                (TRT * ALT) + LOC %in%
                (TRT * ALT), data=divdat)
C.test(my_model) # significant

my_model = lm(FourRTtransf_ACE ~ TRT + ALT + FRT + TRT * ALT + ALT * FRT + TRT * ALT + TRT * FRT *
                ALT + FRT * LOC %in%
                (TRT * ALT) + LOC %in%
                (TRT * ALT), data=divdat)
C.test(my_model) # significant

ace_model = lm(LNtransf_ACE ~ TRT + ALT + FRT + TRT * ALT + ALT * FRT + TRT * ALT + TRT * FRT *
                 ALT + FRT * LOC %in%
                 (TRT * ALT) + LOC %in%
                 (TRT * ALT), data=divdat)
C.test(ace_model) # significant

gad(ace_model)
snk.test(ace_model, term='TRT', among = NULL, within = NULL)
# significant differences treatment
snk.test(ace_model, term="TRT:ALT", among = "TRT", within = "ALT")


## Beta diversity --------------------
## use PERMANOVA+ function in PRIMER-e 

# use frequency data, CLR data or genera data: 

## convert to correct dataset: assign NA genera to their respective families ------

GENUS_abundances<-read.csv("outputs/Complete_filtered_data_NA_converted.csv")

dim(GENUS_abundances)
GENUS_abundances$Genus<-as.factor(GENUS_abundances$Genus)
length(levels(GENUS_abundances$Genus))
GENUS_abundances$Family<-as.factor(GENUS_abundances$Family)
length(levels(GENUS_abundances$Family))

ncol(GENUS_abundances)
colnames(GENUS_abundances)
temp_GENUS_abundances<-rowsum(GENUS_abundances[,9:72], GENUS_abundances$Genus, na.rm=TRUE)
dim(temp_GENUS_abundances)
write.csv(temp_GENUS_abundances, file="outputs/GENERA_aggregated.csv")



## Aldex2 --------------------------
### combine on genus --------------------------
Red_gen <- read.csv("outputs/Aldex/Genera_to_combine_high_altitude.csv", row.names=1) 

Red_gen <- Red_gen %>% 
  group_by(Genus) %>% 
  summarise(across(everything(), ~ sum(.x)))

write.xlsx(Red_gen, file = "outputs/Aldex/Input_Aldex_genera_High_altitude.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

filt_data <- read.csv("outputs/Aldex/Input_Aldex_genera_High_altitude.csv", row.names = 1)
aldexmeta <- read.csv("outputs/Aldex/Aldex_metadata_High_altitude.csv", row.names = 1)

Aldex_analysis_f <- aldex(filt_data, conditions = aldexmeta$Fruit, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)
Aldex_analysis_t <- aldex(filt_data, conditions = aldexmeta$Treatment, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)
Aldex_analysis_a <- aldex(filt_data, conditions = aldexmeta$Altitude, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)

write.table(Aldex_analysis_f,file="outputs/Aldex_analysis_fruit.txt",sep="\t",col.names=TRUE,row.names=TRUE)
write.table(Aldex_analysis_t,file="outputs/Aldex_analysis_treatment.txt",sep="\t",col.names=TRUE,row.names=TRUE)
write.table(Aldex_analysis_a,file="outputs/Aldex_analysis_altitude.txt",sep="\t",col.names=TRUE,row.names=TRUE)


### combine on phylum --------------------------

# Red_phyla <- read.csv("outputs/Aldex/ASVs_to_reduce_phylum.csv", row.names=1) 
Red_phyla <- read.csv("outputs/Aldex/Phyla_to_combine_high_altitude.csv", row.names=1)

Red_phyla <- Red_phyla %>% 
  group_by(Phylum) %>% 
  summarise(across(everything(), ~ sum(.x)))
# other classifications (kingdom, family,... must be removed because the sum
# option can only deal with numeric values outside of their group_by selection)

write.xlsx(Red_phyla, file = "outputs/Aldex/Input_aldex_phyla.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

aldexmeta <- read.csv("outputs/Aldex/Metadata_aldex_f_p.csv", row.names = 1)
phyla_dat <- read.csv("outputs/Aldex/Input_aldex_phyla.csv", row.names = 1)


Aldex_analysis_t <- aldex(phyla_dat, conditions = aldexmeta$Management, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)

write.xlsx(Aldex_analysis_t, file = "outputs/Aldex/Output_aldex_phyla_HIGH.xlsx" , sheetName = "Management_HIGH", 
           col.names = TRUE, row.names = TRUE, append = FALSE)


### combine on family --------------------------


Red_fam <- read.csv("outputs/Aldex/Families_to_combine_high_altitude.csv", row.names=1)

Red_fam <- Red_fam %>% 
  group_by(Family) %>% 
  summarise(across(everything(), ~ sum(.x)))
# other classifications (kingdom, family,... must be removed because the sum
# option can only deal with numeric values outside of their group_by selection)

write.xlsx(Red_fam, file = "outputs/Aldex/Input_aldex_families.xlsx" , sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

aldexmeta <- read.csv("outputs/Aldex/Metadata_aldex_f_p.csv", row.names = 1)
phyla_dat <- read.csv("outputs/Aldex/Input_aldex_families.csv", row.names = 1)


Aldex_analysis_t <- aldex(phyla_dat, conditions = aldexmeta$Management, mc.samples = 1000, test="t", effect=TRUE, denom="all", verbose=TRUE)

write.xlsx(Aldex_analysis_t, file = "outputs/Aldex/Output_aldex_family_HIGH.xlsx" , sheetName = "Management_HIGH", 
           col.names = TRUE, row.names = TRUE, append = FALSE)


