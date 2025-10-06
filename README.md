# newt_microbiome
Analysis of the microbiome of Calotriton arnoldi

**Preparing data for Dada**

The first step taken was to install all necessary programs with the following commands: 

**Package installing**

Installing cutadapt through the terminal with:
```
conda install -c bioconda cutadapt
```
Also created a short script to run through the terminal with Rscript:
```
#Package_install.R
if (!requireNamespace("BiocManager", quietly = TRUE))
 	install.packages("BiocManager", repos = "http://cran.us.r-project.org") #repository was necessary in this case
 	BiocManager::install()

if (!requireNamespace("BiocManager", quietly = TRUE))
 	install.packages("BiocManager", repos = "http://cran.us.r-project.org")
	BiocManager::install (c("dada2", "DECIPHER", "phangorn", "phyloseq", "dplyr"))

install.packages("ggplot2", repos = "http://cran.us.r-project.org")
```

The package "dada2" was giving problems to install, so the following alternative way was utilised

```
#Package_install_2.R
install.packages("devtools", repos = "http://cran.us.r-project.org")
library("devtools")
devtools::install_github("benjjneb/dada2", ref="v1.16")
```

Once this was done, the next step was eliminating adaptors and primers.
Before trimming, a _multiqc_ was done for all the files to check that everything was in order before starting with the command(the program _fastqc_ was already installed before doing this): 

`multiqc *` 

It turnes out that the sequencing company had already eliminated adaptors, so only primers were left to trim.

This was done with _cutadapt_ through the terminal with the following script:

**cutadapt.sh**

```
#Create a list with names
ls *"_L001_R1_001.fastq" | cut -f 1-2 -d "_" > samples

#Eliminating primers
for sample in $(cat samples)
do

    echo "On sample: $sample"

    cutadapt -g ^GTGYCAGCMGCCGCGGTAAC \ #Forward primer. "-g" was used instead of "-a" (and same with "-G" instead of "-A") because primers weren't linked
    -G ^GGACTACNVGGGTWTCTAAT \ #Reverse primer
    -m 200 -M 300 --discard-untrimmed \
    -o ${sample}_L001_R1_001_trimmed.fastq.gz -p ${sample}_L001_R2_001_trimmed.fastq.gz \
    ${sample}_L001_R1_001.fastq.gz ${sample}_L001_R2_001.fastq.gz \
    >> cutadapt_primer_trimming_stats.txt 2>&1

done
```

The generated file ([cutadapt_primer_trimming_stats.txt](https://github.com/SergiTullochJimenez/Newt_microbiome/blob/main/cutadapt_primer_trimming_stats.txt)) is available also, but a summary can be obtained with the following command:


```
paste samples <(grep "passing" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") <(grep "filtered" cutadapt_primer_trimming_stats.txt | cut -f3 -d "(" | tr -d ")") `
 
```

The result for the trimmed sequences per sample was the following:

![image](https://github.com/SergiTullochJimenez/Newt_microbiome/assets/129175590/0a6cbe78-0e41-4f76-aaa7-90089022bf01)



**Sequence processing**

_Dada2.R_

Create working directory and set it up
```
setwd("~/Desktop/Sergi/microbiome")

library(dada2); packageVersion("dada2")
library(ggplot2)

path <- "~/Desktop/Sergi/microbiome/data" # CHANGE ME to the directory containing the fastq files after unzipping.
list.files(path)
```

Read files into R
Check for sample patterns. We need to have available only the names that we are looking for (eliminate .html and .zip)

```
fnFs <- sort(list.files(path, pattern="_L001_R1_001_trimmed.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_L001_R2_001_trimmed.fastq", full.names = TRUE))
```
Extract sample names in R
Carefull, look at how the file comes out
```
sample.names <- sapply(strsplit(basename(fnFs), "_L001"), `[`, 1)
sample.names
```
Plot quality profiles. I chose 3 random samples
```
[plotQualityProfile(fnFs[51:53])](https://github.com/SergiTullochJimenez/Newt_microbiome/blob/main/Quality_plot_F.pdf)
[plotQualityProfile(fnRs[51:53])](https://github.com/SergiTullochJimenez/Newt_microbiome/blob/main/Quality_plot_R.pdf)
```
Filtering and Trimming: According to the plots we decided to go for 210 length for forward and reverse

```
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, trimLeft=5, truncLen=c(210,210),
                     maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=2, matchIDs=TRUE)
head(out)
```

Errors
```
set.seed(100)
errF <- learnErrors(filtFs, nbases=1e8, multithread=2, randomize=TRUE)
errR <- learnErrors(filtRs, nbases=1e8, multithread=2, randomize=TRUE)
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

```
[Forward plot](https://github.com/SergiTullochJimenez/Newt_microbiome/blob/main/Errors_F.pdf)
[Reverse plot](https://github.com/SergiTullochJimenez/Newt_microbiome/blob/main/Errors_R.pdf)


Infer sequences
```
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
sam.names <- sapply(strsplit(basename(filtFs), "_"), `[`, 1)
names(derepFs) <- sam.names
names(derepRs) <- sam.names
ddFs <- dada(derepFs, err=NULL, selfConsist=TRUE)
ddRs <- dada(derepRs, err=NULL, selfConsist=TRUE)
plotErrors(ddFs)
plotErrors(ddRs)
dadaFs <- dada(derepFs, err=ddFs[[1]]$err_out, pool=TRUE, multithread=2)
dadaRs <- dada(derepRs, err=ddRs[[1]]$err_out, pool=TRUE, multithread=2)
dadaFs[[1]]
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
```
[Plot for derepF](https://github.com/SergiTullochJimenez/Newt_microbiome/blob/main/Infer_F.pdf)
[Plot for derepR](https://github.com/SergiTullochJimenez/Newt_microbiome/blob/main/Infer_R.pdf)

Inspect to see how it went and export the results ([mergers.xlsx](https://github.com/SergiTullochJimenez/Newt_microbiome/blob/main/mergers.xlsx))

```
head(mergers[[1]])

library ("xlsx")

i<-1
while (i <= length(sample.names)){
  write.xlsx(mergers[[i]], file = "mergers.xlsx",
             sheetName = sample.names[i], append = TRUE)
  i<-i+1
}
```
Construct table, remove chimeras and save the created R objects
```
seqtab.all <- makeSequenceTable(mergers)
seqtab <- removeBimeraDenovo(seqtab.all)
dim(seqtab)
table(nchar(getSequences(seqtab))) # Inspect distribution of sequence lengths

asv_seqs <- colnames(seqtab)
asv_headers <- vector(dim(seqtab)[2], mode="character")
for (i in 1:dim(seqtab)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "~/ASVs.fa")

asv_tab <- t(seqtab)
write.table(asv_tab, "~/Desktop/Sergi/microbiome/ASVs_counts-seqs.tsv", sep="\t", quote=F, col.names=NA)
```
**Taxonomy**

It is important to download the latest SILVA database and unzip the file
```
setwd("~/Desktop/Sergi/microbiome") #REMEMBER TO LOOK AT IT!

library (dplyr)
library (DECIPHER)
library(dada2); packageVersion("dada2")

#Consulting SILVA database
ref_fasta <- "/Users/imac2016/Desktop/Sergi/TFM/Microbiome_analysis/BonacoltaEMPV4/filtered/silva_nr99_v138.1_train_set.fa.gz" #Check what the path is
taxa <- assignTaxonomy(seqtab, refFasta=ref_fasta, multithread=2)
colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")

#Species assignment

taxa <- addSpecies(taxa, "/Users/imac2016/Desktop/Sergi/TFM/Microbiome_analysis/BonacoltaEMPV4/filtered/silva_species_assignment_v138.1.fa.gz")
taxa.print <- taxa

write.xlsx(taxa.print, file = "taxed.xlsx", 
           sheetName= "taxed", append=TRUE)

asv_tax <- taxa
row.names(asv_tax) <- sub(">", "", asv_headers)
write.table(asv_tax, "~/Desktop/Sergi/microbiome/ASVs_taxonomy.tsv", sep="\t", quote=F, col.names=NA)

asv_tabtax <- cbind(asv_tab,asv_tax)
row.names(asv_tabtax) <- sub(">", "", asv_headers)
write.table(asv_tabtax, "~/Desktop/Sergi/microbiome/ASVs_counts_taxonomy.tsv", sep="\t", quote=F, col.names=NA)

```
At this point a tree would also be done, but it still hasn't been done due to the computational intensity and unavailablity (up to this day) of powerful enough computers.

**Data filtering**

```
setwd("~/Desktop/Sergi/microbiome")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq", version = "3.9")

library("phyloseq"); packageVersion("phyloseq")
library("ggplot2"); packageVersion("ggplot2")
library("viridis"); packageVersion("viridis")
library("tidyr"); packageVersion("tidyr")
library("plyr"); packageVersion("plyr")
library("dplyr"); packageVersion("dplyr")
library("tidyr"); packageVersion("tidyr")
```
Combinig data into an object
```
SV <- as.matrix(read.delim("~/Desktop/Sergi/microbiome/analysis/ASVs_counts.tsv", row.names = 1, check.names= TRUE, sep = "\t", header = TRUE))
tax <-as.matrix(read.delim("~/Desktop/Sergi/microbiome/analysis/ASVs_taxonomy.tsv", row.names = 1, header = TRUE, sep = "\t"))
map <- read.delim("~/Desktop/Sergi/microbiome/analysis/DATABASE_SEQUENCED.csv", sep =",", header = TRUE)

#tree <- ape::read.tree(file = "~/Desktop/Sergi/microbiome/tree.nw")
```
Row names weren't right in "map", changed them in this way
```
row.names(map)
samples<-map$EXTRACTION.CODE
rownames(map)<-samples

rownames(map) = make.names(samples, unique=TRUE)

ps = phyloseq(otu_table(SV, taxa_are_rows=TRUE), 
              sample_data(map), 
              tax_table(tax) 
              #phy_tree(tree)
              )
```



Now filtering by eliminating unwanted taxa (for now). 
As said, it might be interesting to explore the chloroplast sequences or the mitochondrial ones in the future
```
ps_filtered <- subset_taxa(ps, (Order!="Chloroplast"| is.na(Order)))
ps_filtered <- subset_taxa(ps_filtered, (Family!="Mitochondria"| is.na(Family)))
ps_filtered <- subset_taxa(ps_filtered, (Kingdom!="Eukaryota"| is.na(Kingdom)))
ps_filtered <- subset_taxa(ps_filtered, (Class!="Embryophyceae"| is.na(Class)))
```

Now to filter for minimum 250 reads and for the ASV that are present in the 99% of abundance in at least one sample:
```
ps_filtered_2 <- prune_samples(sample_sums(ps_filtered)>=250, ps_filtered)
f1<- filterfun_sample(topf(0.99))
wh1 <- genefilter_sample(ps_filtered, f1, A=2)
ps_filtered_2 <- prune_taxa(wh1, ps_filtered)
```

This object was saved in an RDS file 
```
saveRDS(ps_filtered_2, 'ps_filtered_2.rds')
```

**Decontamination**
```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("decontam")

library("decontam")
```
Following the guides in the Del Campo git-hub, the decontamination was done with the following code in two steps:
The first one included extraction controls, applicable to all samples
```
decontam_result_0.099 <- isContaminant(
  ps_filtered_2,
  neg = sample_data(ps_filtered_2)$Negative,
  method = "prevalence",
  detailed = TRUE
)
write.csv(decontam_result_0.099, "~/Desktop/Sergi/microbiome/analysis/decontam_results_0.099.csv")
```
The [decontamination]() file contained probable contamination ASVs, which were copied and pasted into a vector in R with 
```
badTaxa<-read.delim(pipe("pbpaste"))
```
A new phyloseq object was created without possible contaminants and their samples
```
badTaxa <- bt
allTaxa <- taxa_names(ps_filtered_2)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]

ps_filtered_decontaminated <- prune_taxa(allTaxa, ps_filtered_2)
ps_filtered_decontaminated = subset_samples(ps_filtered_decontaminated, EXTRACTION.CODE != "CTRL1" & EXTRACTION.CODE != "CTRL3") #Remove the two blanks
saveRDS(ps_filtered_decontaminated, 'ps_filtered_decontaminated.rds')
```
Which looked like this:

phyloseq-class experiment-level object

otu_table()   OTU Table:         [ 8209 taxa and 190 samples ]

sample_data() Sample Data:       [ 190 samples by 24 sample variables ]

tax_table()   Taxonomy Table:    [ 8209 taxa by 8 taxonomic ranks ]

The second decontamination was done in batches, as each batch had it's own water, which was sequenced and used for this decontamination
```
ps_filtered_decontaminated_2 = subset_samples_no_zero(ps_filtered_decontaminated, Species == "Calotriton arnoldi") #Remove unwanted sequences

otu_table()   OTU Table:         [ 7481 taxa and 138 samples ]
sample_data() Sample Data:       [ 138 samples by 24 sample variables ]
tax_table()   Taxonomy Table:    [ 7481 taxa by 8 taxonomic ranks ]

ps_merged <- merge_samples(ps_filtered_decontaminated_2,
                           group = "Specimen.Code",
                           fun = "sum")
#Get the original metadata
meta <- as.data.frame(sample_data(ps_filtered_decontaminated_2))

#Make sure it's a plain data.frame, not a tibble
meta <- as.data.frame(meta)

#Merge metadata by Specimen.Code (pick first non-NA entry per group)
meta_merged <- meta %>%
  group_by(Specimen.Code) %>%
  summarise(across(everything(), ~ first(na.omit(.))), .groups = "drop") %>%
  as.data.frame()

#Set rownames to Specimen.Code (must match merged phyloseq sample names)
rownames(meta_merged) <- meta_merged$Specimen.Code

#Check that names match
identical(sort(rownames(meta_merged)), sort(sample_names(ps_merged)))
# This must return TRUE

#Attach metadata back
sample_data(ps_merged) <- sample_data(meta_merged)

```
phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 7481 taxa and 82 samples ]
sample_data() Sample Data:       [ 82 samples by 24 sample variables ]
tax_table()   Taxonomy Table:    [ 7481 taxa by 8 taxonomic ranks ]

wh2 <- genefilter_sample(
  ps_merged,
  filterfun_sample(function(x) x >= 2),
  A = 2
)

ps_filtered_2 <- prune_taxa(wh2, ps_merged)

The final phyloseq object will look like this:

phyloseq-class experiment-level object
otu_table()   OTU Table:         [ 5996 taxa and 82 samples ]
sample_data() Sample Data:       [ 82 samples by 24 sample variables ]
tax_table()   Taxonomy Table:    [ 5996 taxa by 8 taxonomic ranks ]

