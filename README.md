# Newt_microbiome
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

The generated file cutadapt_primer_trimming_stats.txt is available also, but a summary can be obtained with the following command:


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
plotQualityProfile(fnFs[51:53])
plotQualityProfile(fnRs[51:53])
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


#Infer sequences
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
Inspect to see how it went and export the results

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








