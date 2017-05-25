## Author: Víctor Taracena
## email: vtaracena@gmail.com

## Based on the work of Maggie Wagner and collaborators and the bioconductor
## workflow by Callaghan et al. (more info in 0.MasterScript.R)

                        ## 1.PreProcessing ##
## This script has the necessary instructions to process raw sequences 
## like the ones you obtain from Illumina MiSeq sequencer. 
## Here you can find commands to trim out barcodes and primers, and low quality sequences,
## merge different samples to one large data base, the elimination of chimeras and 
## the generation of operational taxonomic units (OTUs). 

#### README IMPORTANT ####
## Para fines del proyecto de Bioinformática, no se planea correr
## esta parte del script, dado el peso máximo de los archivos. 
## Sólo es para mostrar parte del workflow que seguí. 


#### Relative paths to the samples ####
seq <- sort(list.files("RawData", full.names = TRUE))
seqFs <- seq[grepl("R1", seq)]
seqRs <- seq[grepl("R2", seq)]

#### Primer trimming ####
## Input
bar <- readFastq("path/to/file.fastq")
## Extract bar code
code <- narrow(sread(bar), 1, 8) 
## Subset one bar code
aBar <- bar[code == "Aquí va el barcode de las muestras"] 
## Remove bar code
noBar <- narrow(aBar, 11, width(aBar))
pcrPrimer <- "Aquí va el primer de las muestras"
trimmed <- trimLRPatterns(pcrPrimer, Lfixed=FALSE)
## Output
writeFastq(trimmed, "pathe/to/output/file/trimmed.fastq")

#### Trim and filter ####
## Plot Quality Profile for each .fastq file
ii <- sample(length(seqFs), 2)
for(i in ii) { print(plotQualityProfile(seqFs[i]) + ggtitle("Fwd")) }
for(i in ii) { print(plotQualityProfile(seqRs[i]) + ggtitle("Rev")) }
## Stablish relative paths to trimmed and filtered outputs
system('mkdir RawData/filt')
filtFs <- file.path("RawData/filt", basename(seqFs))
filtRs <- file.path("RawData/filt", basename(seqRs))
## Export trimmered and filter .fastq files
out<-filterAndTrim(seqFs, filtFs, seqRs, filtRs, truncLen=c(240,160),
                   maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                   compress=TRUE, multithread=TRUE)

#### Infer sequence variants ####
## Dereplicate .fastq files
derepFs <- derepFastq(filtFs)
derepRs <- derepFastq(filtRs)
sam.names <- sapply(strsplit(basename(filtRs), "_"), `[`, 1)
names(derepFs) <- sam.names
names(derepRs) <- sam.names

## Note: [] number of samples in the derep*s objects
ddF <- dada(derepFs[1:2], err=NULL, selfConsist=TRUE) 
ddR <- dada(derepRs[1:2], err=NULL, selfConsist=TRUE)

plotErrors(ddF, nominalQ=T)
plotErrors(ddR, nominalQ=T)

dadaFs <- dada(derepFs, err=ddF[[1]]$err_out, pool=TRUE, multithread=T)
dadaRs <- dada(derepRs, err=ddR[[1]]$err_out, pool=TRUE, multithread=T)

#### Create OTU table and remove chimeras ####
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs)
OTUtable.all <- makeSequenceTable(mergers)
OTUtable <- removeBimeraDenovo(OTUtable.all)

#### Checking how many sequences were taken off in each step ####
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab.all), rowSums(seqtab))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sam.names
head(track)

#### Save OTU file ####
write.table(OTUtable, "RawData/OTUtable.txt", sep="\t")

            ## Please continue to script 2.AssignTaxonomy ##
