## Author: Víctor Taracena
## email: vtaracena@gmail.com

## Based on the work of Maggie Wagner and collaborators and the bioconductor
## workflow by Callaghan et al. (more info in 0.MasterScript.R)

                      ## 2.AssignTaxonomy ##
## In this script you will find commands for assign taxonomy information to each OTU,
## how to build a phylogenetic tree of the samples and create a phyloseq object to relate the taxonomic
## data with the tree, taxonomic and sample meta data. 
## Also, you will find how to eliminate contaminant sequences (if you have a reference for them), 
## and other sequences that belong to animals or plants (like chloroplasts and mitochondrias) or
## unclassified ones. 

#### README IMPORTANT ####
## Para fines del proyecto de Bioinformática, no se planea correr
## una parte de este script, dado el peso máximo de los archivos y el tiempo. 
## Sólo es para mostrar parte del workflow que seguí. Se puede correr a partir
## del README IMPORTANT 2. 

#### NOTE #### 
## If you closed the R session and didnt save your workspace you will need to
## load the OTUtable again to continue with this script. 
OTUtab<-read.delim("RawData/OTUtable.txt", sep="\t", head=T)

#### Assign taxonomy ####
taxa<-assignTaxonomy(OTUtab, "path/to/RefSequenceDataBase.fa", multithread=TRUE)
colnames(taxa)<-c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
write.table(taxa, "RawData/TaxAssignment.txt", sep="\t")
#### Create tree ####
seqs <- getSequences(OTUtab)
names(seqs) <- seqs
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)

fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
write.tree(fitGTR, file="RawData/phylogeny.tre")

#### README IMPORTANT 2 ####
## A partir de aquí se puede correr el script para el proyecto de Bioinfo. 
## Tienen que cargar los siguientes archivos:
OTUtab<-read.delim("RawData/OTUtable.txt", sep="\t", head=T)
tax<-read.delim("RawData/TaxAssignment.txt", sep="\t", head=T)
tree<-read.tree("RawData/phylogeny.tre")
SMD<-read.delim("RawData/SampleMD.txt", sep="\t", head=T)

#### Creating phyloseq object ####
OTUtab<-data.matrix(OTUtab)
tax<-as.matrix(tax)
bacteria<-phyloseq(otu_table(OTUtab,taxa_are_rows=TRUE),sample_data(SMD),tax_table(tax),tree)

#### Removing contaminant OTU sequences ####
contams<-readDNAStringSet("RawData/contaminants.fasta",format="fasta")
## Extract the OTU_IDs of the contaminants:
## Extract sequence names
contams<-names(contams) %>% strsplit(.," ") %>% unlist() 
## Only keep the OTU_IDs
contams<-contams[grep("OTU_",contams)] 
## Remove OTUs that not belong to bacteria
bacteria<-subset_taxa(bacteria, !is.na(Phylum) & 
                             !Phylum %in% c("", "uncharacterized") &
                             Phylum!=c('Unassigned') & 
                             Class!='Chloroplast' & 
                             Family !='mitochondria')
## Remove contaminant OTUs
bacteria<-prune_taxa(setdiff(taxa_names(bacteria),contams),
                          bacteria) 
## Remove taxa with 0 observations
bacteria<-prune_taxa(taxa_sums(bacteria)>0,bacteria) 
## Adding number of usable reads for each samples to sample_metadata
sample_data(bacteria)$UsableReads<-sample_sums(bacteria)

#### Useful data about the sequences ####
## How many observations?
mean(sample_data(bacteria)$UsableReads)
sd(sample_data(bacteria)$UsableReads)

## How many OTUs?
ntaxa(bacteria)

            ## Please continue to 3.DiversityAnalyses.R script ##