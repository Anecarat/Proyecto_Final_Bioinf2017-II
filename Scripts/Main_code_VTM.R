## Author: Víctor Taracena
## email: vtaracena@gmail.com

## This script is based on Maggie Wagner's script 'main_code.R' from her Dryad repository. 

#### Dryad Repository: 
###### Wagner MR, Lundberg DS, del Rio TG, Tringe SG, Dangl JL, Mitchell-Olds T (2016) 
###### Data from: Host genotype and age shape the leaf and root microbiomes of a wild perennial plant. 
###### Dryad Digital Repository. http://dx.doi.org/10.5061/dryad.g60r3

#### Article: 
###### Wagner MR, Lundberg DS, del Rio TG, Tringe SG, Dangl JL, Mitchell-Olds T (2016) 
###### Host genotype and age shape the leaf and root microbiomes of a wild perennial plant. 
###### Nature Communications 7: 12151. http://dx.doi.org/10.1038/ncomms12151


                                    ###### Main Code ######

## Preparing the environment:
#### If you have all the files in one folder you need to create some folders and move the files
#### into them. 
#### Note: This commands works in Linux / MacOS, if you are using Windows you'll have to change
#### to the commands that Windows use or do it manually. 

system('mkdir Raw_Data')
system('mv contaminants.fasta Raw_Data')
system('mv Eco_Field_copynumest_forR.txt Raw_Data')
system('mv Ecotypes_field_glucosinolates.txt Raw_Data')
system('mv OTUrepSeqs97.fa Raw_Data')
system('mv otuTable97.txt.bz2 Raw_Data')
system('mv otuTable99.txt.bz2 Raw_Data')
system('mv phylogeny.tre Raw_Data')
system('mv plant_key.txt Raw_Data')
system('mv site_coords_Ecotypes.txt Raw_Data')
system('mv SMD.txt Raw_Data')
system('mv soildata.txt Raw_Data')
system('mv taxAssignments97.txt Raw_Data')
system('mv taxAssignments99.txt.bz2 Raw_Data')

#### You need to unzip the .bz2 files

system('bzip2 -dk Raw_Data/otuTable97.txt.bz2')


#### Create the Scripts folder and move the scripts to it. 

system('mkdir Scripts')
system('mv foldchange.R Scripts')
system('mv GLS.R Scripts')
system('mv heritability.R Scripts')
system('mv countmodels_vstLMM.R Scripts')
system('mv populate_vstLMM.R Scripts')
system('mv higher_tax_levels.R Scripts')
system('mv LMMs.R Scripts')
system('mv otu99pct.R Scripts')
system('mv rel_abund.R Scripts')
system('mv succession.R Scripts')
system('mv uUF.R Scripts')

#### Make subdirectories for plots, tables, and intermediate data files

system('mkdir plots')
system('mkdir tables')
system('mkdir intermediate_data')

## Clear workspace

rm(list=ls())

## Session Info

sessionInfo()
"R version 3.3.3 (2017-03-06)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 16.04.2 LTS

locale:
[1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
[3] LC_TIME=es_MX.UTF-8        LC_COLLATE=en_US.UTF-8    
[5] LC_MONETARY=es_MX.UTF-8    LC_MESSAGES=en_US.UTF-8   
[7] LC_PAPER=es_MX.UTF-8       LC_NAME=C                 
[9] LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=es_MX.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods  
[7] base     

loaded via a namespace (and not attached):
[1] tools_3.3.3"

## Load source file
source('ecotypes_source_VTM.R')

                ###### Load individual components of Phyloseq object #######

#### You need to unzip the otuTable97.txt.bz2 files

system('bzip2 -dk ../Raw_Data/otuTable97.txt.bz2')

## Import raw OTU table (not copy number corrected yet)

fullOTU<-as.data.frame(read.table('../Raw_Data/otuTable97.txt',sep='\t',header=TRUE))
#### Store OTU_IDs as row names
rownames(fullOTU)<-paste('OTU_',as.character(fullOTU[,1]),sep='')
#### Remove 'OTU_ID' column
fullOTU<-fullOTU[,-1] 
#### Convert to data matrix format
fullOTU<-data.matrix(fullOTU) 

## Load sample metadata

fullSMD<-as.data.frame(read.table('../Raw_Data/SMD.txt',sep='\t',header=TRUE))
#### Drop columns that won't be used:
fullSMD<-fullSMD[,!names(fullSMD) %in% c("oldPlate")]
#### Clean up data: (change Year, Age to factors and place endogenous plants into Ecotypes 
#### experiment)
endog<-filter(fullSMD,Genotype=='endog')
endog$Age<-as.factor('endog')
endog$Cohort<-as.factor('endog')

fullSMD<-filter(fullSMD,Genotype!='endog')%>%rbind(.,endog)
fullSMD$Age<-factor(fullSMD$Age)
fullSMD$Harvested<-factor(fullSMD$Harvested)
fullSMD$Experiment<-plyr::mapvalues(fullSMD$Experiment,c('endog'),c('ecotypes'))

rownames(fullSMD)<-fullSMD$SampleID

## Load taxonomy assignments

taxfile<-as.matrix(read.table('../Raw_Data/taxAssignments97.txt',sep='\t',header=TRUE))
dimnames(taxfile)[[1]]<-taxfile[,1]  # make OTU_IDs the row names
taxfile<-taxfile[,-c(1,2,8)]  # get rid of taxonomy, OTU.ID, Confidence columns
taxfile<- as.data.frame(taxfile)
# Fix problem noticed during downstream data analysis:
filter(taxfile,Family=='Chromatiaceae') # Family Chromatiaceae classified into 2 different orders: Chromatiales (correct) and Alteromonadales (incorrect)
badChromatiaceae<-subset(taxfile,Family=='Chromatiaceae' & Order=='Alteromonadales')
# Best BLAST matches for 'Chromatiaceae-Alteromonadales':
# OTU_4719 = Rheinheimera sp. = Chromatiaceae, Chromatiales
# OTU_29684 = Alishewanella sp. = Alteromonadaceae, Alteromonadales
# OTU_22285 = Rheinheimera sp. = Chromatiaceae, Chromatiales

# Fix manually: 
badChromatiaceae$Order<-ifelse(rownames(badChromatiaceae)%in%c('OTU_4719','OTU_22285'),'Chromatiales','Alteromonadales')
badChromatiaceae$Family<-ifelse(rownames(badChromatiaceae)=='OTU_29684','Alteromonadaceae','Chromatiaceae')
taxfile<-subset(taxfile,!(Family=='Chromatiaceae' & Order=='Alteromonadales')) %>%
  rbind(.,badChromatiaceae) 
taxfile<-as.matrix(taxfile)
rm(badChromatiaceae)

## load phylogeny
phyfile<-read.tree(file="../Raw_Data/phylogeny.tre")

####### Load reference sequences and contaminants ######
## load reference sequences
refseq<-readDNAStringSet("../Raw_Data//OTUrepSeqs97.fa",format="fasta")
contams<-readDNAStringSet("../Raw_Data//contaminants.fasta",format="fasta")
# Extract the OTU_IDs of the contaminants:
contams<-names(contams) %>% strsplit(.," ") %>% unlist() # extract sequence names
contams<-contams[grep("OTU_",contams)] # only keep the OTU_IDs

