## Author: Víctor Taracena
## email: vtaracena@gmail.com

#### This script is based on Maggie Wagner's script 'main_code.R' from her Dryad repository. 

#### Dryad Repository: 
###### Wagner MR, Lundberg DS, del Rio TG, Tringe SG, Dangl JL, Mitchell-Olds T (2016) 
###### Data from: Host genotype and age shape the leaf and root microbiomes of a wild perennial plant. 
###### Dryad Digital Repository. http://dx.doi.org/10.5061/dryad.g60r3

#### Article: 
###### Wagner MR, Lundberg DS, del Rio TG, Tringe SG, Dangl JL, Mitchell-Olds T (2016) 
###### Host genotype and age shape the leaf and root microbiomes of a wild perennial plant. 
###### Nature Communications 7: 12151. http://dx.doi.org/10.1038/ncomms12151

## Preparing the environment
#### If you have all the files in one folder you need to create some folders and move the files
#### into them. 
#### Note: This commands works in Linux / MacOS, if you are using Windows you'll have to change
#### to the commands that Windows use or do it manually. 

system('mkdir Raw_Data')
system('mv contaminants.fasta raw_data')
system('mv Eco_Field_copynumest_forR.txt raw_data')
system('mv Ecotypes_field_glucosinolates.txt raw_data')
system('mv OTUrepSeqs97.fa raw_data')
system('mv otuTable97.txt.bz2 raw_data')
system('mv otuTable99.txt.bz2 raw_data')
system('mv phylogeny.tre raw_data')
system('mv plant_key.txt raw_data')
system('mv site_coords_Ecotypes.txt raw_data')
system('mv SMD.txt raw_data')
system('mv soildata.txt raw_data')
system('mv taxAssignments97.txt raw_data')
system('mv taxAssignments99.txt.bz2 raw_data')

#### You need to unzip the .bz2 files
############################## Aquí me quedé!!! ###############
system('bzip2 -dk Raw_Data/otuTable97.txt.bz2')


####### Separate code sections into subdirectories #######
system('mkdir foldchange')
system('mkdir GLS')
system('mkdir heritability')
system('mkdir higher_tax_levels')
system('mkdir LMMs')
system('mkdir otu99pct')
system('mkdir rel_abund')
system('mkdir succession')
system('mkdir uUF')
# Move scripts into respective folders:
system('mv foldchange.R foldchange')
system('mv GLS.R GLS')
system('mv heritability.R heritability')
system('mv countmodels_vstLMM.R heritability')
system('mv populate_vstLMM.R heritability')
system('mv higher_tax_levels.R higher_tax_levels')
system('mv LMMs.R LMMs')
system('mv otu99pct.R otu99pct')
system('mv rel_abund.R rel_abund')
system('mv succession.R succession')
system('mv uUF.R uUF')
####### Make subdirectories for plots, tables, and intermediate data files #######
system('mkdir plots')
system('mkdir tables')
system('mkdir intermediate_data')