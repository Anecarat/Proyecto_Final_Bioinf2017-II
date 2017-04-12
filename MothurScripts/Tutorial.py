#Mothur tutorial
#Considering we have a directory called data with the directories called mothur, process, raw and references.  

#Get into process directory
cd data/process

#Start mothur
mothur

#Set input and output directories
set.dir(input=../raw) #raw directory that contains all our raw files
set.dir(output=.) #actual directory

#Make the file with the information of our reads and rename it.  
make.file(inputdir=../raw, type=gz)
system(mv fileList.paired.file stability.files)

#The file format is a tab delimited file with the sample name in the first column, the forward read in the second column and the paired end read in the third column.
