# README
## Workflow: From 16S rRNA gene sequencing raw data to community analyses. 

The present workflow contains the scripts with the necessary commands to load, trimm, clean, generate OTUs, create phyloseq objects and some exploratory analysis for 16S rRNA gene sequences of bacterial communities using only R. 
This workflow (and the data used for running the examples) is based in Wagner's and collaborators' (2016) research of the effect of host genotype, age and site in the microbial communities of leaves and roots of *Boecheria stricta* ([Dryad repository](http://dx.doi.org/10.1038/ncomms12151)) and the [bioconductor workflow](https://f1000research.com/articles/5-1492/v1) published by Callaghan and his collaborators to process raw sequencing data in R.

## Content

Here you will find the **README** that you just opened and you are reading at the moment and a directory called **OPENME**.
Obviously, you have to open the **OPENME** directory where you will find the next files:
- **0.MasterScript:** The first script you will need to run. This one will organize and create the environment to work the consequent scripts. Also, you will find the session information of the PC where I ran this workflow.  
- **Initial:** This directory contains the 'raw data' and scripts you will use later. 

## What do you need to do first?

You need to open the **OPENME** directory, load the **0.MasterScript.R** into R, set your working directory to source file location and run it. This will create all the directories you need and move the files to where they belong. 

## And then?

Well, after that you now have the environment to run everything. Inside the **bin** directory you will find these four scripts:

- **1.PreProcessing:** Everything you need to do first with your raw sequences; from barcode, primer and quality trimming,to OTU generation. 
- **2.AssignTaxonomy:** Taxa assignment to the OTUs, creation of phyloseq objects and data base cleaning. 
- **3.DiversityAnalyses:** Alpha diversity exploratory analyses and other useful plots. Also, the map of the site. 
- **4.MultivariateAnalyses:** PCoA of the bacteria communities (analysis and plot).

You can read more about what each script does inside each script, at the beginning. 

## README IMPORTANT

Due to the size of the files, I don't provide the raw sequences to run the first script (1.PreProcessing.R) and half of the second (2.AssignTaxonomy.R), so if you want to run the scripts with the data for an example of how this workflow works you need to run the **0.MasterScript.R** first (don't forget it) and then run the **2.AssignTaxonomy.R** script after the **README IMPORTANT 2**. 

The first and second scripts also have this information to know where to start with the data I provide. The data I used to run the exampes are a subset of the samples taken by Wagner and collaborators ([2016](http://dx.doi.org/10.1038/ncomms12151)) and don't belong to me.


