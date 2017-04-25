#MiSeq 16S rRNA analyzed

#Start mothur
mothur

#make.contigs reads a forward fastq file and a reverse fastq file and outputs new fasta and report files.
#ffastq and rfastq are forward and reverse fatsq archives.
#The findex and rindex parameters are used to provide a forward/reverse index files to process. If you use an index file, you must provide an oligos file.The index file is a fastq file containing barcodes for the reads.
#oligos is the file that has the primer and barcode information.
#pdiffs is maximum number of differences to the primer sequence, default=0. bdiffs is maximum number of differences to the barcode sequence, default=0. tdiffs is maximum total number of differences to the barcode and primer.
#processors tells how many processors you would like to use.
make.contigs(ffastq=Reads_R1.fastq, rfastq=Reads_R2.fastq, findex=Indx_I1.fastq, oligos=oligos_16SMiSeqv3.oligos, pdiffs=2, bdiffs=1,processors=4)

#summary.seqs command will summarize the quality of sequences in an unaligned or aligned fasta-formatted sequence file.
# fasta is the input file in .fasta format. 
summary.seqs(fasta=Reads_R1.trim.contigs.fasta)

#screen.seqs command enables you to keep sequences that fulfill certain user defined criteria.

screen.seqs(fasta=Reads_R1.trim.contigs.fasta, group=Reads_R1.contigs.groups, maxambig=0, maxlength=300)

unique.seqs(fasta=Reads_R1.trim.contigs.good.fasta)

count.seqs(name=current, group=current)

pcr.seqs(fasta=silva.bacteria.fasta, start=11894, end=25319, keepdots=F, processors=4)

system(mv silva.bacteria.pcr.fasta silva.v4.fasta)
align.seqs(fasta=Reads_R1.trim.contigs.good.unique.fasta, reference=silva.v4.fasta, flip=t)
summary.seqs(fasta=current,count=current)
align.seqs(fasta=Reads_R1.trim.contigs.good.unique.fasta, reference=silva.v4.fasta)
summary.seqs(fasta=current,count=current)
screen.seqs(fasta=Reads_R1.trim.contigs.good.unique.align, count=Reads_R1.trim.contigs.good.count_table, summary=Reads_R1.trim.contigs.good.unique.summary, start=1968, end=11550, maxhomop=8)
summary.seqs(fasta=current,count=current)
filter.seqs(fasta=Reads_R1.trim.contigs.good.unique.good.align, vertical=T, trump=.)
unique.seqs(fasta=Reads_R1.trim.contigs.good.unique.good.filter.fasta, count=Reads_R1.trim.contigs.good.good.count_table)
pre.cluster(fasta=Reads_R1.trim.contigs.good.unique.good.filter.unique.fasta, count=Reads_R1.trim.contigs.good.unique.good.filter.count_table, diffs=2)
chimera.uchime(fasta=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.fasta, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.count_table, dereplicate=t)
#dereplicate set to TRUE bc default is FALSE and Schloss SOP thinks this is too aggressive
remove.seqs(fasta=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.fasta, accnos=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.accnos)
summary.seqs(fasta=current,count=current)
classify.seqs(fasta=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, reference=trainset9_032012.pds.fasta, taxonomy=trainset9_032012.pds.tax, cutoff=80)
system(cp Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy silva.wang.taxonomy)



#this is getting rid of Eukaryota and Archaea
remove.lineage(fasta=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, taxonomy=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
summary.seqs(fasta=current,count=current)
#Preparing inputs for analysis
dist.seqs(fasta=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, processors=4,cutoff=0.2)
cluster(column=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.dist, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table)

#OTU-based files
make.shared(list=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, label=0.03)
classify.otu(list=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, taxonomy=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=0.03)



#Phylotype-based files
phylotype(taxonomy=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy)
make.shared(list=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, label=5)
system(cp Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.shared phylotype.5.phylum.shared)
classify.otu(list=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, taxonomy=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=5)
remove.groups(shared=phylotype.5.phylum.shared, groups=Coc_T4_491-Coc_T5_447-Coc_T5_546-Coc_T6_520-Coc_T8_435-Coc_TG57_Z241_1-Coc_TG57_Z241_2-Coc_TP1_Z254-Coc_Trip11_383_1-Coc_Trip11_383_2-Coc_Trip12_390_1-Coc_Trip12_390_2-Coc_Trip16_387-Coc_Trip7_392)

#Order level:
make.shared(list=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, label=3)
system(cp Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.shared phylotype.3.shared)
classify.otu(list=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, taxonomy=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=3)
remove.groups(shared=phylotype.3.shared, groups=Coc_T4_491-Coc_T5_447-Coc_T5_546-Coc_T6_520-Coc_T8_435-Coc_TG57_Z241_1-Coc_TG57_Z241_2-Coc_TP1_Z254-Coc_Trip11_383_1-Coc_Trip11_383_2-Coc_Trip12_390_1-Coc_Trip12_390_2-Coc_Trip16_387-Coc_Trip7_392)

#class level
make.shared(list=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, label=4)
system(cp Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.shared MiSeq.phylotype.class.shared)
classify.otu(list=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, taxonomy=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=4)

#FAMILY LEVEL SHARED FILES
make.shared(list=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, label=2)
system(cp Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.shared MiSeq_familyOTUs.shared)
classify.otu(list=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, taxonomy=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=2)

#genus level 
make.shared(list=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, label=1)
system(cp Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.shared MiSeq_genusOTUs.shared)
classify.otu(list=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.tx.list, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, taxonomy=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.pick.taxonomy, label=1)

#OTU shared file
system(cp Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared OTU03.an.shared)
system(cp Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy OTU03.an.cons.taxonomy)
count.groups(shared=OTU03.an.shared)

#for alpha diversity, getting rid of nonCrypt scales
remove.groups(shared=OTU03.an.shared, groups=Coc_T4_491-Coc_T5_447-Coc_T5_546-Coc_T6_520-Coc_T8_435-Coc_TG57_Z241_1-Coc_TG57_Z241_2-Coc_TP1_Z254-Coc_Trip11_383_1-Coc_Trip11_383_2-Coc_Trip12_390_1-Coc_Trip12_390_2-Coc_Trip16_387-Coc_Trip7_392)
remove.groups(shared=MiSeqOTU03.noSingle.shared, groups=Coc_T4_491-Coc_T5_447-Coc_T5_546-Coc_T6_520-Coc_T8_435-Coc_TG57_Z241_1-Coc_TG57_Z241_2-Coc_TP1_Z254-Coc_Trip11_383_1-Coc_Trip11_383_2-Coc_Trip12_390_1-Coc_Trip12_390_2-Coc_Trip16_387-Coc_Trip7_392)

#8/28/16 alpha diversity with non-subsampled shared file
summary.single(shared=OTU03.an.0.03.pick.shared, calc=nseqs-coverage-sobs-invsimpson, subsample=1100)
rarefaction.single(shared=OTU03.an.0.03.pick.shared, calc=sobs, freq=100)
rarefaction.single(shared=Bac16S_Crypt_only_OTUs_groupedbysampletype.shared, calc=sobs, freq=10)
rarefaction.single(shared=Bac16S_Crypt_only_nosingleOTUs_groupedbysampletype.shared, calc=sobs, freq=10)

#subsample at 1100, excludes 13 samples including all leaves
sub.sample(shared=OTU03.an.shared, size=1100)
remove.groups(shared=OTU03.an.0.03.subsample.shared, groups=Coc_T4_491-Coc_T5_447-Coc_T5_546-Coc_T6_520-Coc_T8_435-Coc_TG57_Z241_1-Coc_TG57_Z241_2-Coc_TP1_Z254-Coc_Trip11_383_1-Coc_Trip11_383_2-Coc_Trip12_390_1-Coc_Trip12_390_2-Coc_Trip16_387-Coc_Trip7_392)


##OTU-level venn diagrams 
merge.groups(shared=OTU03.an.0.03.subsample.shared,design=MiSeq_familyOTUs.2.subsample.design.txt)
##Azteca
venn(shared=OTU03.an.0.03.subsample.merge.shared,groups=Azteca-Azteca_dom-Scale_Azteca)
remove.groups(shared=OTU03.an.0.03.subsample.merge.shared, groups=Ceph_dom-Cephalotes-Scale_2007-Scale_Ceph-Scale_Chamela-Scale_Triplaris)
system(cp OTU03.an.0.03.subsample.merge.0.03.pick.shared AztecaOTU_only_subsample.shared)
get.sharedseqs(shared=AztecaOTU_only_subsample.shared, sharedgroups=Azteca_dom-Azteca-Scale_Azteca, output=accnos, label=0.03)#AztecaOTU_only_subsample.0.03.Azteca_dom-Azteca-Scale_Azteca.accnos
get.sharedseqs(shared=AztecaOTU_only_subsample.shared, uniquegroups=Azteca_dom-Azteca, output=accnos, label=0.03)
##Cephalotes
venn(shared=OTU03.an.0.03.subsample.merge.shared,groups=Cephalotes-Ceph_dom-Scale_Ceph)
remove.groups(shared=OTU03.an.0.03.subsample.merge.shared, groups=Azteca-Azteca_dom-Scale_Azteca-Scale_2007-Scale_Chamela-Scale_Triplaris)
system(cp OTU03.an.0.03.subsample.merge.0.03.pick.shared CephOTU_only_subsample.shared)
get.sharedseqs(shared=CephOTU_only_subsample.shared, sharedgroups=Cephalotes-Ceph_dom-Scale_Ceph, output=accnos, label=0.03)
get.sharedseqs(shared=CephOTU_only_subsample.shared, uniquegroups=Ceph_dom-Scale_Ceph, output=accnos, label=0.03)

#filter out singleton OTUs
split.abund(fasta=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta, list=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.list, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table, cutoff=1, label=0.03)
make.shared(list=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.abund.list, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.0.03.abund.count_table, label=0.03)
system(cp Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.abund.shared MiSeqOTU03.noSingle.shared)
count.groups(shared=MiSeqOTU03.noSingle.shared)
sub.sample(shared=MiSeqOTU03.noSingle.shared, size=1100)



#run classify.seqs again with greengenes for PICRUSt
classify.seqs(fasta=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, reference=gg_13_8_99.fasta, taxonomy=gg_13_8_99.gg.tax)
remove.lineage(fasta=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.fasta, count=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.count_table, taxonomy=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.gg.wang.taxonomy, taxon=Chloroplast-Mitochondria-unknown-Archaea-Eukaryota)
system(cp Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.gg.wang.pick.taxonomy greengenes.pick.taxonomy)
system(cp Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.fasta greengenes.pick.fasta)
system(cp Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.denovo.uchime.pick.pick.count_table greengenes.pick.count_table)
dist.seqs(fasta=greengenes.pick.fasta, processors=4,cutoff=0.2)
cluster(column=greengenes.pick.dist, count=greengenes.pick.count_table)
classify.otu(list=greengenes.pick.an.unique_list.list, count=greengenes.pick.count_table, taxonomy=Reads_R1.trim.contigs.good.unique.good.filter.unique.precluster.pick.gg.wang.pick.taxonomy, label=0.03)
make.shared(list=greengenes.pick.an.unique_list.list, count=greengenes.pick.count_table, label=0.03)
count.groups(shared=greengenes.coccoids.shared)
sub.sample(shared=greengenes.coccoids.shared, size=1100)
make.biom(shared=greengenes.coccoids.0.03.subsample.shared, label=0.03, reftaxonomy=gg_13_8_99.gg.tax, constaxonomy=greengenes.pick.an.unique_list.0.03.cons.taxonomy, picrust=97.gg.otu_map)
make.biom(shared=greengenes.coccoids.0.03.subsample_onlyCrypt.shared, label=0.03, reftaxonomy=gg_13_8_99.gg.tax, constaxonomy=greengenes.pick.an.unique_list.0.03.cons.taxonomy, picrust=97.gg.otu_map)




