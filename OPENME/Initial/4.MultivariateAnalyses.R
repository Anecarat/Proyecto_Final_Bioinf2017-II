## Author: Víctor Taracena
## email: vtaracena@gmail.com

## Based on the work of Maggie Wagner and collaborators and the bioconductor
## workflow by Callaghan et al. (more info in 0.MasterScript.R)

                     ## 4.MultivariateAnalyses ##
## Finally, this script will make a multidimensional scaling analysis 
## (a.k.a. principal coordinates analysis) to visualize the level of similarity 
## between the bacteria samples of leaves and roots between sites. 

#### Stablishing colour palettes ####
sitePalette<-c("#0298f5","#7c30b6")
typePalette<-c("#36c412","#673604")
#### Subsetting samples ####
root.vst<-subset_samples(bacteria,Type=='root' & Site%in%c('Mah','Sil'))
leaf.vst<-subset_samples(bacteria,Type=='leaf' & Site%in%c('Mah','Sil'))
root.vst<-subset_samples(bacteria,Type=='root' & Site%in%c('Mah','Sil')) %>%
  prune_taxa(taxa_sums(.)>0,.) # remove OTUs with 0 observations
leaf.vst<-subset_samples(bacteria,Type=='leaf' & Site%in%c('Mah','Sil')) %>%
  prune_taxa(taxa_sums(.)>0,.) # remove OTUs with 0 observations

#### Weighted UniFrac and PCoA: separately for leaf and root datasets ####

## replace negative values with 0s just for distance calculations
wUF.root.vst<-UniFrac(transform_sample_counts(root.vst,
                      function(x) x<-ifelse(x<0,0,x)),weighted=TRUE,
                      parallel=TRUE)
wUF.leaf.vst<-UniFrac(transform_sample_counts(leaf.vst,
                      function(x) x<-ifelse(x<0,0,x)),weighted=TRUE,
                      parallel=TRUE)

cap.wUF.root.vst<-capscale(wUF.root.vst~1,data=as(sample_data(root.vst),'data.frame'))
cap.wUF.leaf.vst<-capscale(wUF.leaf.vst~1,data=as(sample_data(leaf.vst),'data.frame'))

## Get inertia for top 3 PCoA axes: leaf ##
cap.wUF.leaf.vst$CA$eig[1:3]/sum(cap.wUF.leaf.vst$CA$eig) 
cap.wUF.root.vst$CA$eig[1:3]/sum(cap.wUF.root.vst$CA$eig) 

#### Scree plots: weighted UniFrac ####
png("Plots/ScreePlotsLeaves.png", 800, 800)
barplot(cap.wUF.leaf.vst$CA$eig[1:20]/sum(cap.wUF.leaf.vst$CA$eig),col=typePalette[1],
        main="Leaves: weighted UniFrac",xlab="Principal coordinate",ylab="Proportion Variance",
        cex.lab=2,cex.main=2.5,cex.axis=1.5,font.axis=2,font.main=2,font.lab=2,xaxt='n')
dev.off()
png("Plots/ScreePlotsRoots.png",800,800)
barplot(cap.wUF.root.vst$CA$eig[1:20]/sum(cap.wUF.root.vst$CA$eig),col=typePalette[2],
        main="Roots: weighted UniFrac",xlab="Principal coordinate",ylab="Proportion Variance",
        cex.lab=2,cex.main=2.5,cex.axis=1.5,font.axis=2,font.main=2,font.lab=2,xaxt='n')
dev.off()

#### How much variation is explained by the top 3 PCo axes? ####
sink("Data/ordination_top3_cumulative_PVEs.txt")
print("cumulative percent variance explained by top 3 PCo:")
print("weighted UniFrac, roots:")
sum(cap.wUF.root.vst$CA$eig[1:3])/sum(cap.wUF.root.vst$CA$eig)
print("weighted UniFrac, leaves:")
sum(cap.wUF.leaf.vst$CA$eig[1:3])/sum(cap.wUF.leaf.vst$CA$eig)
print("Individual percent variance explained by top 3 PCo:")
print("weighted UniFrac, roots:")
(cap.wUF.root.vst$CA$eig[1:3])/sum(cap.wUF.root.vst$CA$eig) 
print("weighted UniFrac, leaves:")
(cap.wUF.leaf.vst$CA$eig[1:3])/sum(cap.wUF.leaf.vst$CA$eig) 
sink()

#### Plotting PcoA #### 
png("Plots/PCoAleaves.png")
plot_ordination(leaf.vst,cap.wUF.leaf.vst,type="samples",axes=1:2,color="Site")+
  scale_colour_manual(values=sitePalette)+
  geom_point(size=4,alpha=1)+
  xlab("PCo1 [41.7%]")+ylab("PCo2 [26.4%]")+
  ggtitle("Leaves")+theme_classic()+
  theme(plot.title = element_text(size=44, face="bold"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=30,face="bold"))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=30,face="bold"))+
  theme(legend.title= element_text(size=40),legend.text=element_text(size=36,face="bold"))+
  theme(legend.key.height=unit(2.5,"lines"),legend.key.width=unit(2,"lines"))+
  theme(legend.background = element_rect(fill="gray90", size=.5))
dev.off()

png("Plots/PCoAroots.png")
plot_ordination(root.vst,cap.wUF.root.vst,type="samples",axes=1:2,color="Site")+
  scale_colour_manual(values=sitePalette,guide=FALSE)+
  geom_point(size=4,alpha=1)+
  xlab("PCo1 [34.1%]")+ylab("PCo2 [17.9%]")+
  ggtitle("Roots")+theme_classic()+
  theme(plot.title = element_text(size=44, face="bold"))+
  theme(axis.title.x=element_text(size=36,face="bold"),axis.text.x=element_text(size=30,face="bold"))+
  theme(axis.title.y=element_text(size=36,face="bold"),axis.text.y=element_text(size=30,face="bold"))
dev.off()

