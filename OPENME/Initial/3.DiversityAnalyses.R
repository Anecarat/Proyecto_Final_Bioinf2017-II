## Author: Víctor Taracena
## email: vtaracena@gmail.com

## Based on the work of Maggie Wagner and collaborators and the bioconductor
## workflow by Callaghan et al. (more info in 0.MasterScript.R)

                        ## 3.DiversityAnalyses ##

## In this one you will create the map of the site, calculate and plot different alpha diversity indexes, 
## rarefaction curves, Venn diagrams of the samples, relative abundance of the more abundant phyla and
## density plots of unique leaves and roots OTUs vs. shared ones. 

#### Stablishing colour palettes ####
sitePalette<-c("#0298f5","#7c30b6")
typePalette<-c("#36c412","#673604")
#### Site Map ####
## Load site coordinates
coords<-as.data.frame(read.table("RawData/SiteCoords.txt",sep="\t",header=TRUE))

## Make inset
png("Plots/Map.png")
par(mfrow=c(1,2))
map(database="state",regions=c("Montana","Idaho", "Wyoming","Oregon","Washington"),
    xlim=c(-117.25,-111.05),ylim=c(42,49),fill=TRUE,col="white")
map(database="state",regions=c("Idaho"),add=TRUE,fill=TRUE,col="gray90")
polygon(x=c(-113,-113,-115,-115),y=c(43.7,46,46,43.7),border="black",lwd=2)
text(x=c(-112.3), y=c(48.7),labels="Idaho", font=2, col="black",family="Helvetica", cex=2)
     
## Make regional map
map(database="county",regions=c("Idaho"),xlim=c(-115,-113),ylim=c(43.7,46))
polygon(x=c(-113,-113,-115,-115),y=c(43.7,46,46,43.7),border="black",lwd=8)
points(x=subset(coords,Site=="Mah")$Lon, y=subset(coords,Site=="Mah")$Lat, pch=17,cex=3,col=sitePalette[1])
points(x=subset(coords,Site=="Sil")$Lon, y=subset(coords,Site=="Sil")$Lat, pch=17,cex=3,col=sitePalette[2])

text(x=c(-113.8,-114.6, -113.4),
     y=c(44.05,45.02, 45.9),
     labels=c("Mah","Sil","Idaho"),
     col=c(sitePalette, "black"),cex=2,font=2,family="Helvetica")
dev.off()

#### Alpha diversity plots ####
png("Plots/ADiversity.png", 960, 480)
p<-plot_richness(bacteria, x="Site", color="Type")+
  scale_color_manual(values=typePalette)+
  theme_grey()+
  theme(axis.title = element_text(size=24))
p + geom_boxplot(data = p$data, aes(x = Site, color = Type), alpha=0.1)
dev.off()

#### Rarefaction Curves ####

## Note: Please be sure that the function calculate_rarefaction_curves is loaded
## if not load it through the 0.MasterScript.R

## Calculating alpha diversity this could take a while, 
## so take a break and go for some coffee :)
rarecurvedata<-calculate_rarefaction_curves(bacteria, 
               c('Observed'), rep(c(1, 10, 100, 1000, 1:100 * 10000),
                                             each = 10))
## Summarizing alpha diversity
rarecurvedatasum <- ddply(rarecurvedata, 
                          c('Depth', 'Sample', 'Measure'), 
                          summarise, Alpha_diversity_mean = mean(Alpha_diversity), 
                          Alpha_diversity_sd = sd(Alpha_diversity))
## Add sample data
rarecurvedatasum_verbose <- merge(rarecurvedatasum, data.frame(sample_data(bacteria)),
                                  by.x = 'Sample', by.y = 'row.names')
## Plot Rarefaction curves
png("Plots/Rarecurve.png", 600, 600)
ggplot(
  data = rarecurvedatasum_verbose,
  mapping = aes(x = Depth, y = Alpha_diversity_mean, ymin = Alpha_diversity_mean - Alpha_diversity_sd,
    ymax = Alpha_diversity_mean + Alpha_diversity_sd, colour = Type, group = Sample))+ 
  geom_line(size=1) + 
  ylab("Species") +
  xlab("Sample Size")+
  theme(axis.title = element_text(size=24))+
  theme_gray()+
  scale_color_manual(values=typePalette)+
  scale_x_continuous(labels = scales::comma)+
  facet_wrap(facets = ~ Site, scales = 'free_y')+
  theme(axis.title.x=element_text(size=22,face='bold'))+
  theme(axis.title.y=element_text(size=22,face='bold'))
dev.off()

#### Venn Diagrams ####
## Subsetting samples by type
roots<-subset_samples(bacteria,Type=='root' & Site%in%c('Mah','Sil')) %>%
  prune_taxa(taxa_sums(.)>0,.) # remove OTUs with 0 observations
leaves<-subset_samples(bacteria,Type=='leaf' & Site%in%c('Mah','Sil')) %>%
  prune_taxa(taxa_sums(.)>0,.) # remove OTUs with 0 observations

## Core across sites
## Leaves
sitePalette<-c("#0298f5","#7c30b6")
png("Plots/VennDLeaves.png")
grid.newpage()
VDleaves<-venn.diagram(x=list(Mah=taxa_names(subset_samples(leaves,Site=='Mah')%>%prune_taxa(taxa_sums(.)>0,.)),
                    Sil=taxa_names(subset_samples(leaves,Site=='Sil')%>%prune_taxa(taxa_sums(.)>0,.))),
             cat.col=sitePalette[1:2],fill=sitePalette[1:2],alpha=c(0.5),cat.cex=2,cat.fontface="bold",cat.fontfamily="Helvetica",
             label.col="white",cex=1.7,fontface="bold",fontfamily="Helvetica",
             main="leaf-associated OTUs",main.col="forest green",main.cex="2",main.fontface="bold",main.fontfamily="Helvetica",
             filename=NULL) %>% grid.draw()
dev.off()
## Roots
png("Plots/VennDRoots.png")
grid.newpage()
VDroots<-venn.diagram(x=list(Mah=taxa_names(subset_samples(roots,Site=='Mah')%>%prune_taxa(taxa_sums(.)>0,.)),
                    Sil=taxa_names(subset_samples(roots,Site=='Sil')%>%prune_taxa(taxa_sums(.)>0,.))),
             cat.col=sitePalette[1:2],fill=sitePalette[1:2],alpha=c(0.5),cat.cex=2,cat.fontface="bold",cat.fontfamily="Helvetica",
             label.col="white",cex=1.7,fontface="bold",fontfamily="Helvetica",
             main="root-associated OTUs",main.col="grey",main.cex="2",main.fontface="bold",main.fontfamily="Helvetica",
             filename=NULL) %>% grid.draw()
dev.off()
#### Relative abundance of major phyla ####
## Agglomerate taxa by phylum
bac.phy<-tax_glom(bacteria, taxrank="Phylum")
taxa_names(bac.phy)<-as.data.frame(tax_table(bac.phy))$Phylum
## Find top 12 phyla by abundance
sort(taxa_sums(bac.phy), decreasing=TRUE)[1:12] %>% names() -> top12phyla
setdiff(taxa_names(bac.phy),top12phyla) -> rarephyla
## Consolidate rare phyla into "Low abundance" category
merged.bac.phy<-merge_taxa(bac.phy, rarephyla, "TM7") %>%
  prune_taxa(taxa_sums(.)>0,.)
ntaxa(merged.bac.phy) # 13 taxa remain, as expected

## Extract phylum counts & sample metadata into a data frame
phydata<-as(otu_table(merged.bac.phy),'matrix') %>% t() %>% as.data.frame() %>% 
  mutate('SampleID'=factor(row.names(.))) %>% 
  merge(.,as(sample_data(merged.bac.phy),'data.frame')%>%
          select(SampleID,Type,Site),by='SampleID') %>% # merge with sample metadata
  gather(key=Phylum,value=Abundance,2:14) %>% # melt data frame
  group_by(SampleID) %>% mutate(RelAbund=Abundance/sum(Abundance)) %>%   # transform into relative abundance in each sample:
  ungroup %>% as.data.frame %>%
  mutate(Phylum=ifelse(Phylum=='TM7','Low abundance',Phylum)) %>% # rename TM7 as 'Low abundance'
  mutate(Phylum=ordered(Phylum,levels=rev(c('Low abundance',top12phyla[12:1]))))

phyPalette<-c("#000000","#B15928","#33A02C","#B2DF8A","#1F78B4","#FB9A99",
              "#FDBF6F","#E31A1C","#6A3D9A","#CAB2D6","#FF7F00","#FFFF99",
              "#A6CEE3")


png("Plots/AbunRel.png", 960,480)
  ggplot(phydata,aes(x=SampleID, y=100*RelAbund, fill=Phylum, color=Phylum)) +
  geom_bar(stat='identity') +
  facet_wrap(Type~Site,scales='free_x',strip.position='bottom')+
  theme_classic()+ 
  ylab("Relative abundance (%)")+xlab("Site")+
  scale_fill_manual(values=phyPalette[13:1])+
  scale_color_manual(values=phyPalette[13:1],guide=FALSE)+
  scale_x_discrete(breaks=NULL)+
  theme(strip.background=element_blank(),strip.text.x=element_text(size=26,face='bold',vjust=1.5))+
  theme(axis.title.x=element_text(size=22,face='bold'))+
  theme(axis.title.y=element_text(size=22,face='bold'),axis.text.y=element_text(size=24,face='bold'))
  dev.off()

#### Unique OTUs to leaves and roots ####
leafonlyOTUs<-setdiff(taxa_names(leaves),taxa_names(roots)) 
rootonlyOTUs<-setdiff(taxa_names(roots),taxa_names(leaves))
  
RelAbund.otu<-data.frame("Taxon"=taxa_names(bacteria),"Level"='otu') %>%
mutate("leaves"=taxa_sums(subset_samples(bacteria,Site %in% c('Mah','Sil') &  Type=='leaf'))/sum(taxa_sums(subset_samples(bacteria,Site %in% c('Mah','Sil') &  Type=='leaf')))) %>%
mutate("roots"=taxa_sums(subset_samples(bacteria,Site %in% c('Mah','Sil') &  Type=='root'))/sum(taxa_sums(subset_samples(bacteria,Site %in% c('Mah','Sil') &  Type=='root'))))

## How common are these?
filter(RelAbund.otu, Taxon%in%leafonlyOTUs) %>% select(leaves) %>% sum 
filter(RelAbund.otu, Taxon%in%rootonlyOTUs) %>% select(roots) %>% sum 
  
png("Plots/leafonlyOTU_Abun.png", 700,700)
  ggplot(RelAbund.otu,aes(x=log(leaves)))+
    geom_density(aes(fill="all leaf OTUs\n",colour="all leaf OTUs\n"),alpha=0.3,position="identity")+
    geom_density(data=subset(RelAbund.otu,Taxon%in%leafonlyOTUs),aes(x=log(leaves),colour="OTUs exclusively\nfound in leaves",fill="OTUs exclusively\nfound in leaves"),alpha=0.5,position="identity")+
    scale_x_continuous(labels=function(x) format(exp(x),digits=4),breaks=c(log(0.0000001),log(0.00001),log(0.001),log(0.1)))+
    scale_fill_manual(values=c("#1b6209", typePalette[1]),name="")+
    scale_colour_manual(values=c("#1b6209", typePalette[1]),name="")+
    theme_classic()+xlab("Relative abundance")+theme(legend.position=c(0.8,0.9))+
    theme(legend.title = element_text(size=20, face="bold"), legend.text = element_text(size=28,face="bold"))+
    theme(axis.title.x = element_text(size=36,face="bold"), axis.text.x = element_text(size=30,face="bold"))+
    theme(axis.title.y = element_text(size=36,face="bold"), axis.text.y = element_text(size=30,face="bold"))
  dev.off()
  
png("Plots/rootonlyOTU_Abund.png",700,700)
  ggplot(RelAbund.otu,aes(x=log(roots)))+
    geom_density(aes(fill="all root OTUs\n",colour="all root OTUs\n"),alpha=0.3,position="identity")+
    geom_density(data=subset(RelAbund.otu,Taxon%in%rootonlyOTUs),aes(x=log(roots),colour="OTUs exclusively\nfound in roots",fill="OTUs exclusively\nfound in roots"),alpha=0.5,position="identity")+
    scale_x_continuous(labels=function(x) format(exp(x),digits=4),breaks=c(log(0.0000001),log(0.00001),log(0.001),log(0.1)))+
    scale_fill_manual(values=c("#1e1001", typePalette[2]),name="")+
    scale_colour_manual(values=c("#1e1001", typePalette[2]),name="")+
    theme_classic()+xlab("Relative abundance")+theme(legend.position=c(0.8,0.9))+
    theme(legend.title = element_text(size=20, face="bold"), legend.text = element_text(size=28,face="bold"))+
    theme(axis.title.x = element_text(size=36,face="bold"), axis.text.x = element_text(size=30,face="bold"))+
    theme(axis.title.y = element_text(size=36,face="bold"), axis.text.y = element_text(size=30,face="bold"))
  dev.off()
  
  
  
          ## Please, continue to 4.MultivariateAnalyses.R script ## 