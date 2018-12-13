## Rhoades Final project

Main goals 

1.Generate taxonomic abundance shotgun metagenomic data.<br />
2. Plot this data using principal cordinate analysis.<br />
3.Look at the 10 most abundant bacterial species using a heatmap. <br />

##Project info
This data is from a cohort of infant rhesus macaques that had chronic diarrhea and had to be Euthanized. We have shotgun metagenomic data from colonic luminal contents along with some healthy controls


Also important I am in Dr. Messaoudis' lab and we actully work on the UC Riverside cluster so some of the intial sequence processing steps will likely not work on the UCI cluster. To rememedy this I have made processed data along the way downloadable.

These libraries are each about 20-30 million 100bp sequences each so I reduced them to 100 thousand  sequences so this all runs a bit faster. 

The data table I indroduce later is from the full dataset.<br />

##Setting up the project folder


```
# Get to a directory where you can start a new project

cd /pub/jje/ee282/$USER

mkdir Final_project_NR
cd Final_project_NR
mkdir trim decon species
cp /pub/jje/ee282/rhoadesn/Final_project_NR/results ./
cp /pub/jje/ee282/rhoadesn/Final_project_NR/sequences ./
cp /pub/jje/ee282/rhoadesn/Final_project_NR/Scripts ./

cd /pub/jje/ee282/rhoadesn/Final_project/sequences
```

## How I reduced down to 100k reads.

```
#module load seqtk
#seqtk sample -s100 Sample1.fq.gz 100000 > Sample1_100k.fq  
```

#Example Scripts Do not run
I wrote out how to run each step on a single sample, it would be inefficent to run these all one at a time for 16 samples. The reall scripts are in the immplmentations section.

##Trimming low quality bases using trimmomatic. <br />
___We will actutally implment this using a perl loop, but below is a single sample example___


###Sequence QC using trimmomatic

module load trimmomatic

###Example single sample script (.jar file probaly lacking on UCI server)
java -j /path_to_jar_file/trimmomatic-0.36.jar SE -threads 1 -phred33 Sample1_100k.fq.gz ../trim/Sample1.qc.fq.gz ILLUMINACLIP:NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

###Remove host sequence. Since these metagenomes are from luminal contents we need to remove any read that map to the host genome. Using Bowtie2.

___Again there is a perl loop further down that we will actually implment___


Change directory to the newly trimmed sequence

cd ../trim

###Example Single sample Script

module load bowtie2

bowtie2 -q -p 1 -x /path_to_bowtie_index/ -U Sample1.qc.fq.gz --quiet --un ../decon/Sample1.decon.qc.fq.gz 2> ../decon alignment_stats.txt 



###Asigning taxonomy using MetaPhlan2.
This assigns taxonomy to short reads.

___Again there is a perl loop further down that we will actually implment___


Change directory to the newly host decontaminated sequences

cd ../decon

module load metaphlan2

metaphlan2.py Sample1.decon.qc.fq.gz --input_type fastq --nproc 1 --tax_lev s --ignore_viruses --sample_id Sample1 -o ../species/SampleID.species.txt


#Implementation (run these)


This all runs pretty  fast on the scalled down data set even with 1 core


###Trimming

```
qrsh -q class
module load perl
module load trimmomatic
module load bowtie2
module load metaphlan2

cd /pub/jje/ee282/$USER/Final_project_NR/sequences
perl ../Scripts/run_trimmomatic_SE.pl
```
###Decontamination<br />
Had to nohup the output because the script was flooding my terminal with dialog. (sorry)

```
cd /pub/jje/ee282/$USER/Final_project_NR/trim
nohup perl ../Scripts/run_bowtie2NR.pl
rm nohup.out

```
###Speices identification
## Metaphlan module is broken on UCI Cluster

```
#This is what you would normally run
cd /pub/jje/ee282/$USER/Final_project_NR/decon
perl run_metaphlan2.pl

#Alternitive just copy the data i generated on the UCR cluster
cd /pub/jje/ee282/$USER/Final_project_NR
cp /pub/pub/jje/ee282/rhoadesn/Final_project_NR/species
```
###Merge Data
Merge into a single table
and remove the 1st line of the table that we dont need

```
cd /pub/jje/ee282/$USER/Final_project_NR/species
../Scripts/merge_metaphlan_tables.py *.species.txt | sed '1d' > ../results/Colon_species.txt

head ../results/Colon_species.txt


```
#Part 2 (plotting)
###Next we will plot the species level data in R this is using the output from the full dataset not the 100k version.
###First we generate a principal component analysis

```
cd /pub/jje/ee282/$USER/Final_project_NR/results

module load R/3.5.1

R


library(vegan)
library(ggplot2)
library(grid)

#Read in our data
inTable = read.delim("Species_colon_full.txt", header=TRUE, row.names=1, sep='\t', check.names=FALSE, comment.char="")

#This map is a 2 column text file that has the health status of each animal
map = read.table("map.txt", header=T, sep='\t', row.names=1, check.names=F, comment.char='')

#generates a dissimilarity matrix from our input table
distMat <- vegdist(inTable, method="bray", binary=FALSE, diag=TRUE, upper=TRUE, na.rm = FALSE )

#The actual principal coordinate analysis
pcoa_cord = data.frame(cmdscale(distMat,k=10))
names(pcoa_cord) = c('PCoA1', 'PCoA2','PCoA3')
pa = cmdscale(distMat, eig = TRUE)

#using a factor from our map to color the plot
factorColor = map$Status

##makes PCoA plot##
png(filename="Colon_species_pcoa.png", 
    units="in", 
    width=8, 
    height=8, 
    pointsize=10, 
    res=300)
 ##GGplot has a ton of options, the ones listed below just serve to change the apperance of the plot##    
ggplot(pcoa_cord,aes(x=PCoA1, y=PCoA2, color=factorColor)) +
  
  geom_point(na.rm=T, alpha=0.8, size=8) +
  #changes theme to black and white#
  theme_bw()+
  labs(title ="Colon species", x = "PCoA1", y = "PCoA2 ", legend="aa") +
  scale_color_manual(name='Status', 
                     breaks=c('Sick', 'Healthy'),values=c('green','red')) +
  
  #removes tick marks on axis, can toggle on and off#
  theme(axis.ticks=element_blank()) +
  theme(axis.text=element_text(size = 20)) +
 
  
  #changes appearance of plot border#
  theme(panel.border=element_rect(colour="black", size=1.1)) +
  theme(axis.title= element_text(size =18)) +
  
  ##changes apperance of text in legend##
  theme(legend.text = element_text(size = 12))+

  ##changes appearnce of major grid lines in plot##
 theme(panel.grid.major=element_line(linetype="dashed", colour="white", size=0.0)) +
  
  ##changes appearance of minor grid lines in plot##
 theme(panel.grid.minor=element_line(linetype="dashed", colour="white", size=0.0)) +
  
  ##changes legend appearance for shape element
  scale_shape_discrete(name='Location')     
  dev.off()

```
## The next script will make a heatmap of the 10 most abundant bacterial species still in R

```
args <- commandArgs(TRUE)
library("RColorBrewer")
library(gplots)
library(lattice)

#Define the colors
hmcol <- colorRampPalette(brewer.pal(9, "RdBu"))(100)

# file1 is the list of merged genes
file1  <- args[1]
# Output name - heatmap pdf prefix
outprefix <- args[2]

#Open file handles and save the first columns
f1 <- read.delim("Species_colon_full.txt", sep="\t", header=TRUE, row.names=1)


f2 <- colMeans(f1, na.rm=TRUE)
f3 <- sort(f2, decreasing = FALSE)
f1 <- f1[,order(f3)]
f4 <- subset(f1, select=c(1,2,3,4,5,6,7,8,9,10))
y <- (t(as.matrix(f4)))
```
## plotting untransformed heatmap
```
png(filename="Colon_species_HM.png", 
    units="in", 
    width=8, 
    height=10, 
    pointsize=10, 
    res=300)
levelplot(t(y), height=0.3, col.regions=rev(hmcol), scales=list(x=list(rot=90)), main="", colorkey=list(space="top"), xlab="", ylab="", pretty=TRUE, width=0.75, cexRow=0.1, cexCol=0.1, aspect=2.5)
dev.off()

```
## log transformed heatmap

```
f4 <- log10(f4+1)

y <- (t(as.matrix(f4)))

png(filename="Colon_species_log_HM.png", 
    units="in", 
    width=8, 
    height=10, 
    pointsize=10, 
    res=300)
levelplot(t(y), height=0.3, col.regions=rev(hmcol), scales=list(x=list(rot=90)),main="", colorkey=list(space="top"), xlab="", ylab="", pretty=TRUE, width=0.5, cexRow=0.1, cexCol=0.1, aspect=2.5)
dev.off()
```
