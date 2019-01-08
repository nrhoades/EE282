# EE282 HW4

**To download the all chromosomes fasta file into a new folder:**

```
mkdir melanogaster_data
cd melanogaster_data
wget ftp://ftp.flybase.net/genomes/Drosophila_melanogaster/current/fasta/dmel-all-chromosome-r6.24.fasta.gz
```

**To calculate total number of nucleotides, total number of Ns, and total number of all sequences ≤ 100kb and all sequences > 100kb:**

```
module load jje/jjeutils
module load jje/kent

bioawk -c fastx 'length($seq) > 100000{ print ">"$name; print $seq }' dmel-all-chromosome-r6.24.fasta.gz > 100kb.fa

bioawk -c fastx 'length($seq) <= 100000{ print ">"$name; print $seq }' dmel-all-chromosome-r6.24.fasta.gz > 99kb.fa


faSize 100kb.fa > summary_100kb.txt
faSize 99kb.fa > summary_99kb.txt

less summary_100kb.txt
less summary_99kb.txt

bioawk -c fastx ' { print length($seq) }' dmel-all-chromosome-r6.24.fasta.gz \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nWG_Ctg\t0" } { print "WG_Ctg\t" $1 } ' \
> Whole_genome_sort.txt

bioawk -c fastx ' { print length($seq) }' 100kb.fa \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\n100kb_Ctg\t0" } { print "100kb_Ctg\t" $1 } ' \
> 100kb_sort.txt


bioawk -c fastx ' { print length($seq) }' 99kb.fa  \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\n99kb_Ctg\t0" } { print "99kb_Ctg\t" $1 } ' \
> 99kb_sort.txt

plotCDF2 {Whole_genome,100kb,99kb}_sort.txt /dev/stdout \
| tee WG_v_100_99.png \
| display
```
#Combined CDF#
![CDF](https://github.com/nrhoades/EE282/blob/master/WG_v_100_99.png)

```
bioawk -c fastx ' { print gc($seq) }' dmel-all-chromosome-r6.24.fasta | sort -rn | awk ' BEGIN { print "Assembly\tGC" } { print "WG_Ctg\t" $1 } ' > Whole_genome_gc_sort.txt

bioawk -c fastx ' { print gc($seq) }' 100kb.fa | sort -rn | awk ' BEGIN { print "Assembly\tGC" } { print "WG_Ctg\t" $1 } ' > 100kb_gc_sort.txt

bioawk -c fastx ' { print gc($seq) }' 99kb.fa | sort -rn | awk ' BEGIN { print "Assembly\tGC" } { print "WG_Ctg\t" $1 } ' > 99kb_gc_sort.txt

module load R
R
kb99=read.table("99kb_sort.txt", header=T, dec=".", sep="\t")
kb100=read.table("100kb_sort.txt", header=T, dec=".", sep="\t")
WG=read.table("Whole_genome_sort.txt", header=T, dec=".", sep="\t")
kb99_gc=read.table("99kb_gc_sort.txt", header=T, dec=".", sep="\t")
kb100_gc=read.table("100kb_gc_sort.txt", header=T, dec=".", sep="\t")
WG_gc=read.table("Whole_genome_gc_sort.txt", header=T, dec=".", sep="\t")
library(ggplot2)

kb99_plot <- ggplot(data = kb99)
kb99_plot + geom_histogram(mapping = aes(x = Length), bins=10)
pdf("kb99_size.pdf")
dev.off()

kb100_plot <- ggplot(data = kb100)
kb100_plot + geom_histogram(mapping = aes(x = Length), bins=10)
pdf("kb100_size.pdf")
dev.off()

WG_plot <- ggplot(data = WG)
WG_plot + geom_histogram(mapping = aes(x = Length), bins=10)
pdf("WG_size.pdf")
dev.off()

kb100_gc_plot <- ggplot(data = kb100_gc)
kb100_gc_plot + geom_histogram(mapping = aes(x = GC), bins=10)
pdf("kb100_gc.pdf")
dev.off()

kb99_gc_plot <- ggplot(data = WG_gc)
kb99_gc_plot + geom_histogram(mapping = aes(x = GC), bins=10)
pdf("kb99_gc.pdf")
dev.off()

WG_gc_plot <- ggplot(data = WG_gc)
WG_gc_plot + geom_histogram(mapping = aes(x = GC), bins=10)
pdf("WG_gc.pdf")
dev.off()

```

#Whole genome size distribution
![Whole_genome_Size](https://github.com/nrhoades/EE282/blob/master/WG_size.png)

#Greater than 100kb size distribution
![<100kb_Size](https://github.com/nrhoades/EE282/blob/master/kb100_size.png)

#Less than 100kb size distribution
![>100kb_Size](https://github.com/nrhoades/EE282/blob/master/kb99_size.png)

#Whole genome gc distribution
![Whole_genome_gc](https://github.com/nrhoades/EE282/blob/master/WG_gc.png)

#Greater than 100kb gc distribution
![<100kb_gc](https://github.com/nrhoades/EE282/blob/master/kb100_gc.png)

#Less than 100kb gc distribution
![>100kb_gc](https://github.com/nrhoades/EE282/blob/master/kb99_gc.png)

### Comments on "Summarize partitions of a genome assembly"

Good job. I know your results are faithfully recorded in the txt files redirected from the ```faSize``` utility, but in general, please record those answes in this document. Part of this class is encouraging good habits. Otherwise, great job.

**Assemble a genome from MinION reads**

```
module load jje/jjeutils

n50 () {
  bioawk -c fastx ' { print length($seq); n=n+length($seq); } END { print n; } ' $1 \
  | sort -rn \
  | gawk ' NR == 1 { n = $1 }; NR > 1 { ni = $1 + ni; } ni/n > 0.5 { print $1; exit; } '
}

minimap=$(which minimap)
miniasm=$(which miniasm)
basedir=/pub/jje/ee282/$USER
projname=nanopore_assembly
basedir=$basedir/$projname
raw=$basedir/$projname/data/raw
processed=$basedir/$projname/data/processed
figures=$basedir/$projname/output/figures
reports=$basedir/$projname/output/reports

createProject $projname $basedir
ln -sf /bio/share/solarese/hw4/rawdata/iso1_onp_a2_1kb.fastq $raw/reads.fq

$minimap -t 32 -Sw5 -L100 -m0 $raw/reads.fq{,} \
| gzip -1 \
> $processed/onp.paf.gz

$miniasm -f $raw/reads.fq $processed/onp.paf.gz \
> $processed/reads.gfa

awk ' $0 ~/^S/ { print ">" $2" \n" $3 } ' $processed/reads.gfa \
| tee >(n50 /dev/stdin > $reports/n50.txt) \
| fold -w 60 \
> $processed/unitigs.fa


cp unitig.fa ../../../../melanogaster_data

```

**Assembly assessment**
_N50 was calculated above, to see results Reference N50= 21,485,538 Our N50= 4,494,246_

```
less ./nanopore_assembly/nanopore_assembly/output/reports/n50.txt

```

**Compare your assembly to the contig assembly from Drosophila melanogaster on FlyBase using Mummer**

```
module load jje/jjeutils 
module load jje/kent
faSplitByN dmel-all-chromosome-r6.24.fasta  dmel-all-chromosome-r6.24.contigs.fasta  10

###Loading of binaries via module load or PATH reassignment
source /pub/jje/ee282/bin/.qmbashrc
module load gnuplot/4.6.0

###Query and Reference Assignment. State my prefix for output filenames
REF="dmel-all-chromosome-r6.24.contig.fasta"
PREFIX="flybase"
SGE_TASK_ID=1
QRY=$(ls u*.fa | head -n $SGE_TASK_ID | tail -n 1)
PREFIX=${PREFIX}_$(basename ${QRY} .fa)

###please use a value between 75-150 for -c. The value of 1000 is too strict.
nucmer -l 100 -c 1000 -d 10 -banded -D 5 -prefix ${PREFIX} ${REF} ${QRY}
mummerplot --fat --layout --filter -p ${PREFIX} ${PREFIX}.delta \
  -R ${REF} -Q ${QRY} --png

```
![Mummer_plot](https://github.com/nrhoades/EE282/blob/master/flybase_unitigs.png)

**Compare your assembly to both the contig assembly and the scaffold assembly from the Drosophila melanogaster on FlyBase using a contiguity plot**

```
module load jje/jjeutils 
module load jje/kent

bioawk -c fastx ' { print length($seq) }' dmel-all-chromosome-r6.24.fasta \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nScaff_Ctg\t0" } { print "Scaff_Ctg\t" $1 } ' \
> Scaffold_sort.txt

bioawk -c fastx ' { print length($seq) }' dmel-all-chromosome-r6.24.contig.fasta \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nContig_Ctg\t0" } { print "Contig_Ctg\t" $1 } ' \
> Contigs_sort.txt 


bioawk -c fastx ' { print length($seq) }' unitigs.fa \
| sort -rn \
| awk ' BEGIN { print "Assembly\tLength\nMine_Ctg\t0" } { print "Mine_Ctg\t" $1 } ' \
> My_assembly_sort.txt


plotCDF2 {Scaffold,Contigs,My_assembly}_sort.txt /dev/stdout \
| tee Scaf_Contig_Mine_CDF.png \
| display
```
![CDF2](https://github.com/nrhoades/EE282/blob/master/Scaf_Contig_Mine_CDF.png)


**Calculate BUSCO scores of both assemblies and compare them**

_Your working directory should contain the assembly fasta you want the BUSCO scores for._

```
module load augustus/3.2.1
module load blast/2.2.31 hmmer/3.1b2 boost/1.54.0
source /pub/jje/ee282/bin/.buscorc

INPUTTYPE="geno"
MYLIBDIR="/pub/jje/ee282/bin/busco/lineages/"
MYLIB="diptera_odb9"
OPTIONS="-l ${MYLIBDIR}${MYLIB}"
QRY="unitigs.fa"
MYEXT=".fasta" 


#you can change the value after -c to tell busco how many cores to run on. 
BUSCO.py -c 40 -i ${QRY} -m ${INPUTTYPE} -o $(basename ${QRY} ${MYEXT})_${MYLIB}${SPTAG} ${OPTIONS}

##run for contigs assembly

INPUTTYPE="geno"
MYLIBDIR="/pub/jje/ee282/bin/busco/lineages/"
MYLIB="diptera_odb9"
OPTIONS="-l ${MYLIBDIR}${MYLIB}"
QRY="dmel-all-chromosome-r6.24.contig.fasta"
MYEXT=".fasta" 


#you can change the value after -c to tell busco how many cores to run on. 
BUSCO.py -c 40 -i ${QRY} -m ${INPUTTYPE} -o $(basename ${QRY} ${MYEXT})_${MYLIB}${SPTAG} ${OPTIONS}


```
_Contig assembly Busco  C:98.3%[S:97.8%,D:0.5%],F:0.9%,M:0.8%,n:2799_


_My Busco C:0.5%[S:0.5%,D:0.0%],F:1.1%,M:98.4%,n:2799_
