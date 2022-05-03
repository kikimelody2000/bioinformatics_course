#!/bin/bash

cd ~
mkdir assessment
mkdir assessment/dnaseq
cd ~/assessment/dnaseq 
mkdir data results
cd ~/assessment/dnaseq/data
mkdir untrimmed_fastq
mkdir trimmed_fastq

#download the input files4
cd ~/assessment

wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R1.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/NGS0001.R2.fastq.qz
wget https://s3-eu-west-1.amazonaws.com/workshopdata2017/annotation.bed
wget http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz

# need to move the files and convert files to gz files 
mv ~/assessment/NGS0001.R1.fastq.qz ~/assessment/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq.gz
mv ~/assessment/NGS0001.R2.fastq.qz ~/assessment/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq.gz

#perform fastqc 
fastqc -t 4 *.fastq.gz

#perform trimming 
trimmomatic PE -threads 4 -phred33 ~/assessment/dnaseq/data/untrimmed_fastq/NGS0001.R1.fastq.gz ~/assessment/dnaseq/data/untrimmed_fastq/NGS0001.R2.fastq.gz  -baseout ~/assessment/dnaseq/data/trimmed_fastq/trimmed_data ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 TRAILING:25 MINLEN:50

fastqc -t 4 /home/ubuntu/ assessment/dnaseq/data/trimmed_fastq/trimmed_data_1P /home/ubuntu/assessment/dnaseq/data/trimmed_fastq/trimmed_data_2P

#make a new folder to put the fastqc trimmed reads output
mkdir ~/assessment/dnaseq/results/fastqc_trimmed_reads 

#move the trimmed fastqc outputs to the results fastqc trimmed reads folder
mv ~/assessment/dnaseq/data/trimmed_fastq/ trimmed_data_1P  ~/ assessment/dnaseq /results/fastqc_trimmed_reads/ 
mv ~/assessment/dnaseq/data/trimmed_fastq/ trimmed_data_2P  ~/ assessment/dnaseq /results/fastqc_trimmed_reads/ 

#make the BWA index
mkdir ~/assessment/dnaseq/data/reference
mv ~/assessment/dnaseq/data/reference/hg19.fa.gz ~/assessment/dnaseq/data/reference/
bwa index ~/assessment/dnaseq/data/reference/hg19.fa.gz 

#run BWA mem
mkdir ~/assessment/dnaseq/data/aligned_data

bwa mem -t 4 -v 1 -R '@RG\tID:11V6WR1:111:D1375ACXX:1.NGS0001\tSM:NGS0001\tPL:ILLUMIN\tLB:NGS0001-blood\tDT:2017-02-23\tPU:11V6WR1' -I 250,50 ~/assessment/dnaseq/data/reference/hg19.fa.gz ~/assessment/dnaseq/data/trimmed_fastq/trimmed_data_1P ~/assessment/dnaseq/data/trimmed_fastq/trimmed_data_2P > ~/assessment/dnaseq/data/aligned_data/NGS0001.sam

#convert sam to bam
cd ~/assessment/dnaseq/data/aligned_data
samtools view -h -b NGS0001.sam > NGS0001.bam 

#sort the bam file
samtools sort NGS0001.bam > NGS0001_sorted.bam

#generate bai file
samtools index NGS0001_sorted.bam

#mark duplicates
picard MarkDuplicates I= NGS0001_sorted.bam O= NGS0001_sorted_marked.bam M=marked_dup_metrics.txt

#index and filter the bam file
samtools index NGS0001_sorted_filtered.bam
samtools view -F 1796  -q 20 -o NGS0001_sorted_filtered.bam NGS0001_sorted_marked.bam

#generate stats 
samtools flagstat NGS0001_sorted_filtered.bam

samtools idxstats NGS0001_sorted_filtered.bam

picard CollectInsertSizeMetrics -H NGS0001_histogram.pdf -I NGS0001_sorted_filtered.bam -O NGS0001_insertsize.txt

bedtools coverage -a NGS0001_sorted_filtered.bam -b ~/assessment/annotation.bed > NGS0001_coverageoutput.txt

#calling variants with freebayes
zcat ~/assessment/dnaseq/data/reference/hg19.fa.gz > ~/assessment/dnaseq/data/reference/hg19.fa  #converts file to fa

samtools faidx ~/assessment/dnaseq/data/reference/hg19.fa

freebayes --bam ~/assessment/dnaseq/data/aligned_data/NGS0001_sorted_filtered.bam --fasta-reference ~/assessment/dnaseq/data/reference/hg19.fa --vcf ~/assessment/dnaseq/results/NGS0001.vcf
bgzip ~/assessment/dnaseq/results/NGS0001.vcf #zips the file

tabix -p vcf ~/assessment/dnaseq/results/NGS0001.vcf.gz #indexes the file

#Apply filters to VCF
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ~/assessment/dnaseq/results/NGS0001.vcf.gz >  /assessment/dnaseq/results/NGS0001_filtered.vcf

bedtools intersect -header -wa -a ~/assessment/dnaseq/results/NGS0001_filtered.vcf -b ~/assessment/dnaseq/data/annotation.hg19.bed > ~/assessment/dnaseq/results/ NGS0001_filtered_a.vcf

bgzip ~/assessment/dnaseq/results/ NGS0001_filtered_a.vcf
tabix -p vcf ~/assessment/dnaseq/results/ NGS0001_filtered_a.vcf.gz

#set up annovar 
tar -zxvf annovar.latest.tar.gz
cd ~/annovar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/

#VCF to annovar
cd ~/annovar
./convert2annovar.pl -format vcf4 ~/assessment/dnaseq/results/ NGS0001_filtered_a.vcf.gz > ~/assessment/dnaseq/results/ NGS0001_filtered_a.avinput
#convert to csv output

./table_annovar.pl ~/assessment/dnaseq/results/ NGS0001_filtered_a.avinput humandb/ -buildver hg19  -out ~/assessment/dnaseq/results/ NGS0001_filtered_a -remove  -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout
