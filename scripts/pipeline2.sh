#!/bin/bash

'''
trimmomatic PE  \
  -threads 4 \
  -phred33 \
  $1 $2 \   #these are the two fast.gz files
  -baseout ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data \
  ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic-0.39-2/adapters/NexteraPE-PE.fa:2:30:10 \
  TRAILING:25 MINLEN:50

fastqc -t 4 /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_1P \
        /home/ubuntu/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_2P

mkdir ~/ngs_course/dnaseq_pipeline/results/fastqc_trimmed_reads #make a new folder to put the fastqc trimmed reads output

mv ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/*fastqc* ~/ngs_course/dnaseq_pipeline/results/fastqc_trimmed_reads/ #move the trimmed fastqc outputs to the results fastqc trimmed reads folder


#make the BWA index
mkdir ~/ngs_course/dnaseq_pipeline/data/reference
mv ~/ngs_course/dnaseq/data/reference/hg19.fa.gz ~/ngs_course/dnaseq_pipeline/data/reference/
bwa index ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz 
'''

#run BWA mem
mkdir ~/ngs_course/dnaseq_pipeline/data/aligned_data

bwa mem -t 4 \
-v 1 \
-R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.WES01\tSM:WES01\tPL:ILLUMINA\tLB:nextera-wes01-blood\tDT:2017-02-23\tPU:HWI-D00119' \
-I 250,50 ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz \
~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_1P ~/ngs_course/dnaseq_pipeline/data/trimmed_fastq/trimmed_data_2P > ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m.sam


#convert sam to bam
cd ~/ngs_course/dnaseq_pipeline/data/aligned_data
samtools view -h -b WES01_chr22m.sam > WES01_chr22m.bam

#sort the bam file
samtools sort WES01_chr22m.bam > WES01_chr22m_sorted.bam

#generate bai file
samtools index WES01_chr22m_sorted.bam

#mark duplicates
picard MarkDuplicates I=WES01_chr22m_sorted.bam O=WES01_chr22m_sorted_marked.bam M=marked_dup_metrics.txt

#index the bam file
samtools index WES01_chr22m_sorted_marked.bam

#filter the BAM file
samtools view -F 1796  -q 20 -o WES01_chr22m_sorted_filtered.bam WES01_chr22m_sorted_marked.bam

samtools index WES01_chr22m_sorted_filtered.bam

#calling variants with freebayes
zcat ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa.gz > ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa  #converts file to fa
samtools faidx ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa

freebayes --bam ~/ngs_course/dnaseq_pipeline/data/aligned_data/WES01_chr22m_sorted_filtered.bam --fasta-reference ~/ngs_course/dnaseq_pipeline/data/reference/hg19.fa --vcf ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf

bgzip ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf #zips the file

tabix -p vcf ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf.gz #indexes the file

#Apply filters to VCF
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m.vcf.gz > ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered.vcf


bedtools intersect -header -wa -a ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered.vcf -b /home/ubuntu/ngs_course/dnaseq_pipeline/data/chr22.genes.hg19.bed > ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf

bgzip ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf

tabix -p vcf ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf.gz

'''
#set up annovar 
tar -zxvf annovar.latest.tar.gz

cd ~/annovar
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar ensGene humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20180603 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
./annotate_variation.pl -buildver hg19 -downdb -webfrom annovar dbnsfp31a_interpro humandb/
'''

#VCF to annovar
cd ~/annovar
./convert2annovar.pl -format vcf4 ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.vcf.gz > ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.avinput

#convert to csv output
./table_annovar.pl ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22.avinput humandb/ -buildver hg19  \
   -out ~/ngs_course/dnaseq_pipeline/results/WES01_chr22m_filtered_chr22 -remove   \
      -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring . -csvout
