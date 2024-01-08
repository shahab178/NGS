#!/bin/bash

#Unzip the files
gunzip *.gz
#Concatanate all R1
cat *R1*.fastq > seq_R1.fastq
#Concatanate all R2
cat *R2*.fastq > seq_R2.fastq

chmod 755 /.../FastQC/fastqc

mkdir fastqc_before
mkdir fastqc_after

/.../FastQC/fastqc -o ./fastqc_before -t 4 --nogroup *.fastq

chmod 755 /.../Trimmomatic-x.xx/trimmomatic-x.xx.jar

java -jar /.../Trimmomatic-x.xx/trimmomatic-x.xx.jar PE -threads 2 -phred33 ./seq_R1.fastq ./seq_R2.fastq seq_PE_1.fq seq_SR_1.fq seq_PE_2.fq seq_SR_2.fq LEADING:20 TRAILING:20 SLIDINGWINDOW:4:25 MINLEN:36

/.../FastQC/fastqc -o ./fastqc_after -t 4 --nogroup *.fq
#Here, we have to check if the quality is accaptable. Ifnot, it would be necessary to change the options on the Trimmomatics command

fastx_quality_stats -Q 33 -i seq_PE_1.fq -o seq_PE_1.txt

fastx_quality_stats -Q 33 -i seq_PE_2.fq -o seq_PE_2.txt
#make a txt file of the bases quality 
cat *.txt | awk '{sum += $2} END {print sum}' > seq.quality_bases.txt
#put the reference at th esame folder or should be provided the full path on the command bellow
bowtie2-build reference.fasta reference

bowtie2 -p 2 --no-unal -x reference -1 seq_PE_1.fq -2 seq_PE_2.fq -S seq.sam
#Convert sam to bam (sam can be deleted after)
samtools view -bS seq.sam > seq.bam

samtools sort seq.bam seq.sorted

samtools index seq.sorted.bam
#Having a report of the sequences
samtools depth seq.sorted.bam | awk '{sum+=$3} END {print "Average= ", sum/NR}' > seq.average.txt
#To make the assembled sequence. It can be done by opening with Ugene to save the concensus and coverage.
/.../sam2fasta.py reference.fasta seq.sam consensus.fasta
#De novo assembley
/.../SPAdes-x.xx.x-Linux/bin/spades.py --pe1-1 seq_PE_1.fq --pe1-2 seq_PE_2.fq -o spades_output

cd spades_output
#The result will be save at spade_out folder.
grep '^>' contigs.fasta -c > number_of_contigs.txt








