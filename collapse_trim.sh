#!/bin/bash
#PBS -P FFbigdata
#PBS -N collapse_trim
#PBS -l select=1:ncpus=24:mem=48GB
#PBS -l walltime=12:00:00
#PBS -q alloc-op

#load modules
module load bmftools
module load bwa
module load picard
module load java
module load anaconda
module load python/3.5.1
module load R/3.3.2
module load macs
module load fastqc
module load trimgalore

cd /project/RDS-SMS-FFbigdata-RW/Genetics/umi_c9/BaseCalls/

# gunzip all files
# gunzip *fastq.gz

#bmftools - collapse .fastq read based on UMI. Note - fastq name needs to be in format '%_R*.fastq'
for read1 in $(ls *R1.fastq);do

#Inputs - Sample
read2=${read1%_R1.fastq}_R2.fastq
read_umi=${read1%_R1.fastq}_I2.fastq
out_prefix=${read1%_R1.fastq}_collapsed
tmp_prefix=tmp

echo $read1 $read2 $read_umi $out_prefix

# # bmftools collapse cant handle >299bp PE reads. Therefore we trim the reads to 299bp before collapsing
# cut -c -299 $read1 > trunc_${read1}
# cut -c -299 $read2 > trunc_${read2}

# # run
# bmftools collapse secondary -s 5 -p 12 -o $tmp_prefix -f $out_prefix -i $read_umi trunc_${read1} trunc_${read2}

# # fastqc
# fastqc ${read1%_R1.fastq}_collapsed.R1.fq.gz ${read1%_R1.fastq}_collapsed.R2.fq.gz

# trimming
## trim_galore -q 30 --paired --illumina --fastqc ${read1%_R1.fastq}_collapsed.R1.fq.gz ${read1%_R1.fastq}_collapsed.R2.fq.gz
# NOTE - trimming parameters were changed due to 1) bmftoools annotating "#" in quality score (Phred=2) for any base that may be different between collapsed reads i.e. 1 base is N.
# Therefore quality trimming steps were removed an only adapter trimming steps were kept.
trim_galore -q 0 --paired ${read1%_R1.fastq}_collapsed.R1.fq.gz ${read1%_R1.fastq}_collapsed.R2.fq.gz

#rename for STRetch. Note: need to be in format '%_R*.fastq.gz'
rename ${read1%_R1.fastq}_collapsed.R1_val_1.fq.gz ${read1%_R1.fastq}_collapsed_trimmed_R1.fastq.gz ${read1%_R1.fastq}_collapsed.R1_val_1.fq.gz
rename ${read1%_R1.fastq}_collapsed.R2_val_2.fq.gz ${read1%_R1.fastq}_collapsed_trimmed_R2.fastq.gz ${read1%_R1.fastq}_collapsed.R2_val_2.fq.gz

done

