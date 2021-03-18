#!/bin/bash
#PBS -P FFbigdata
#PBS -N STRetch_hm
#PBS -l select=1:ncpus=24:mem=98GB
#PBS -l walltime=2:00:00
#PBS -M zacchatt@gmail.com

#load modules
module load anaconda
module load python/3.5.1
module load bmftools
module load bwa
module load picard
module load java
module load macs
module load fastqc
module load trimgalore
source activate STR
module load R/3.3.2
STR_location=/project/RDS-FSC-Frontier-RW/local_lib/STRetch

cd /project/RDS-SMS-FFbigdata-RW/Genetics/umi_c9/BaseCalls/STRetch_analysis_retrim
BED=/project/RDS-SMS-FFbigdata-RW/zacc/c9_test/get_bed/art_pulldown.bed #bed file of artificial pulldown

# # STRetch - whole genome
${STR_location}/tools/bin/bpipe run ${STR_location}/pipelines/STRetch_wgs_pipeline.groovy *fastq.gz

# # STRetch - exome
# ${STR_location}/tools/bin/bpipe run -p EXOME_TARGET="/project/RDS-SMS-FFbigdata-RW/zacc/c9_test/get_bed/art_pulldown.bed" \
# 	${STR_location}/pipelines/STRetch_exome_pipeline.groovy *fastq.gz

#mark and sort reads and get family stats on aligned file
for out_prefix in *bam;do

bmftools mark ${out_prefix} > ${out_prefix%%.bam}_marked.bam
bmftools sort --output-fmt BAM ${out_prefix%%.bam}_marked.bam > ${out_prefix%%.bam}_sm.bam

#Get MID usage metrics, family size, and target depth:
bmftools famstats frac ${out_prefix%%.bam}_sm.bam
bmftools depth –sb $BED ${out_prefix%%.bam}_sm.bam

#call peaks from aligned bam file
BAM=${out_prefix%%.bam}_sm.bam
name=${out_prefix%%.bam}_macs_out

macs2 callpeak -t $BAM -g 2.7e9 -B --call-summit -n $name 2> ${name}_macs2.log
macs2 callpeak -t $BAM -g 2.7e9 -B  --broad --broad-cutoff 0.1 -n $name 2> ${name}_macs2.log2

done

