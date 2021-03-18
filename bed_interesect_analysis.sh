#!/bin/bash
#PBS -P FFbigdata
#PBS -N read_gen_stretch
#PBS -l select=1:ncpus=8:mem=32GB
#PBS -l walltime=20:00:00
#PBS -M zacchatt@gmail.com

# load modules
module load python/3.6.5
module load samtools
module load bedtools

# setwd
cd /project/RDS-SMS-FFbigdata-RW/Genetics/umi_c9/BaseCalls/STRetch_analysis_retrim


# get intersect with artificial regions
for BAM in *.bam; do
echo $BAM
bedtools intersect -wb -a tmp -b $BAM > ${BAM}_intersect_artpulll.sam
done

for SAM in *artpulll.sam;do
echo $SAM 
grep "2757" $SAM | grep "chr9" | wc -l
done

# create unique BED files from aligned bam files for genomic region enrichment with LOLA
for BAM in *.bam; do
echo $BAM
bedtools bamtobed -i $BAM > ${BAM%%.bam}.bed
done

# get intersect between repeSeq and sWGS 
bedtools intersect -bed -a UMB1-2.bam -b UMB2-2.bam UMB3-2.bam UMB4-2.bam -u > tmp
grep -v "STR" tmp > inter_all_UMB1-2.bed
bedtools intersect -bed -a UMB1-2.bam -b UMB1-1.bam -u > tmp
grep -v "STR" tmp > inter_repewgs_UMB1-2.bed

bedtools intersect -bed -a UMB2-2.bam -b UMB1-2.bam UMB3-2.bam UMB4-2.bam -u > tmp
grep -v "STR" tmp > inter_all_UMB2-2.bed
bedtools intersect -bed -a UMB2-2.bam -b UMB2-1.bam -u > tmp
grep -v "STR" tmp > inter_repewgs_UMB2-2.bed

bedtools intersect -bed -a UMB3-2.bam -b UMB2-2.bam UMB1-2.bam UMB4-2.bam -u > tmp
grep -v "STR" tmp > inter_all_UMB3-2.bed
bedtools intersect -bed -a UMB3-2.bam -b UMB3-1.bam -u > tmp
grep -v "STR" tmp > inter_repewgs_UMB3-2.bed

bedtools intersect -bed -a UMB4-2.bam -b UMB2-2.bam UMB3-2.bam UMB1-2.bam -u > tmp
grep -v "STR" tmp > inter_all_UMB4-2.bed
bedtools intersect -bed -a UMB4-2.bam -b UMB4-1.bam -u > tmp
grep -v "STR" tmp > inter_repewgs_UMB4-2.bed

# get intersect between sWGS 
bedtools intersect -bed -a UMB1-1.bam -b UMB2-1.bam UMB3-1.bam UMB4-1.bam -u > tmp
grep -v "STR" tmp > tmp2
bedtools merge -i tmp2 > inter_all_UMB1-1_merged.bed

bedtools intersect -bed -a UMB2-1.bam -b UMB1-1.bam UMB3-1.bam UMB4-1.bam -u > tmp
grep -v "STR" tmp > tmp2
bedtools merge -i tmp2 > inter_all_UMB2-1_merged.bed

bedtools intersect -bed -a UMB3-1.bam -b UMB2-1.bam UMB1-1.bam UMB4-1.bam -u > tmp
grep -v "STR" tmp > tmp2
bedtools merge -i tmp2 > inter_all_UMB3-1_merged.bed

bedtools intersect -bed -a UMB4-1.bam -b UMB2-1.bam UMB3-1.bam UMB1-1.bam -u > tmp
grep -v "STR" tmp > tmp2
bedtools merge -i tmp2 > inter_all_UMB4-1_merged.bed

#qsub -I -P FFbigdata -q alloc-op -l select=1:ncpus=4:mpiprocs=4:mem=96gb,walltime=2:00:00

# merge all regions
for i in *inter_*.bed;do
bedtools merge -i $i > test
mv test ${i%%.bed}_merged.bed
done

for i in UMB*.bed;do
echo $i
bedtools merge -i $i > tmp
grep -v "STR" tmp > ${i%%.bed}_merged.bed
done

# count the bases covered 
for i in *_merged.bed;do
echo $i
cat $i | awk -F'\t' 'BEGIN{SUM=0}{ SUM+=$3-$2 }END{print SUM}'
done

# combine scwgs to create universe for LOLA analysis
cat UMB1-1.bed UMB2-1.bed UMB3-1.bed UMB4-1.bed > tmp
grep -v "STR" tmp > scwgs.bed
bedtools merge -i scwgs.bed > scwgs_merged.bed

cat UMB1-2.bed UMB2-2.bed UMB3-2.bed UMB4-2.bed > tmp
grep -v "STR" tmp > repe.bed
bedtools merge -i repe.bed > repe_merged.bed


