#!/bin/bash
#PBS -P FFbigdata
#PBS -N read_gen_stretch
#PBS -l select=1:ncpus=8:mem=32GB
#PBS -l walltime=20:00:00
#PBS -M zacchatt@gmail.com

Load required modules
module load blast+/2.7.1
module load bedtools
module load artfastgen
module load art
module load python
module load trimmomatic
module load parallel
module load seqtk


# make new output directory
indir=/project/RDS-SMS-FFbigdata-RW/zacc/c9_test/insilico_c9_150321_zc
code_location=/project/RDS-SMS-FFbigdata-RW/zacc/c9_test/
cd $indir

####################
### HOUSEKEEPING ###
####################

REPEAT_VAR='GGCCCC'
REPEATS_PULLDOWN=4
REPEATS_PATHTEST=2000 # hard-coded as 2000
REPEATS_INTERVALTEST=20 # no longer needed
FASTQ_REPS=3
PATH_CHR='chr9'
PATH_START='27573526'
PATH_END='27573544'
EXT_BP=1000
hg19_genome=/project/RDS-FSC-Frontier-RW/local_lib/STRetch/reference-data/hg19.STRdecoys.sorted.fasta
ArtFastGen_location=/project/RDS-SMS-FFbigdata-RW/zacc/c9_test/ArtificialFastqGenerator_19_05_2015
njobs=4

${code_location}/read_generator_1.2.sh \
    $REPEAT_VAR \
    $REPEATS_PULLDOWN \
    $REPEATS_PATHTEST \
    $REPEATS_INTERVALTEST \
    $FASTQ_REPS \
    $PATH_CHR \
    $PATH_START \
    $PATH_END \
    $EXT_BP \
    $hg19_genome \
    $ArtFastGen_location \
    $njobs

# run STRetch
# load modules
module load anaconda
module load python/3.5.1
source activate STR
module load R/3.3.2
STR_location=/project/RDS-FSC-Frontier-RW/local_lib/STRetch

# change to directory containing .fastq files
cd ${indir}/fastq_files

# ## whole - genome ##
gzip *fastq

# create controls - Uncomment the line in ${STR_location}/pipelines/pipeline_config.groovy use a set of WGS samples as controls, or specify your own
for rep_fq in $(ls *_R1.fastq.gz| cut -c 1-3 | sort | uniq);do
    echo ${rep_fq}*fastq.gz
    ${STR_location}/tools/bin/bpipe run ${STR_location}/pipelines/STRetch_wgs_pipeline.groovy ${rep_fq}*fastq.gz
done

mkdir tmp
for i in *locus_counts;do
echo $i
/project/RDS-SMS-FFbigdata-RW/zacc/c9_test/mv_empty.sh $i
done

# recalculate median coverage from .bam file
for input_bam in *bam;do
input_bed=/project/RDS-FSC-Frontier-RW/local_lib/STRetch/reference-data/hg19.simpleRepeat_period1-6_dedup.sorted.bed
/project/RDS-FSC-Frontier-RW/local_lib/STRetch/tools/bin/goleft covmed $input_bam $input_bed | cut -f 1 > ${input_bam%%.bam}.median_cov
done
# check and then divide by 100 to get fraction
for i in *.median_cov;do
a=$(awk '{print $1}' $i)
echo "$a / 100" | bc -l > $i
done

# run STRetch
rm STRs.tsv
${STR_location}/tools/bin/bpipe run ${STR_location}/pipelines/STRetch_wgs_pipeline.groovy *.fastq.gz




