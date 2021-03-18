#!/bin/bash

# define inputs
REPEAT_VAR=${1}
echo "repetitive element = $REPEAT_VAR"

REPEATS_PULLDOWN=${2}
echo "number of repetitive elements = $REPEATS_PULLDOWN"

REPEATS_PATHTEST=${3}
echo "$REPEATS_PATHTEST"

REPEATS_INTERVALTEST=${4}
echo "$REPEATS_INTERVALTEST"

FASTQ_REPS=${5}
echo "$FASTQ_REPS"

PATH_CHR=${6}
echo "$PATH_CHR"

PATH_START=${7}
echo "$PATH_START"

PATH_END=${8}
echo "$PATH_END"

EXT_BP=${9}
echo "$EXT_BP"

hg19_genome=${10}
echo "$hg19_genome"

ArtFastGen_location=${11}
echo "$ArtFastGen_location"

njobs=${12}
echo "$njobs"

######################
### PULLDOWN FASTA ###
######################

# Create FASTA file with repetitive sequence for use in BLAST
echo ">rg_${REPEATS_PULLDOWN}_${REPEAT_VAR}" > rg_${REPEATS_PULLDOWN}_${REPEAT_VAR}.fsa
printf -v v "%-*s" "${REPEATS_PULLDOWN}" ""; echo "${v// /${REPEAT_VAR}}" >> rg_${REPEATS_PULLDOWN}_${REPEAT_VAR}.fsa

# Make human genome database from hg19 genome
makeblastdb -in $hg19_genome -dbtype nucl -title hg19_genome -out hg19_genome.db

# Local hg19 query
ws=$(tail -1 rg_${REPEATS_PULLDOWN}_${REPEAT_VAR}.fsa | wc -c)
word_size=$(expr $ws - 1)
blastn -task megablast -word_size $word_size -query rg_${REPEATS_PULLDOWN}_${REPEAT_VAR}.fsa -evalue 1000 -dust no -db hg19_genome.db -out rg_${REPEATS_PULLDOWN}_${REPEAT_VAR}_hg19local30_out -outfmt "7 sacc sstart send evalue qacc qstart qend"
# remote hg38 query - this isnt returning the exact reuslt we want yet
# blastn -task megablast -word_size 16 -query 5C9Repeats.fsa -evalue 1000 -dust no -db GPIPE/9606/109/ref_top_level -out 5C9Repeats_hg38remote_out -remote -outfmt "7 qacc sacc evalue qstart qend sstart send"

# Set expected value threshold (default is 10)
awk '!/^#/ && ($4 < 10)' rg_${REPEATS_PULLDOWN}_${REPEAT_VAR}_hg19local30_out > temp1.bed

# Reformat bed file to ensure smallest (start) value on the left
awk -F'\t' '$2 > $3 { $0 = $1 FS $3 FS $2 } 1' temp1.bed > temp2.bed

# Sort and Merge overlapping features
sort -k1,1 -k2,2n temp2.bed | awk -v OFS='\t' '{print $1, $2, $3}' > sorted.bed
bedtools merge -i sorted.bed > sorted_merged.bed

# Extend output +/- bp (set to 150bp)
awk -v x=$EXT_BP {'printf ("%s\t%s\t%s\n", $1, $2 - x, $3 + x)'} sorted_merged.bed > temp3.bed

# Get fasta sequences - This represents the C9 pulldown library.
bedtools getfasta -fi $hg19_genome -bed temp3.bed -fo rg_${REPEATS_PULLDOWN}_${REPEAT_VAR}.fa

# Split to make reads from each enriched region
mkdir split_dir
csplit -s -z rg_${REPEATS_PULLDOWN}_${REPEAT_VAR}.fa '/>/' '{*}'
for i in xx* ; do \
  n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ; \
  mv "$i" "split_dir/$n.fa" ; \
done

# cleanup temporary files
rm temp1.bed temp2.bed sorted.bed

########################
### ARTIFICIAL FASTA ###
########################

# Create bed file of expansion location
echo $PATH_CHR $PATH_START $PATH_END > temp4.bed

# Expand bed file to include +/- basepairs
awk -v x=$EXT_BP {'printf ("%s\t%s\t%s\n", $1, $2 - x, $3 + x)'} temp4.bed > temp5.bed

# Get fasta sequences
bedtools getfasta -fi $hg19_genome -bed temp5.bed -fo path_repeat.fa

# make total .bed for downstream analysis
cat temp5.bed temp3.bed > art_pulldown.bed

# select start and end FASTA sequence +/- base-pairs added
string_length=$(sed -n '2p' path_repeat.fa | wc -m)
COUNT=`expr $string_length - $EXT_BP`
start=$(sed -n '2p' path_repeat.fa | cut -c1-${EXT_BP})
end=$(sed -n '2p' path_repeat.fa | cut -c${COUNT}-${string_length})

# Create artificial FASTA files with pathogenic repeats to be tested, place into split directory for cleanliness
rm repeat_lengths_list
for i in $(seq 0 1 7); do
echo "scale=2;2^$i" | bc >> repeat_lengths_list
done
seq 200 100 1000 >> repeat_lengths_list
seq 1200 200 2000 >> repeat_lengths_list

for i in $(cat repeat_lengths_list); do
    echo ">repeat_element" > split_dir/rg_pathrep${i}_${REPEAT_VAR}.fsa
    middle=$(printf -v v "%-*s" "${i}" ""; echo "${v// /${REPEAT_VAR}}")
    echo $start$middle$end >> split_dir/rg_pathrep${i}_${REPEAT_VAR}.fsa
done

########################
### ARTIFICIAL FASTQ ###
########################
# Create artificial fastq files for each FASTA file. Note: variables contained within single inverted commas cause errors in bash and so '"${ used for location variable to expand
#parallel -a fa_list --xapply -j $njobs --eta 'artfastgen -O {1.}_R -R {1} -S ">repeat_element" -RL 300 -TLM 400 -URQS true -SE true -F1 '"${ArtFastGen_location}"'/test1.fastq -F2 '"${ArtFastGen_location}"'/test2.fastq' 
#parallel -a fa_list --xapply -j $njobs --eta 'artfastgen -O {1.}_R -R {1} -S ">repeat_element" -CMP 100 -CSD 0.1 -RL 301 -TLM 400 -CMGCS 1 -CMPGC 1 -URQS true -SE true -F1 '"${URQS_read1}"' -F2 '"${URQS_read2}"''
rm fa_list fsa_list
ls split_dir/*fa >> fa_list
ls split_dir/*fsa >> fsa_list
#parallel -a fa_list --xapply -j $njobs --eta 'artfastgen -O {1.}_R -R {1} -S ">" -CMP 1000000 -RL 301 -TLM 600 -GCR 20 -CMPGC 1'

# run fa and fsa seperately as artificial c9 will generate ~ 10X more reads at same depth due to GC richness
#parallel -a fa_list --xapply -j $njobs --eta 'artfastgen -O {1.}_R -R {1} -S ">" -CMP 1000000 -RL 301 -TLM 600 -GCR 20 -CMPGC 0.8'
#parallel -a fsa_list --xapply -j $njobs --eta 'artfastgen -O {1.}_R -R {1} -S ">" -CMP 100000 -RL 301 -TLM 600 -GCR 20 -CMPGC 0.8'
parallel -a fa_list --xapply -j $njobs --eta 'artfastgen -O {1.}_R -R {1} -S ">" -CMP 100000 -RL 301 -TLM 600'
parallel -a fsa_list --xapply -j $njobs --eta 'artfastgen -O {1.}_R -R {1} -S ">" -CMP 100000 -RL 301 -TLM 600'

# rename and gzip .fastq files
rename _R.1.fastq _R1.fastq split_dir/*_R.1.fastq
rename _R.2.fastq _R2.fastq split_dir/*_R.2.fastq

# sample reads to create replicate files with a portion of .fastq from each .fasta location
mkdir fastq_files
for rep_fq in $(cat repeat_lengths_list);do
    rep_fq2=split_dir/rg_pathrep${rep_fq}_GGCCCC
for FQ in $(ls split_dir/chr*_R1.fastq);do
    echo $FQ
for reads in {0.1,0.2,0.5,1}; do
#for reads in {5000,25000,50000,100000,200000}; do
#reads=1000
seqtk sample -s 100 ${rep_fq2}_R1.fastq $reads > fastq_files/${rep_fq}R${reads}1_R1.fastq
seqtk sample -s 100 ${rep_fq2}_R2.fastq $reads > fastq_files/${rep_fq}R${reads}1_R2.fastq
seqtk sample -s 50 ${rep_fq2}_R1.fastq $reads > fastq_files/${rep_fq}R${reads}2_R1.fastq
seqtk sample -s 50 ${rep_fq2}_R2.fastq $reads > fastq_files/${rep_fq}R${reads}2_R2.fastq
seqtk sample -s 10 ${rep_fq2}_R1.fastq $reads > fastq_files/${rep_fq}R${reads}3_R1.fastq
seqtk sample -s 10 ${rep_fq2}_R2.fastq $reads > fastq_files/${rep_fq}R${reads}3_R2.fastq

   # reads=$(echo "$cell_reads * $reads_input" | bc | awk '{print int($1+0.5)}')
seqtk sample -s 100 $FQ $reads >> fastq_files/${rep_fq}R${reads/./}1_R1.fastq
seqtk sample -s 100 ${FQ%%_R1.fastq}_R2.fastq $reads >> fastq_files/${rep_fq}R${reads/./}1_R2.fastq
seqtk sample -s 50 $FQ $reads >> fastq_files/${rep_fq}R${reads/./}2_R1.fastq
seqtk sample -s 50 ${FQ%%_R1.fastq}_R2.fastq $reads >> fastq_files/${rep_fq}R${reads/./}2_R2.fastq
seqtk sample -s 10 $FQ $reads >> fastq_files/${rep_fq}R${reads/./}3_R1.fastq
seqtk sample -s 10 ${FQ%%_R1.fastq}_R2.fastq $reads >> fastq_files/${rep_fq}R${reads/./}3_R2.fastq
done
done
done



