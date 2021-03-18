# repeSeq_scripts
_repeseq_scripts_ is a collection of scripts used for the 
1) Artificial generation of sequencing reads replicating a capture experiment of a 4 repeat “GGCCCC” probe within the human genome.
2) Alignment and repeat expansion quantification from repeSeq capture

## Instructions
	code_location=/Users/zacc/USyd/UMI_repeatexpansions/manuscript_scripts

## 1) Articifical read generation and alignment
	${code_location}/launch_read_generator_align.pbs
	# NOTE - artificial reads can be generated independently by ${code_location}/bash/read_generator_1.2.sh

## 2a) Collapse and trimming of sWGS & repeSeq 
	${code_location}/collapse_trim.sh

## 2b) Alignment and repeat expansion quantification of repeSeq using STRetch
	${code_location}/repeSeq_STRetch.sh

## 2c) Intersections (overlap) between repeSeq and sWGS
	${code_location}/bed_interesect_analysis.sh

## 2d) script for counting STR within reads
	for i in *collapsed_trimmed_R1.fastq.gz;do
	echo $i
	zcat $i | awk 'NR%4==2' | grep "CCGGGG" | wc -l
	zcat $i | awk 'NR%4==2' | grep "CCGGGGCCGGGG" | wc -l
	done

## 2e) script for plotting sequencing metrics
	${code_location}/general_plotting.R

## 2f) script for analysing enrichment of genomic features using LOLA
	${code_location}/lola_analysis_repeseq.R