#!/bin/bash

mkdir "cal_ov_res_R"

for file in $(ls "cal_ov_R")
do
    echo cal_ov_R/$file
    sed -i -e '$a\' cal_ov_R/$file
    calSOV -f 2 cal_ov_R/$file > cal_ov_res_R/$file
done





	# this removes the file ending from the fastq file to be used as a folder name
	#y=${file%.bar}
	#subdir=${y##*/}


	#mkdir ${STARfiles}/mm_ncbi_index/${subdir}
	# working on ncbi gtf
	#$Stardir --runThreadN 64 --runMode genomeGenerate --genomeDir mm39_ncbi_index/ --sjdbGTFfile mm39.ncbiRefSeq.gtf --sjdbOverhang 100 --outfileNamePrefix "${STARfiles}/mm_ncbi_index/${subdir}" > out_mm_ncbi_index.txt
