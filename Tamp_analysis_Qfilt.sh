#!/bin/bash
#$ -S /bin/bash
#$ -l mfree=32G
#$ -cwd
#$ -N Barseq
#$ -o /net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool/StdOut/
#$ -e /net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool/StdErr/

# Feb. 9, 2020

# load modules
module load modules modules-init modules-gs
module load fastx-toolkit/0.0.14
module load python/2.7.3

# run w/: python run_qsub.py radicicol_samples.txt Tamp_analysis.sh

PRIMER=$1
# adjust based on filenames for a particular run
READ1=${PRIMER}_S*_R1_001.fastq.gz
READ2=${PRIMER}_S*_R2_001.fastq.gz

DIR=/net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool
SCRIPTDIR=/net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool/Scripts

# set up folder structure
cd ${DIR}
mkdir -p ${DIR}/merge_files
MERGEDIR=${DIR}/merge_files

echo ***filtering reads***
# isolate barcode from forward read
zcat ${DIR}/${READ1} | fastx_trimmer -l 12 -Q33 -z -o ${DIR}/${PRIMER}_BC_R1.fastq.gz
# isolate UMI from reverse read
zcat ${DIR}/${READ2} | fastx_trimmer -l 20 -Q33 -z -o ${DIR}/${PRIMER}_BC_R2.fastq.gz

# merge files with paste, then fix formatting
paste -d "" <(zcat ${DIR}/${PRIMER}_BC_R2.fastq.gz) <(zcat ${DIR}/${PRIMER}_BC_R1.fastq.gz) | awk '{if (NR % 4 == 1) print $1; else if (NR % 4 == 3) print "+"; else print $0;}' > ${DIR}/${PRIMER}_both.fastq

# filter based on quality (Q20)
fastq_quality_filter -i ${DIR}/${PRIMER}_both.fastq -q 20 -p 100 -Q33 -z -o ${DIR}/${PRIMER}_both_Q20_filtered.fastq.gz

gzip ${DIR}/${PRIMER}_both.fastq
rm ${DIR}/${PRIMER}_BC_R1.fastq.gz
rm ${DIR}/${PRIMER}_BC_R2.fastq.gz

# merge 50bp R2 and 12bp R1 with awk command
echo ***making merged barcode files***
zcat ${PRIMER}_both_Q20_filtered.fastq.gz | awk '{if(NR%4==2)print substr($0,1,20)"\t"substr($0,21)}' > ${MERGEDIR}/F${PRIMER}_merge.merge


######################################
# for merging with additional files, run basic qfilter to merge on those separately cat to append to main files
######################################


# count barcodes in each merge.merge file
echo ***counting barcodes***

# if masterlist in scripts directory
# MASTERLIST=${SCRIPTDIR}/180327_gene_master_list.txt

MASTERLIST=/net/dunham/vol2/Abby/Projects/Tamp_data/Tamp_summary_files/180327_gene_master_list.txt
python ${SCRIPTDIR}/190311_count_barcodes.py ${MERGEDIR}/F${PRIMER}_merge.merge ${MASTERLIST} F${PRIMER}_counts.txt





