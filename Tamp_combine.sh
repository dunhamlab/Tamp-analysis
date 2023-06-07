#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool/
#$ -o /net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool/OutStd/
#$ -e /net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool/ErrStd/
#$ -N Barseq_format
#$ -l mfree=12G

# load modules
module load modules modules-init modules-gs
module load python/2.7.3

# group is flask ID plus up- or downtags (MN1D_up)
# run w/ python run_qsub.py radicicol_groups.txt Tamp_combine.sh
GROUP=$1
SAMPFILE=${GROUP}.txt
DIR=/net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool
SCRIPTDIR=${DIR}/Scripts
MERGEDIR=${DIR}/merge_files

# combine timepoint files into single file per flask (up or down)
python ${SCRIPTDIR}/combine_counts.py ${SAMPFILE} ${DIR}
