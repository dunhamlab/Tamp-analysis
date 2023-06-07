#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool/
#$ -o /net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool/OutStd/
#$ -e /net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool/ErrStd/
#$ -N Barseq_freq
#$ -l mfree=4G

# load modules
module load modules modules-init modules-gs
module load python/2.7.3 numpy/1.13.3 biopython/1.70
module load gcc/8.1.0
module load R/latest

# group is flask ID plus up- or downtags (MN1D_up)
# run w/ python run_qsub.py radicicol_exps.txt Tamp_frequencies.sh

GROUP=$1
DIR=/net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool
SCRIPTDIR=${DIR}/Scripts
MERGEDIR=${DIR}/merge_files

python ${SCRIPTDIR}/190225_calc_frequencies.py ${GROUP} ${DIR}

python ${SCRIPTDIR}/190225_calc_log2.py ${GROUP} ${DIR}

# make directory for slope plots
mkdir -p ${MERGEDIR}/${GROUP}_plots
Rscript ${SCRIPTDIR}/190227_getSlopes.R ${SCRIPTDIR} ${MERGEDIR} ${GROUP}

python ${SCRIPTDIR}/190227_remove_NAs.py ${GROUP} ${MERGEDIR}

Rscript ${SCRIPTDIR}/190227_dataAnalysis.R ${SCRIPTDIR} ${MERGEDIR} ${GROUP}

python ${SCRIPTDIR}/190228_format_files.py ${GROUP} ${DIR}

