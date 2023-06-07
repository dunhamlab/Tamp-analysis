#!/bin/bash
#$ -S /bin/bash
#$ -wd /net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool/
#$ -o /net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool/OutStd/
#$ -e /net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool/ErrStd/
#$ -N Barseq_plotting
#$ -l mfree=4G

## run with: qsub Tamp_plots.sh MN_D MN_R DMSOvRadicicol

# load modules
module load modules modules-init modules-gs
module load python/2.7.3 numpy/1.13.3 biopython/1.70
module load gcc/8.1.0
module load R/latest

CON1=$1
CON2=$2
EXPNAME=$3
DIR=/net/dunham/vol2/Abby/Projects/Tamp_data/fastq/radicicol/Tamp_pool
SCRIPTDIR=${DIR}/Scripts
MERGEDIR=${DIR}/merge_files

Rscript ${SCRIPTDIR}/190228_fitness_plots.R ${SCRIPTDIR} ${MERGEDIR} ${CON1} ${CON2} ${EXPNAME}

Rscript ${SCRIPTDIR}/190301_plot_chromosomes.R ${SCRIPTDIR} ${MERGEDIR} ${CON1} ${CON2}

