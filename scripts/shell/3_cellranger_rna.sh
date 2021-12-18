#!/bin/sh

#PBS -q LARGE
#PBS -l select=1:ncpus=64:mem=512gb
#PBS -N cellranger_rna
#PBS -e ${PBS_O_WORKDIR}
#PBS -m e
#PBS -o ${PBS_O_WORKDIR}
#PBS -M johannes.nicolaus@gmail.com

# make directories
dir_out="${PBS_O_WORKDIR}/../../results/data/cellranger_rna"

mkdir -p ${dir_out}
cd ${dir_out}

#cellranger count --id=neg \
#    --localcores=64 \
#    --localvmem=512 \
#    --localmem=512 \
#    --fastqs=${PBS_O_WORKDIR}/../../data/fastq/rna_neg/Dlxneg_0_1_unknow_flowcell \
#    --transcriptome=${PBS_O_WORKDIR}/../../data/reference/refdata-gex-mm10-2020-A

cellranger count --id=pos \
    --localcores=64 \
    --localvmem=512 \
    --localmem=512 \
    --fastqs=${PBS_O_WORKDIR}/../../data/fastq/rna_pos/Dlxpos_0_1_unknow_flowcell \
    --transcriptome=${PBS_O_WORKDIR}/../../data/reference/refdata-gex-mm10-2020-A


# create csv file to aggregate both data
echo -e "sample_id,molecule_h5" > ${dir_out}/aggr.csv
echo -e "negative,${dir_out}/neg/outs/molecule_info.h5" >> ${dir_out}/aggr.csv
echo -e "positive,${dir_out}/pos/outs/molecule_info.h5" >> ${dir_out}/aggr.csv

# perform cellranger aggr
cellranger aggr --id=neg_pos --csv=aggr.csv
