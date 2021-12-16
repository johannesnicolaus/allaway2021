#!/bin/sh

#PBS -q SMALL
#PBS -l select=1:ncpus=1:mem=64gb
#PBS -N download_data
#PBS -e ${PBS_O_WORKDIR}
#PBS -m e
#PBS -o ${PBS_O_WORKDIR}
#PBS -M johannes.nicolaus@gmail.com

# download files

mkdir ${PBS_O_WORKDIR}/../../data/bam
cd ${PBS_O_WORKDIR}/../../data/bam


# download raw bam files
wget -N https://sra-download.ncbi.nlm.nih.gov/traces/sra78/SRZ/013402/SRR13402802/E13_MGE_Dlx6pos_possorted_genome_bam.bam
wget -N https://sra-download.ncbi.nlm.nih.gov/traces/sra77/SRZ/013402/SRR13402803/E13_MGE_Dlx6neg_possorted_genome_bam.bam
wget -N https://sra-download.ncbi.nlm.nih.gov/traces/sra44/SRZ/013477/SRR13477673/E13_MGE_dlxpos_2_ATAC_possorted_bam.bam
wget -N https://sra-download.ncbi.nlm.nih.gov/traces/sra73/SRZ/013477/SRR13477674/E13_MGE_dlxneg_2_ATAC_possorted_bam.bam

# download index
mkdir ${PBS_O_WORKDIR}/../../data/reference
wget -N https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz -P ${PBS_O_WORKDIR}/../../data/reference
tar -zxvf ${PBS_O_WORKDIR}/../../data/reference/refdata-gex-mm10-2020-A.tar.gz -C ${PBS_O_WORKDIR}/../../data/reference

