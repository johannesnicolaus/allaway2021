#!/bin/sh

#PBS -q LARGE
#PBS -l select=1:ncpus=64:mem=256gb
#PBS -N bam_to_fastq
#PBS -e ${PBS_O_WORKDIR}
#PBS -m e
#PBS -o ${PBS_O_WORKDIR}
#PBS -M johannes.nicolaus@gmail.com

# download files
mkdir ${PBS_O_WORKDIR}/../../data/fastq
cd ${PBS_O_WORKDIR}/../../data/bam

cellranger bamtofastq --nthreads=64 --reads-per-fastq=10000000000 E13_MGE_Dlx6neg_possorted_genome_bam.bam ../fastq/rna_neg

cellranger bamtofastq --nthreads=64 --reads-per-fastq=10000000000 E13_MGE_Dlx6pos_possorted_genome_bam.bam ../fastq/rna_pos
