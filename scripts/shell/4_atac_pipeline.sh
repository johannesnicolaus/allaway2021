#!/bin/sh

#PBS -q LARGE
#PBS -l select=1:ncpus=64:mem=512gb
#PBS -N findpeaks
#PBS -e ${PBS_O_WORKDIR}
#PBS -m e
#PBS -o ${PBS_O_WORKDIR}
#PBS -M johannes.nicolaus@gmail.com

# make directories
dir_out="${PBS_O_WORKDIR}/../../results/data/atac"
dir_bam="${PBS_O_WORKDIR}/../../data/bam"
dir_ref="${PBS_O_WORKDIR}/../../data/reference"

mkdir -p ${dir_out}
cd ${dir_out}

dir_homer="${HOME}/apps/homer/bin"


# make homer directory
#${dir_homer}/makeTagDirectory ${dir_out}/neg \
#    ${dir_bam}/E13_MGE_dlxneg_2_ATAC_possorted_bam.bam -sspe \
#    -single \
#    -tbp 1

#${dir_homer}/makeTagDirectory ${dir_out}/pos \
#    ${dir_bam}/E13_MGE_dlxpos_2_ATAC_possorted_bam.bam -sspe \
#    -single \
#    -tbp 1

# find peaks using homer
#for tagdir in neg pos
#do
#    ${dir_homer}/findPeaks ${dir_out}/${tagdir} \
#        -style super -o auto \
#        -typical ${dir_out}/${tagdir}/typicalEnhancers.txt \
#        -minDist 5000 \
#        -L 0 -fdr 0.0001
#done

# convert to homer pos file to bed file
for tag_dir in neg pos
do
    mkdir ${dir_out}/bed
    ${dir_homer}/pos2bed.pl ${dir_out}/${tag_dir}/superEnhancers.txt > ${dir_out}/bed/${tag_dir}_superEnhancers.bed
    ${dir_homer}/pos2bed.pl ${dir_out}/${tag_dir}/typicalEnhancers.txt > ${dir_out}/bed/${tag_dir}_typicalEnhancers.bed
done

# merge bed file and remove comment lines
for type in superEnhancers typicalEnhancers
do
cat ${dir_out}/bed/pos_${type}.bed ${dir_out}/bed/neg_${type}.bed | \
    sed '/^#/d' | sort -k1,1 -k2,2n | mergeBed > ${dir_out}/bed/merged_${type}.bed
done


# Annotate stitched peaks with intensity (perform for both tag directories, stimulated and unstimulated cells)
mkdir ${dir_out}/annotated

gunzip ${dir_ref}/gencode.vM25.annotation.gtf.gz

for i in merged_superEnhancers merged_typicalEnhancers
do
    ${dir_homer}/annotatePeaks.pl ${dir_out}/bed/${i}.bed \
        ${dir_ref}/GRCm38.primary_assembly.genome.fa.gz \
        -gtf ${dir_ref}/gencode.vM25.annotation.gtf \
        -norm 1000000 \
        -size given \
        -d ${dir_out}/neg ${dir_out}/pos > ${dir_out}/annotated/${i}_annotated_intensity.txt
done


