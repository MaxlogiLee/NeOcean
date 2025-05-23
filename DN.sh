#!/bin/bash
# Define variables for sample name and directory
source ~/miniconda3/etc/profile.d/conda.sh #"~" need to be replaced with the actual path.
SAMPLE_NAME="<sample_name>"
SAMPLE_DIR="<sample_directory>"
mkdir -p $SAMPLE_DIR/reference
mkdir -p $SAMPLE_DIR/documents
mkdir -p $SAMPLE_DIR/STAR_ensemble

# Enter the original directory where the samples are stored
cd $SAMPLE_DIR
mkdir -p $SAMPLE_NAME

# De novo analysis
cd $SAMPLE_DIR
mkdir -p $SAMPLE_NAME/de_novo
cd $SAMPLE_NAME/de_novo

stringtie $SAMPLE_NAME/bam/star/$SAMPLE_NAME/${SAMPLE_NAME}.sorted.2passAligned.out.bam -o ${SAMPLE_NAME}_2pass.gtf -m 100 -c 1 -p 20

echo ${SAMPLE_NAME}_2pass.gtf > ${SAMPLE_NAME}_cuffmergegtf.list

cuffmerge -o ./merged_asm -g $SAMPLE_DIR/reference/gencode.v42.primary_assembly.annotation.gtf ${SAMPLE_NAME}_cuffmergegtf.list

mv ./merged_asm/merged.gtf ${SAMPLE_NAME}_de_novo_merged_candidates.gtf

gffread ${SAMPLE_NAME}_de_novo_merged_candidates.gtf -g $SAMPLE_DIR/reference/GRCh38.primary_assembly.genome.fa -w ${SAMPLE_NAME}_transcript_de_novo.fa

getorf -find 1 -minsize 24 -sequence ${SAMPLE_NAME}_transcript_de_novo.fa -outseq ${SAMPLE_NAME}_de_novo.getorf.fa
