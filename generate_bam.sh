#!/bin/bash
# Define variables for sample name and directory
source ~/miniconda3/etc/profile.d/conda.sh #"~" need to be replaced with the actual path.
SAMPLE_NAME="<sample_name>"
SAMPLE_DIR="<sample_directory>"
mkdir -p $SAMPLE_DIR/reference
mkdir -p $SAMPLE_DIR/documents
mkdir -p $SAMPLE_DIR/STAR_ensemble
STAR --runMode genomeGenerate --runThreadN 20 --genomeDir $SAMPLE_DIR/STAR_ensemble --genomeFastaFiles $SAMPLE_DIR/reference/Homo_sapiens.GRch38.dna_sm.primary_assembly.fa --sjdbGTFfile $SAMPLE_DIR/reference/Homo_sapiens.GRch38.77.gtf

# Enter the original directory where the samples are stored
cd $SAMPLE_DIR
mkdir -p $SAMPLE_NAME

# Create bam and gtf files for 1pass and 2pass (used for TElocal, TEprof2, and DN)
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/bam
# 1pass
mkdir -p $SAMPLE_NAME/bam/1pass
bash $SAMPLE_NAME/bam/1pass/star_1pass_n.sh $SAMPLE_NAME

# Generate annotation based on 1pass results
mkdir -p $SAMPLE_NAME/bam/genomedir
STAR --runMode genomeGenerate --genomeDir $SAMPLE_DIR/$SAMPLE_NAME/bam/genomedir/ --genomeFastaFiles $SAMPLE_DIR/reference/GRCh38.primary_assembly.genome.fa --sjdbOverhang 100 --runThreadN 40 --limitSjdbInsertNsj 7000000 --sjdbFileChrStartEnd $SAMPLE_DIR/$SAMPLE_NAME/bam/star/*/*SJ.out.tab

# 2pass
mkdir -p $SAMPLE_NAME/bam/2pass
bash $SAMPLE_NAME/bam/2pass/star_2pass_unsorted_rv.sh $SAMPLE_NAME
 
# Sort BAM files
samtools sort -@ 40 -o $SAMPLE_NAME/bam/star/$SAMPLE_NAME/${SAMPLE_NAME}.sorted.2passAligned.out.bam $SAMPLE_NAME/bam/star/$SAMPLE_NAME/${SAMPLE_NAME}.unsorted.2passAligned.out.bam 

# Index BAM files
samtools index $SAMPLE_NAME/bam/star/$SAMPLE_NAME/${SAMPLE_NAME}.sorted.2passAligned.out.bam 

# Generate GTF with StringTie
stringtie $SAMPLE_NAME/bam/star/$SAMPLE_NAME/${SAMPLE_NAME}.sorted.2passAligned.out.bam -o $SAMPLE_NAME/bam/star/$SAMPLE_NAME/${SAMPLE_NAME}.stringtie.gtf --rf -m 100 -c 1 -p 20 -l ${SAMPLE_NAME} -v 1>$SAMPLE_NAME/bam/star/$SAMPLE_NAME/${SAMPLE_NAME}.stringtie.lg 2>&1
