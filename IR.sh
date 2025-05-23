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

# Intron retention
RAW_DATA_forward="/path/to/raw_data_forward"
RAW_DATA_reverse="/path/to/raw_data_reverse"
cd $SAMPLE_DIR
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/intron_retention
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/intron_retention/iread
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/intron_retention/mapping

cd $SAMPLE_DIR/$SAMPLE_NAME/intron_retention/iread

STAR --twopassMode Basic --quantMode GeneCounts --runThreadN 20 --genomeDir $SAMPLE_DIR/STAR_ensemble --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $SAMPLE_DIR/$SAMPLE_NAME/intron_retention/mapping/${SAMPLE_NAME} --readFilesCommand gunzip -c --readFilesIn "${RAW_DATA_forward}" "${RAW_DATA_reverse}"

cd $SAMPLE_DIR/$SAMPLE_NAME/intron_retention/mapping
grep "Uniquely mapped reads number" ${SAMPLE_NAME}Log.final.out | awk '{print $NF/2}' > result.txt

cd $SAMPLE_DIR/$SAMPLE_NAME/intron_retention/iread
mkdir -p result
mkdir -p translation

# Read the number from result.txt and pass it to -t
python $SAMPLE_DIR/documents/iread.py $SAMPLE_DIR/$SAMPLE_NAME/intron_retention/mapping/${SAMPLE_NAME}Aligned.sortedByCoord.out.bam $SAMPLE_DIR/documents/intron_human_3877.bed -o $SAMPLE_DIR/$SAMPLE_NAME/intron_retention/iread/results -t $(cat result.txt)

cd $SAMPLE_DIR/$SAMPLE_NAME/intron_retention/iread/translation
awk '$7 == "yes" {gsub("-", "\t", $1); print $1}' $SAMPLE_DIR/$SAMPLE_NAME/intron_retention/iread/result/${SAMPLE_NAME}Aligned.sortedByCoord.out.ir.txt > $SAMPLE_DIR/$SAMPLE_NAME/intron_retention/iread/result/${SAMPLE_NAME}_yes.bed

bedtools getfasta -fi $SAMPLE_DIR/reference/GRCh38.primary_assembly.genome.fa -bed $SAMPLE_DIR/$SAMPLE_NAME/intron_retention/iread/result/${SAMPLE_NAME}_yes.bed -fo ${SAMPLE_NAME}_intron_retention.fa -name

getorf -sequence ${SAMPLE_NAME}_intron_retention.fa -outseq ${SAMPLE_NAME}_IR.getorf.fa -find 1 -minsize 24
