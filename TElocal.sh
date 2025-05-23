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

# TElocal
## Step 1: Generate BAM files needed for TElocal
### Ensure star_1pass_n.sh and star_2pass_unsorted.sh are prepared and modified
cd $SAMPLE_DIR  # Enter the original directory where samples are stored

mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/TElocal
cd $SAMPLE_DIR/$SAMPLE_NAME/TElocal

TElocal -b $SAMPLE_DIR/$SAMPLE_NAME/bam/star/$SAMPLE_NAME/${SAMPLE_NAME}.sorted.2passAligned.out.bam --GTF $SAMPLE_DIR/reference/gencode.v42.annotation.gtf --TE $SAMPLE_DIR/reference/GRCh38_GENCODE_rmsk_TE.gtf.locInd --sortByPos --project ./${SAMPLE_NAME}_TElocal

head -n 1 ${SAMPLE_NAME}_TElocal.cntTable > ${SAMPLE_NAME}_header.txt
sed '1d' ${SAMPLE_NAME}_TElocal.cntTable | grep ":" | cat ${SAMPLE_NAME}_header.txt - > ${SAMPLE_NAME}_TE.txt
sed '1d' ${SAMPLE_NAME}_TElocal.cntTable | grep -v ":" | cat ${SAMPLE_NAME}_header.txt - > ${SAMPLE_NAME}_genes.txt

### Filter non-zero TE
awk '$2 != 0' ${SAMPLE_NAME}_TE.txt > ${SAMPLE_NAME}_filtered_TE.txt

### Filter TE counts greater than 10
awk '$2 > 10' ${SAMPLE_NAME}_TE.txt > ${SAMPLE_NAME}_filtered_10_TE.txt
awk -F: '{print $1}' ${SAMPLE_NAME}_filtered_10_TE.txt > ${SAMPLE_NAME}_TE_name.txt

tail -n +2 ${SAMPLE_NAME}_TE_name.txt | while read id;
do
    grep -w "${id}" $SAMPLE_DIR/documents/GRCh38_GENCODE_rmsk_TE.gtf.locInd.locations >> ${SAMPLE_NAME}_TE_name_location.txt
done

awk 'BEGIN {OFS="\t"} {split($2, a, "[:-]"); print a[1], a[2], a[3], $1, 0, a[4]}' ${SAMPLE_NAME}_TE_name_location.txt > ${SAMPLE_NAME}_TE_name_location.bed

awk 'BEGIN {OFS="\t"} $6 == "" {$6="-"} {print}' ${SAMPLE_NAME}_TE_name_location.bed > ${SAMPLE_NAME}_TE_name_location_new.bed

bedtools getfasta -fi $SAMPLE_DIR/reference/GRCh38.p13.genome.fa -bed ${SAMPLE_NAME}_TE_name_location_new.bed -fo ${SAMPLE_NAME}_TE_name_location_new.fa

awk 'NR==FNR {header[FNR]=$4; next} /^>/ {sub(/>.*/, ">" header[++count])} 1' ${SAMPLE_NAME}_TE_name_location_new.bed ${SAMPLE_NAME}_TE_name_location_new.fa > ${SAMPLE_NAME}_TE.fa

getorf -find 1 -minsize 24 -sequence ${SAMPLE_NAME}_TE.fa -outseq ${SAMPLE_NAME}_TElocal.getorf.fa
