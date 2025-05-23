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

# Fusion detection
# Prepare the sample_list.txt file in advance

# Step 1
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/fusion
cd $SAMPLE_DIR/$SAMPLE_NAME/fusion
perl -c config.txt -fd $SAMPLE_DIR/$SAMPLE_NAME/fusion/SOAPfuse -l $SAMPLE_DIR/$SAMPLE_NAME/fusion/SOAPfuse/${SAMPLE_NAME}_list.txt -o $SAMPLE_DIR/$SAMPLE_NAME/fusion/SOAPfuse/Fusion_detection

# Step 2
input_file="$SAMPLE_NAME/fusion/SOAPfuse/Fusion_detection/final_fusion_genes/$SAMPLE_NAME/analysis/For_peptides_analysis/${SAMPLE_NAME}.trans.fusion.peptide.chain"
output_file="$SAMPLE_NAME/fusion/SOAPfuse/Fusion_detection/final_fusion_genes/$SAMPLE_NAME/analysis/For_peptides_analysis/${SAMPLE_NAME}_fusion.fa"

# Read and process each line of the input file
while IFS= read -r line; do
    col1=$(echo "$line" | awk '{print $1}')
    col2=$(echo "$line" | awk '{print $2}')
    echo ">$col1" >> "$output_file"
    echo "$col2" >> "$output_file"
done < "$input_file"

# Step 3
sed 's/|//g' $SAMPLE_DIR/$SAMPLE_NAME/fusion/SOAPfuse/Fusion_detection/final_fusion_genes/$SAMPLE_NAME/analysis/For_peptides_analysis/${SAMPLE_NAME}_fusion.fa > $SAMPLE_DIR/$SAMPLE_NAME/fusion/SOAPfuse/Fusion_detection/final_fusion_genes/$SAMPLE_NAME/analysis/For_peptides_analysis/${SAMPLE_NAME}_fusion_new.fa

cd $SAMPLE_DIR/$SAMPLE_NAME/fusion/SOAPfuse/Fusion_detection/final_fusion_genes/$SAMPLE_NAME/analysis/For_peptides_analysis
getorf -sequence ${SAMPLE_NAME}_fusion_new.fa -outseq ${SAMPLE_NAME}_fusion.getorf.fa -find 1 -minsize 24
