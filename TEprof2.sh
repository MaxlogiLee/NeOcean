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

#TEProf2
# Step 1: Run annotation on each GTF file
cd $SAMPLE_DIR
mkdir -p $SAMPLE_NAME/TEProf2
cd $SAMPLE_NAME/TEProf2
cp $SAMPLE_NAME/bam/star/$SAMPLE_NAME/${SAMPLE_NAME}.stringtie.gtf ./
cp $SAMPLE_NAME/bam/star/$SAMPLE_NAME/${SAMPLE_NAME}.sorted.2passAligned.out.bam  /${SAMPLE_NAME}.stringtie.bam
cp $SAMPLE_NAME/bam/star/$SAMPLE_NAME/${SAMPLE_NAME}.sorted.2passAligned.out.bam.bai ./${SAMPLE_NAME}.stringtie.bam.bai

python $SAMPLE_DIR/documents/rmskhg38_annotate_gtf_update_test_tpm.py ${SAMPLE_NAME}.stringtie.gtf $SAMPLE_DIR/documents/arguments.txt

# Step 2: Process annotation files to get rough estimate of relative expression of transcript versus all other transcripts of gene
python $SAMPLE_DIR/documents/annotationtpmprocess.py ${SAMPLE_NAME}.stringtie.gtf_annotated_filtered_test_all

# Step 3: Aggregate annotation samples across samples
$SAMPLE_DIR/documents/aggregateProcessedAnnotation.R -a $SAMPLE_DIR/documents/arguments.txt -e T

# Step 4: Calculate Read Information
mkdir filterreadstats
$SAMPLE_DIR/documents/commandsmax_speed.py filter_combined_candidates.tsv ./

# Run the commands (requires bedtools)
parallel -j 10 < filterreadcommands.txt

# Combine all the read information files
find ./filterreadstats/ -name "*.stats" -type f -maxdepth 1 -print0 | xargs -0 -n128 -P1 grep e > ./resultgrep_filterreadstatsdone.txt
cat ./resultgrep_filterreadstatsdone.txt | sed 's/\:/\t/g' > ./filter_read_stats.txt

# Step 5: Filter Candidates based on read information
$SAMPLE_DIR/documents/filterReadCandidates.R
rm Step4.RData  # Optionally retain for data robustness

# Step 6: Merge with Reference GTF
gffread -E candidate_transcripts.gff3 -T -o candidate_transcripts.gtf
echo candidate_transcripts.gtf > cuffmergegtf.list
cuffmerge -o ./merged_asm_full -g $SAMPLE_DIR/reference/gencode.v42.primary_assembly.annotation.gtf cuffmergegtf.list
mv ./merged_asm_full/merged.gtf reference_merged_candidates.gtf
gffread -E reference_merged_candidates.gtf -o- > reference_merged_candidates.gff3

# Step 7: Annotate Merged GTF
python2.7 $SAMPLE_DIR/documents/rmskhg38_annotate_gtf_update_test_tpm_cuff.py reference_merged_candidates.gff3 $SAMPLE_DIR/documents/arguments.txt

# Step 8: Calculate Transcript-Level Expression
find ./ -maxdepth 1 -name "*bam" | while read file ; do
    xbase=${file##*/}
    echo "samtools view -q 255 -h \"$file\" | stringtie - -o \"${xbase%.*}.gtf\" -e -b \"${xbase%.*}_stats\" -p 2 -m 100 -c 1 -G reference_merged_candidates.gtf" >> quantificationCommands2.txt
done
cat quantificationCommands2.txt | parallel -j 20

# Step 9: Process and map expression output
$SAMPLE_DIR/documents/mergeAnnotationProcess.R

# Process stringtie transcript annotation files to get relevant information and aggregate
## Obtaining intron coverage information
find . -maxdepth 2 -name "*i_data.ctab" > ctab_i.txt
cat ctab_i.txt | while read ID ; do 
    fileid=$(echo "$ID" | awk -F "/" '{print $2}'); 
    cat <(printf 'chr\tstrand\tstart\tend\t'${fileid/_stats/}'\n') <(grep -F -f candidate_introns.txt $ID | awk -F'\t' '{ print $2"\t"$3"\t"$4"\t"$5"\t"$6 }') > ${ID}_cand ; 
done

cat <(find . -maxdepth 2 -name "*i_data.ctab_cand" | head -1 | while read file ; do cat $file | awk '{print $1"\t"$2"\t"$3"\t"$4}' ; done;) > table_i_all
find . -maxdepth 2 -name "*i_data.ctab_cand" | while read file ; do paste -d'\t' <(cat table_i_all) <(cat $file | awk '{print $5}') > table_i_all_temp; mv -f table_i_all_temp table_i_all; done ;

## Obtaining the transcript-level expression information for candidates
ls ./*stats/t_data.ctab > ctablist.txt
cat ctablist.txt | while read file ; do echo "stringtieExpressionFrac.py $file" >> stringtieExpressionFracCommands.txt ; done;
cat stringtieExpressionFracCommands.txt | parallel -j 20

### Aggregate the stats
ls ./*stats/t_data.ctab_frac_tot > ctab_frac_tot_files.txt
ls ./*stats/t_data.ctab_tpm > ctab_tpm_files.txt

cat <(echo "TranscriptID") <(find . -name "*ctab_frac_tot" | head -1 | while read file ; do sort $file | awk '{print $1}' ; done;) > table_frac_tot
cat ctab_frac_tot_files.txt | while read file ; do 
    fileid=$(echo "$file" | awk -F "/" '{print $2}') ; 
    paste -d'\t' <(cat table_frac_tot) <(cat <(echo ${fileid/_stats/}) <(sort $file | awk '{print $2}')) > table_frac_tot_temp; 
    mv table_frac_tot_temp table_frac_tot; 
done ;

cat <(echo "TranscriptID") <(find . -name "*ctab_tpm" | head -1 | while read file ; do sort $file | awk '{print $1}' ; done;) > table_tpm
cat ctab_tpm_files.txt | while read file ; do 
    fileid=$(echo "$file" | awk -F "/" '{print $2}') ; 
    paste -d'\t' <(cat table_tpm) <(cat <(echo ${fileid/_stats/}) <(sort $file | awk '{print $2}')) > table_tpm_temp; 
    mv table_tpm_temp table_tpm; 
done ;

cat <(head -1 table_frac_tot) <(grep -Ff candidate_names.txt table_frac_tot) > table_frac_tot_cand
cat <(head -1 table_tpm) <(grep -Ff candidate_names.txt table_tpm) > table_tpm_cand

# Step 10: Quantification processing, sample identification, and final table creation
$SAMPLE_DIR/documents/finalStatisticsOutput_new.R -e T -a $SAMPLE_DIR/documents/arguments.txt
rm Step10.RData

# Step 11: Translating transcripts with Kozak method and generating FASTA of RNA sequences
$SAMPLE_DIR/documents/translationPart1.R -g 'BSgenome.Hsapiens.UCSC.hg38'

# Step 12: Get ORF
getorf -find 1 -minsize 24 -sequence candidates.fa -outseq $SAMPLE_NAME\_TEprof2.getorf.fa
