# Define variables for sample name and directory
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
samtools sort -@ 40 -o $SAMPLE_NAME/bam/${SAMPLE_NAME}.sorted.2passAligned.out.bam $SAMPLE_NAME/bam/${SAMPLE_NAME}.unsorted.2passAligned.out.bam 

# Index BAM files
samtools index $SAMPLE_NAME/bam/${SAMPLE_NAME}.sorted.2passAligned.out.bam 

# Generate GTF with StringTie
stringtie $SAMPLE_NAME/bam/${SAMPLE_NAME}.sorted.2passAligned.out.bam -o $SAMPLE_NAME/bam/${SAMPLE_NAME}.stringtie.gtf --rf -m 100 -c 1 -p 20 -l ${SAMPLE_NAME} -v 1>$SAMPLE_NAME/bam/${SAMPLE_NAME}.stringtie.lg 2>&1

{
# Fusion detection
# Prepare the sample_list.txt file in advance

# Step 1
cd $SAMPLE_DIR
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
} &

{
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
} &

{
# TElocal
## Step 1: Generate BAM files needed for TElocal
### Ensure star_1pass_n.sh and star_2pass_unsorted.sh are prepared and modified
cd $SAMPLE_DIR  # Enter the original directory where samples are stored

mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/TElocal
cd $SAMPLE_DIR/$SAMPLE_NAME/TElocal

TElocal -b $SAMPLE_DIR/$SAMPLE_NAME/bam/${SAMPLE_NAME}.sorted.2passAligned.out.bam --GTF $SAMPLE_DIR/reference/gencode.v42.annotation.gt f --TE $SAMPLE_DIR/reference/GRCh38_GENCODE_rmsk_TE.gtf.locInd --sortByPos --project ./${SAMPLE_NAME}_TElocal

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
} &

{
# De novo analysis
cd $SAMPLE_DIR
mkdir -p $SAMPLE_NAME/de_novo
cd $SAMPLE_NAME/de_novo

stringtie $SAMPLE_NAME/bam/${SAMPLE_NAME}.sorted.2passAligned.out.bam -o ${SAMPLE_NAME}_2pass.gtf -m 100 -c 1 -p 20

echo ${SAMPLE_NAME}_2pass.gtf > ${SAMPLE_NAME}_cuffmergegtf.list

cuffmerge -o ./merged_asm -g $SAMPLE_DIR/reference/gencode.v42.primary_assembly.annotation.gtf ${SAMPLE_NAME}_cuffmergegtf.list

mv ./merged_asm/merged.gtf ${SAMPLE_NAME}_de_novo_merged_candidates.gtf

gffread ${SAMPLE_NAME}_de_novo_merged_candidates.gtf -g $SAMPLE_DIR/reference/GRCh38.primary_assembly.genome.fa -w ${SAMPLE_NAME}_transcript_de_novo.fa

getorf -find 1 -minsize 24 -sequence ${SAMPLE_NAME}_transcript_de_novo.fa -outseq ${SAMPLE_NAME}_de_novo.getorf.fa
} &

{
#TEProf2
# Step 1: Run annotation on each GTF file
cd $SAMPLE_DIR
mkdir -p $SAMPLE_NAME/TEProf2
cd $SAMPLE_NAME/TEProf2
cp $SAMPLE_NAME/bam/${SAMPLE_NAME}.stringtie.gtf ./
cp $SAMPLE_NAME/bam/${SAMPLE_NAME}.sorted.2passAligned.out.bam  /${SAMPLE_NAME}.stringtie.bam
cp $SAMPLE_NAME/bam/${SAMPLE_NAME}.sorted.2passAligned.out.bam.bai ./${SAMPLE_NAME}.stringtie.bam.bai

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
} &

##Peptide screening
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/code
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/data
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/res
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/res/$SAMPLE_NAME
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/res/peptide
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/data/MSdata
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/data/BN
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/data/NCBI
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/data/Uniport
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/data/1.mutation
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/data/2.Fusion
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/data/3.IR
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/data/4.TE_local
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/data/5.TEprof2
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/data/6.Denovo
mkdir -p $SAMPLE_DIR/$SAMPLE_NAME/grep/data/Decoy
cp $SAMPLE_DIR/$SAMPLE_NAME/fusion/SOAPfuse/Fusion_detection/final_fusion_genes/$SAMPLE_NAME/analysis/For_peptides_analysis/${SAMPLE_NAME}_fusion.getorf.fa $SAMPLE_DIR/$SAMPLE_NAME/grep/data/2.Fusion
cp $SAMPLE_DIR/$SAMPLE_NAME/intron_retention/iread/translation/${SAMPLE_NAME}_IR.getorf.fa $SAMPLE_DIR/$SAMPLE_NAME/grep/data/3.IR
cp $SAMPLE_DIR/$SAMPLE_NAME/TElocal/${SAMPLE_NAME}_TElocal.getorf.fa $SAMPLE_DIR/$SAMPLE_NAME/grep/data/4.TE_local
cp $SAMPLE_DIR/$SAMPLE_NAME/de_novo/${SAMPLE_NAME}_de_novo.getorf.fa $SAMPLE_DIR/$SAMPLE_NAME/grep/data/6.Denovo
cp $SAMPLE_DIR/$SAMPLE_NAME/TEProf2/$SAMPLE_NAME\_TEprof2.getorf.fa $SAMPLE_DIR/$SAMPLE_NAME/grep/data/5.TEprof2
#Step1 Run it in R. The mass spectrometry results in this code include variants, and these variants are also included in the peptide of grep later
cd $SAMPLE_DIR/$SAMPLE_NAME/grep
nohup Rscript code/MS_results.R & #The table in the MSdata folder is named $SAMPLE_NAME.xlsx

#Step2 Peptide extraction
nohup Rscript code/peptide_grep_code.R $SAMPLE_NAME res/peptide/$SAMPLE_NAME\_peptide_list.txt BN data/BN/BN_processed.fa FALSE res/$SAMPLE_NAME &

nohup Rscript code/peptide_grep_code.R $SAMPLE_NAME res/peptide/$SAMPLE_NAME\_peptide_list NCBI data/NCBI/ncbi_protein.fa FALSE res/$SAMPLE_NAME &

nohup Rscript code/peptide_grep_code.R $SAMPLE_NAME res/peptide/$SAMPLE_NAME\_peptide_list Uniport data/Uniport/uniprot_protein.fa FALSE res/$SAMPLE_NAME &

nohup Rscript code/peptide_grep_code.R $SAMPLE_NAME res/peptide/$SAMPLE_NAME\_peptide_list Mutation data/1.mutation/mutation.fa FALSE res/$SAMPLE_NAME &

nohup Rscript code/peptide_grep_code.R $SAMPLE_NAME res/peptide/$SAMPLE_NAME\_peptide_list Fusion data/2.Fusion/${SAMPLE_NAME}_fusion.getorf.fa FALSE res/$SAMPLE_NAME &

nohup Rscript code/peptide_grep_code.R $SAMPLE_NAME res/peptide/$SAMPLE_NAME\_peptide_list IR data/3.IR/${SAMPLE_NAME}_IR.getorf.fa FALSE res/$SAMPLE_NAME &

nohup Rscript code/peptide_grep_code.R $SAMPLE_NAME res/peptide/$SAMPLE_NAME\_peptide_list TElocal data/4.TE_local/${SAMPLE_NAME}_TElocal.getorf.fa FALSE res/$SAMPLE_NAME &
  
nohup Rscript code/peptide_grep_code.R $SAMPLE_NAME res/peptide/$SAMPLE_NAME\_peptide_list TEprof2 data/5.TEprof2/$SAMPLE_NAME\_TEprof2.getorf.fa FALSE res/$SAMPLE_NAME &
  
nohup Rscript code/peptide_grep_code.R $SAMPLE_NAME res/peptide/$SAMPLE_NAME\_peptide_list DN data/6.Denovo/${SAMPLE_NAME}_de_novo.getorf.fa FALSE res/$SAMPLE_NAME &

#Step3 DN2loc
nohup Rscript code/gtf2refBed.R data/6.Denovo/$SAMPLE_NAME\_de_novo_merged_candidates.gtf data/6.Denovo/$SAMPLE_NAME\_de_novo_merged_candidates.refbed &
nohup Rscript code/DN_or_TEprof2_to_chromosome_Location.R res/$SAMPLE_NAME/DN_grep.txt data/6.Denovo/$SAMPLE_NAME\_de_novo_merged_candidates.refbed res/$SAMPLE_NAME/DN2loc_grep.txt FALSE &

###TEprof2 refbed
#Extract exon first
grep "exon" data/5.TEprof2/gtf/$SAMPLE_NAME.stringtie.gtf > data/5.TEprof2/gtf/$SAMPLE_NAME.stringtie_exon.gtf &

nohup Rscript code/gtf2refBed.R data/5.TEprof2/gtf/$SAMPLE_NAME.stringtie_exon.gtf data/5.TEprof2/refbed/$SAMPLE_NAME.stringtie.refbed &

nohup Rscript code/DN_or_TEprof2_to_chromosome_Location.R res/$SAMPLE_NAME/TEprof2_grep.txt data/5.TEprof2/refbed/$SAMPLE_NAME.stringtie.refbed res/$SAMPLE_NAME/TEprof22loc_grep.txt FALSE &

#Step4 Integrate mass spectrometry data with grep results
nohup Rscript code/Integrate_result.R res/$SAMPLE_NAME.xlsx res/sp_file.txt res/database_file.txt &

##Use the script reverse_fasta.py to generate a decoy FASTA file. Modify the input_file and output_file variables in the script, setting the input file to the UniProt FASTA
python code/reverse_fasta.py

##After obtaining the decoy FASTA, check whether the sequences are multiline. If they are, merge them into a single line using seqkit seq -w 0 input > output
seqkit seq -w 0 /data/Decoy/Decoy_uniprotkb.fasta > data/Decoy/Decoy_merge_uniprotkb.fasta

#Execute this command first to generate the file res/target_decoy_signif_sequences.txt
nohup Rscript target_decoy_Calculate_FDR.R res/Uniport_cut_peptidelist.txt res 0.05 data/Uniport/uniprot_protein.fa data/Decoy/Decoy_merge_uniprotkb.fasta &

#Step5 For mass spectrometry results that were not found by grep, continue to cut the peptides and perform UniProt grep again
nohup Rscript code/peptide_grep_code.R target_decoy res/target_decoy_signif_sequences.txt Uniport_cut data/database/uniprot_protein.fa FLASE res &

#Step6 Integrate the UniProt grep results of the cut peptides, and finally generate five files: cis*, trans*, the mass spectrometry results for peptides not matched by grep, and detailed information tables for both cis and trans
nohup Rscript code/Integrate_Cis_trans.R res/Uniport_cut_grep.txt res/Uniport_cut_peptide.txt res/target_decoy_result.xls res/0_unmatched_MS_result.xlsx res/sp_file.txt&
