#!/bin/bash
id=$1
dir_proj=$SAMPLE_DIR/$SAMPLE_NAME/bam
cd ${dir_proj}
input_fq1=${dir_proj}/${id}_1.fq.gz
input_fq2=${dir_proj}/${id}_2.fq.gz
dir_star=${dir_proj}/star/${id}
mkdir -p ${dir_star}
log_star=${dir_star}/${id}.star.lg
exit_star=${dir_star}/${id}.star.exit
bam_star=${dir_star}/${id}.Aligned.sortedByCoord.out.bam
if [ ! -f ${exit_star} ]
then
echo 1 > ${exit_star}
fi
if [ ! -f ${bam_star} ] || [ `cat ${exit_star}` == 1 ]
then
echo "start STAR for ${id}" `date`
ulimit -v 98765432100
# add --outSAMstrandField intronMotif and --outSAMattributes XS for polyA mRNA unstranded data
# add --outSAMattrIHstart 0 for compatibility with StringTie
STAR --genomeDir $SAMPLE_DIR/STAR_ensemble --readFilesType Fastx --readFilesIn ${input_fq1} ${input_fq2} --readFilesCommand zcat --outFileNamePrefix ${dir_star}/${id}. --outSAMtype None --outSAMmode None --runThreadN 20 --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMattributes NH HI NM MD AS nM XS --outSAMattrIHstart 0 1>${log_star} 2>&1
if [ $? -eq 0 ]
then
echo "SUCCESS STAR for ${id}" `date`
echo 0 > ${exit_star}
else
echo "ERROR STAR for ${id}" `date`
echo 1 > ${exit_star}
fi
else
echo "${bam_star} exist, skip star" `date`
fi
