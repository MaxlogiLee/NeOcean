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
log_star_2pass=${dir_star}/${id}.star.lg_2pass
exit_star_2pass=${dir_star}/${id}.star.exit_2pass
bam_star_2pass=${dir_star}/${id}.Aligned.Unsorted.out.bam_2pass
if [[ ! -f ${exit_star_2pass} ]]
then
echo 1 > ${exit_star_2pass}
fi
if [[ $(cat ${exit_star}) = 0 ]] && [[ $(cat ${exit_star_2pass}) = 1 ]]
then
echo "STAR for ${id} SUCCESS, start STAR_2pass for ${id}" $(date)
ulimit -v 98765432100
# add --outSAMstrandField intronMotif and --outSAMattributes XS for polyA mRNA unstranded data
# add --outSAMattrIHstart 0 for compatibility with StringTie
STAR --runThreadN 20 --genomeDir $SAMPLE_DIR/$SAMPLE_NAME/bam/genomedir/ --readFilesType Fastx --readFilesIn ${input_fq1} ${input_fq2} --readFilesCommand zcat --outFileNamePrefix ${dir_star}/${id}.unsorted.2pass --outFilterMultimapScoreRange 1 --outFilterMultimapNmax 20 --outFilterMismatchNmax 10 --alignIntronMax 500000 --alignMatesGapMax 1000000 --sjdbScore 2 --alignSJDBoverhangMin 1 --genomeLoad NoSharedMemory --limitBAMsortRAM 98765432100 --outFilterMatchNminOverLread 0.33 --outFilterScoreMinOverLread 0.33 --sjdbOverhang 100 --outSAMunmapped Within --outSAMtype BAM Unsorted --outSAMattributes NH HI NM MD AS nM XS --outSAMattrIHstart 0 --limitSjdbInsertNsj 7000000 --outFilterType BySJout 1>${log_star_2pass} 2>&1
if [[ $? -eq 0 ]]
then
echo "SUCCESS STAR_2pass for ${id}" $(date)
echo 0 > ${exit_star_2pass}
else
echo "ERROR STAR_2pass for ${id}" $(date)
echo 1 > ${exit_star_2pass}
fi
else
if [[ $(cat ${exit_star}) = 1 ]]
then 
echo "STAR for ${id} not done"
fi
if [[ $(cat ${exit_star_2pass}) = 0 ]]
then
echo "STAR_2pass for ${id} EXIST"
fi
fi
