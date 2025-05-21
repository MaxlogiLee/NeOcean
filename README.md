# NeOcean
NeOcean is a tumor neoantigen analysis platform that integrates proteogenomic strategies to overcome the limitations of traditional approaches relying solely on DNA mutation burden. By incorporating gene fusions, intron retention, chimeric transcripts, and transposable elements, and validating peptides via HLA-pulldown mass spectrometry, NeOcean enables the high-efficiency identification of non-canonical neoantigens, with a detection rate nearly a thousand times higher than conventional methods. It features an optimized RNA-seq analysis pipeline, machine learning–based neoantigen prediction, and high-throughput mass spectrometry validation. Systematic analyses across tumor types reveal distinct sequence features and conserved expression patterns of non-canonical neoantigens. Ongoing functional studies will further evaluate their immunogenic potential, providing a new framework for neoantigen-based cancer immunotherapy.  
# Outline
# Requirements
## 1.Software
Each step is executed in a separate environment.  
All of the software tools listed above can be installed individually.  
This pipeline has only been tested with the specified versions, and newer versions may cause errors or unexpected behavior.  
### (1) Generate BAM file
STAR >= 2.7.11b  
samtools >= 1.21  
stringtie >= 2.2.1  
### (2) Fusion detection
(https://github.com/Nobel-Justin/SOAPfuse)  
SOAPfuse >= 0.01  
PERL >= 5.32.1  
EMBOSS >= 6.5.7  
### (3) Intron retention
(https://github.com/genemine/iread)  
STAR >= 2.7.11b  
samtools >= 1.2  
bedtools >= v2.30.0  
EMBOSS >= 6.5.7  
Bedops >= 2.4.20 ( available at: https://bedops.readthedocs.io/en/latest/)  
PERL >= v5.32.1 (The PERL module Parallel::ForkManager needs to be installed to support multi-core computing)  
python >= 2.7.18 (python module: argparse. If not installed, run 'pip install argparse' from shell to install)  
### (4) TElocal
(https://github.com/mhammell-laboratory/TElocal)  
TElocal >= 1.1.1  
Python >= 2.7.18  
bedtools >= v2.30.0  
EMBOSS >= 6.5.7  
### (5) De novo
stringtie >= 2.2.1  
cufflinks >= 2.2.1  
gffread >= 0.12.7  
EMBOSS >= 6.5.7  
### (6) TEProf2
(https://github.com/twlab/TEProf2Paper)  
TEProF2 >= v0.1  
stringtie >= 1.3.3  
samtools >= 1.3.1  
cufflinks >= 2.2.1  
python 2.7 (cPickle, pytabix 0.1)  
R >= 3.4.1 (ggplot2, bsgenome.hsapiens.ucsc.hg38 (or genome of your choosing), Xmisc, reshape2)  
**It is recommended to follow the official website for installing TEprof2: https://github.com/twlab/TEProf2Paper.**
## 2.Reference Files
### (1) Defaults
Download Link hg38: [External Download Link](https://wangftp.wustl.edu/~nshah/rnapipeline_public_link/rnapipelinerefhg38.tar.gz)  
Download Link hg19: [External Download Link](https://wangftp.wustl.edu/~nshah/rnapipeline_public_link/rnapipelinerefhg19.tar.gz)  
Download Link mm39: [External Download Link](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz)  
Download Link mm10: [External Download Link](https://wangftp.wustl.edu/~nshah/rnapipeline_public_link/rnapipelinerefmm10.tar.gz)  
Download Link: [Homo_sapiens.GRCh38.113.gtf](https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz)  
Download Link: [Mus_musculus.GRCm39.113.gtf](https://ftp.ensembl.org/pub/release-113/gtf/mus_musculus/Mus_musculus.GRCm39.113.gtf.gz)  
### (2) Fusion detection
```
wget https://sourceforge.net/projects/soapfuse/files/SOAPfuse_Package/SOAPfuse-v1.27.tar.gz  
wget ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.chr.gtf.gz  
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz  
``` 
### (3) Intron retention
```
wget https://raw.githubusercontent.com/genemine/iread/refs/heads/master/iread.py  
wget ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz  
wget ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz  
wget https://raw.githubusercontent.com/genemine/iread/master/meta/intron_human_3877.bed  
wget https://raw.githubusercontent.com/genemine/iread/refs/heads/master/meta/intron_mouse_3875.bed
```
### (4) TElocal
```
wget https://raw.githubusercontent.com/mhammell-laboratory/TElocal/refs/heads/master/TElocal  *You can download the software using wget, or refer to the setup.py file on GitHub(https://github.com/mhammell-laboratory/TElocal/blob/master/setup.py) and install it directly using the command:  
python setup.py install  
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.p13.genome.fa.gz  
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz
```
### (5) De novo
```
wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.primary_assembly.annotation.gtf.gz  
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz  
```
### (6) TEProf2
All of the following files are from:https://github.com/twlab/TEProf2Paper  
#### (A) Gencode Dictionary
1.Download Gencode GTF reference file  
`wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.primary_assembly.annotation.gtf.gz`  
2.Sort the File  
`cat <GENCODE GTF>  | awk '{if($3=="transcript"||$3=="exon"||$3=="start_codon"){print}}' |  awk -F "; " '{print $0"\t"$2}' > <OUTPUT_sorted.gtf>`   
3.Use custom script to create dictionary  
`genecode_to_dic.py <OUTPUT_sorted.gtf>`  
This will generate 2 files: (1) genecode_plus.dic and (2) genecode_minus.dic  
4.Rename as needed  
#### (B) Repeatmasker Files
1.Download Repeatmasker Annotation BED. We use the one from the UCSC Table Browser (For this step download it in BED format):[External Download Link](https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=693256623_kh0RR0o6vajdA2WLLTA8OeaAPNB6&clade=mammal&org=Human&db=hg38&hgta_group=rep&hgta_track=hg38Patch11&hgta_table=0&hgta_regionType=genome&position=chr12%3A20816734-20825794&hgta_outputType=bed&hgta_outFileName=repeats.bed)  
2.Sorted the bed file  
`cat | sort -k1,1 -k2,2n > <SORTED.BED>`  
3.bgzip the file  
`bgzip rmsk.bed > rmsk.bed.gz`  
4.Create tabix index  
`tabix -p bed rmsk.bed.gz`  
5.Both the bgzipped file and the tabix index must be in the same directory  
#### (C) Repeatmasker Subfamily, Family, Class Mapping
`cat repeatmasker.lst | awk 'NR>1{print $11"\t"$12"\t"$13}' | sort | uniq > repeatmasker_hg38_description_uniq.lst`  
#### (D) Intron Annotations
1.download genecode_introns.py  
`wget https://raw.githubusercontent.com/twlab/TEProf2Paper/refs/heads/main/bin/genecode_introns.py`    
2.Use custom script to create intron annotations  
`genecode_introns.py <OUTPUT_sorted.gtf>`  
This will generate 2 files: (1) <OUTPUT_sorted.gtf>_introns_plus and (2) <OUTPUT_sorted.gtf>_introns_minus  
3.Sort the intron annotations bed file  
```
sort -k1,1 -k2,2n -k3,3n <OUTPUT_sorted.gtf>_introns_plus > <OUTPUT_sorted.gtf>_introns_plus_sorted  
sort -k1,1 -k2,2n -k3,3n <OUTPUT_sorted.gtf>_introns_minus > <OUTPUT_sorted.gtf>_introns_minus_sorted  
```
4.bgzip the file (bgzip comes with samtools)  
```
bgzip <OUTPUT_sorted.gtf>_introns_plus_sorted > <OUTPUT_sorted.gtf>_introns_plus_sorted.gz  
bgzip <OUTPUT_sorted.gtf>_introns_minus_sorted > <OUTPUT_sorted.gtf>_introns_minus_sorted.gz
```
5.Create tabix index  
```
tabix -p bed <OUTPUT_sorted.gtf>_introns_plus_sorted.gz  
tabix -p bed <OUTPUT_sorted.gtf>_introns_minus_sorted.gz
```
6.Both the bgzipped file and the tabix index need to be in the same folder  
#### (E) Gene Filter List
For large sets of analysis, it might be computationally advantageous to limit the analysis to only a small set of genes for subsequent analysis. Create a file with a genesymbol per line. Example:   
```
TP63  
GAPDH  
SYT1  
PKIB
```
#### (F) Input files at runtime
*aggregateProcessedAnnotation.R*  
`wget https://raw.githubusercontent.com/twlab/TEProf2Paper/main/bin/aggregateProcessedAnnotation.R`  
*annotationtpmprocess.py*  
`wget https://raw.githubusercontent.com/twlab/TEProf2Paper/refs/heads/main/bin/annotationtpmprocess.py`  
*arguments.txt*  
`wget https://raw.githubusercontent.com/twlab/TEProf2Paper/refs/heads/main/bin/arguments.txt`  
**The reference files in arguments.txt need to be replaced with the actual paths to your own files.**  
*commandsmax_speed.py*  
`wget https://raw.githubusercontent.com/twlab/TEProf2Paper/refs/heads/main/bin/commandsmax_speed.py`  
*filterReadCandidates.R*  
`wget https://raw.githubusercontent.com/twlab/TEProf2Paper/refs/heads/main/bin/filterReadCandidates.R`  
*finalStatisticsOutput.R*  
`wget https://raw.githubusercontent.com/twlab/TEProf2Paper/refs/heads/main/bin/finalStatisticsOutput.R`  
*mergeAnnotationProcess.R*  
`wget https://raw.githubusercontent.com/twlab/TEProf2Paper/refs/heads/main/bin/mergeAnnotationProcess.R`  
*rmskhg38_annotate_gtf_update_test_tpm.py*  
`wget https://raw.githubusercontent.com/twlab/TEProf2Paper/refs/heads/main/bin/rmskhg38_annotate_gtf_update_test_tpm.py`  
*rmskhg38_annotate_gtf_update_test_tpm_cuff.py*  
`wget https://raw.githubusercontent.com/twlab/TEProf2Paper/refs/heads/main/bin/rmskhg38_annotate_gtf_update_test_tpm_cuff.py`  
*stringtieExpressionFrac.py*  
`wget https://raw.githubusercontent.com/twlab/TEProf2Paper/refs/heads/main/bin/stringtieExpressionFrac.py`  
*translationPart1.R*  
`wget https://raw.githubusercontent.com/twlab/TEProf2Paper/refs/heads/main/bin/translationPart1.R`  
*rmsk_annotate_bedpe_speed.py*  
`wget https://raw.githubusercontent.com/twlab/TEProf2Paper/refs/heads/main/bin/rmsk_annotate_bedpe_speed.py`  
# Usage
Once all files and environment configurations are ready, we can proceed to run the pipeline.  
## Quick Run
**It is recommended to use the shell script we provided for direct execution.**  
Run the following command   
`nohup bash NeOcean.sh &`  
**Since the script 'NeOcean.sh' activates the environment multiple times, you need to add source '/path/conda.sh' at the very beginning of the script.  You can find the path to 'conda.sh' by running the command:**  
```
conda info | grep 'base environment'  
Example output:base environment : /public/home/miniconda3  (writable)  
Then your "conda.sh" is located at:/public/home/miniconda3/etc/profile.d/conda.sh  
source /public/home/lijing/miniconda3/etc/profile.d/conda.sh  
```
**You can also split the scripts and run them individually by following these steps**  
## (1) Generate BAM file  
`nohup bash generate_bam.sh &` 
Note:The generate_bam.sh script includes two sub-scripts: star_1pass_n.sh and star_2pass_unsorted_rv.sh.  
> 
### Input File(s)
(1)GRCh38.primary_assembly.genome.fa  
(2)Homo_sapiens.GRCh38.113.gtf  
(3)raw data(fasta)  
### Output File(s)  
(1)STAR_ensemble
> Index directory (contains multiple .sjdb*, .genome* files)

(2)*sample*.SJ.out.tab
> The 1-pass script primarily generates files for use in subsequent 2-pass builds of the splice junction index.

(2)*sample*.sorted.2passAligned.out.bam  
(3)*sample*.sorted.2passAligned.out.bam.bai  
(4)*sample*.stringtie.gtf  
## (2) Fusion detection
`nohup bash fusion.sh &`  
### Input File(s)
(1)config.txt  
(2)*sample*_list.txt
(3)raw data(fasta)  
### Output File(s)
(1)*sample*_fusion.getorf.fa  
## (3) Intron retention
`nohup bash IR.sh &`  
### Input File(s)  
(1)raw data(fasta)  
(2)intron_human_3877.bed(intron_mouse_3875.bed)  
(3)Homo_sapiens.GRCh38.77.gtf.gz  
(4)Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz  
### Output File(s)  
(1)*sample*_IR.getorf.fa  
## (4) TElocal
`nohup bash TElocal.sh &`  
### Input File(s)  
(1)*sample*.sorted.2passAligned.out.bam  
(2)gencode.v42.annotation.gtf  
(3)GRCh38_GENCODE_rmsk_TE.gtf.locInd  
(4)GRCh38_GENCODE_rmsk_TE.gtf.locInd.locations  
(5)GRCh38.p13.genome.fa.gz  
### Output File(s)
(1)*sample*_TElocal.getorf.fa  
## (5) De novo
`nohup bash DN.sh &`  
### Input File(s)  
(1)*sample*.sorted.2passAligned.out.bam  
(2)gencode.v42.primary_assembly.annotation.gtf  
(3)GRCh38.primary_assembly.genome.fa  
### Output File(s)
(1)*sample*_2pass.gtf  
> Transcript structure prediction results generated by StringTie (GTF format)  
*-m 100: minimum transcript length is 100bp;  
*-c 1: minimum 1 read per transcript;   
*-p 20: use 20 threads.  

(2)*sample*_cuffmergegtf.list  
> List file containing the GTF files from the previous step, used as input for cuffmerge.  

(3)merged_asm  
merged.gtf  
> The merged de novo transcript prediction GTF file integrates pre-existing annotations with new transcripts.  

(4)*sample*_de_novo_merged_candidates.gtf  
> Rename merged transcripts as candidate files with sample names for easy follow-up.  

(5)*sample*_transcript_de_novo.fa  
> FASTA sequence files of transcripts extracted from GTF structures.  

(6)*sample*_de_novo.getorf.fa  
> Open reading frame (ORF) sequences extracted from each transcript using getorf, defaults to protein sequence, FASTA format.  
*-find 1: Find ORFs on the justice strand;  
*-minsize 24: the shortest ORF length is 24bp (i.e. 8 amino acids);  
*-outseq: output ORF protein sequence file.  

## (6) TEProf2  
`nohup bash DN.sh &`  
### Input File(s)  
(1)*sample*.sorted.2passAligned.out.bam  
(2)*sample*.stringtie.gtf  
(3)arguments.txt  
Note:
> These are all the arguments that are needed. The following arguments.txt file would work:
```
rmsk	rmskhg38.bed6.gz
rmskannotationfile repeatmasker_description_uniq.lst
gencodeplusdic	genecode_plus_hg38.dic
gencodeminusdic	genecode_minus_hg38.dic
```
Note：
> The following Arguments are available(If you don't need these options, delete)

**focusgenes:** The program has two outputs (1) on a focus set of genes (2) with all genes. This file lists the genes that the user wants to filter for originally (Gene Filter List)  
**plusintron:** Tabix file of all the plus strand introns (Intron Annotations)  
**minusintron:** Tabix file of all the minus strand introns (Intron Annotations)
(4)All TEProF2 Required R Script Files
### Output File(s)  
(1)*sample*_TEprof2.getorf.fa  
***Please go to the official documentation for the intermediate files that run all of the above R scripts.(Except for De novo）***
