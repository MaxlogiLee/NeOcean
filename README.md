# NeOcean
NeOcean is a tumor neoantigen analysis platform that integrates proteogenomic strategies to overcome the limitations of traditional approaches relying solely on DNA mutation burden. By incorporating gene fusions, intron retention, chimeric transcripts, and transposable elements, and validating peptides via HLA-pulldown mass spectrometry, NeOcean enables the high-efficiency identification of non-canonical neoantigens, with a detection rate nearly a thousand times higher than conventional methods. It features an optimized RNA-seq analysis pipeline, machine learning–based neoantigen prediction, and high-throughput mass spectrometry validation. Systematic analyses across tumor types reveal distinct sequence features and conserved expression patterns of non-canonical neoantigens. Ongoing functional studies will further evaluate their immunogenic potential, providing a new framework for neoantigen-based cancer immunotherapy.
# Outline
# Requirements
## 1.Software
Each step is executed in a separate environment.  
All of the software tools listed above can be installed individually. This pipeline has only been tested with the specified versions, and newer versions may cause errors or unexpected behavior.
### 1.1 Generate BAM file
STAR >= 2.7.11b
samtools >= 1.21
stringtie >= 2.2.1
### 1.2 Fusion detection
(https://github.com/Nobel-Justin/SOAPfuse)  
SOAPfuse >= 0.01  
PERL >= 5.32.1  
EMBOSS >= 6.5.7  
### 1.3 Intron retention
(https://github.com/genemine/iread)  
STAR >= 2.7.11b  
samtools >= 1.2  
bedtools >= v2.30.0  
EMBOSS >= 6.5.7  
Bedops >= 2.4.20 ( available at: https://bedops.readthedocs.io/en/latest/)  
PERL >= v5.32.1 (The PERL module Parallel::ForkManager needs to be installed to support multi-core computing)  
python >= 2.7.18 (python module: argparse. If not installed, run 'pip install argparse' from shell to install)  
### 1.4 TElocal
(https://github.com/mhammell-laboratory/TElocal)  
TElocal >= 1.1.1  
Python >= 2.7.18  
bedtools >= v2.30.0  
EMBOSS >= 6.5.7  
### 1.5 De novo
stringtie >= 2.2.1  
cufflinks >= 2.2.1  
gffread >= 0.12.7  
EMBOSS >= 6.5.7  
### 1.6 TEProf2
(https://github.com/twlab/TEProf2Paper)  
TEProF2 >= v0.1  
stringtie >= 1.3.3  
samtools >= 1.3.1  
cufflinks >= 2.2.1  
python 2.7 (cPickle, pytabix 0.1)  
R >= 3.4.1 (ggplot2, bsgenome.hsapiens.ucsc.hg38 (or genome of your choosing), Xmisc, reshape2)  
**It is recommended to follow the official website for installing TEprof2: https://github.com/twlab/TEProf2Paper.**
## Reference Files
### Defaults
Download Link hg38: [External Download Link](https://wangftp.wustl.edu/~nshah/rnapipeline_public_link/rnapipelinerefhg38.tar.gz)  
Download Link hg19: [External Download Link](https://wangftp.wustl.edu/~nshah/rnapipeline_public_link/rnapipelinerefhg19.tar.gz)  
Download Link mm39: [External Download Link](https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz)  
Download Link mm10: [External Download Link](https://wangftp.wustl.edu/~nshah/rnapipeline_public_link/rnapipelinerefmm10.tar.gz)  
### Fusion detection
`wget https://sourceforge.net/projects/soapfuse/files/SOAPfuse_Package/SOAPfuse-v1.27.tar.gz`
`wget ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.chr.gtf.gz`  
`wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz`  
### Intron retention
`wget https://raw.githubusercontent.com/genemine/iread/refs/heads/master/iread.py `  
`wget ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz`  
`wget ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz`  
`wget https://raw.githubusercontent.com/genemine/iread/master/meta/intron_human_3877.bed`  
`wget https://raw.githubusercontent.com/genemine/iread/master/meta/intron_human_3875.bed`  
### TElocal
`wget https://raw.githubusercontent.com/mhammell-laboratory/TElocal/refs/heads/master/TElocal`  *You can download the software using wget, or refer to the setup.py file on GitHub(https://github.com/mhammell-laboratory/TElocal/blob/master/setup.py) and install it directly using the command:  `python setup.py install`  
`wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.p13.genome.fa.gz`  
`wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz`
### De novo
`wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.primary_assembly.annotation.gtf.gz`
### TEProf2
`wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.primary_assembly.annotation.gtf.gz`  
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
