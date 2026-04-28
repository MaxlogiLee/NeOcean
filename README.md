# NeOcean: Tumor Neoantigen Discovery Pipeline v2.0

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

NeOcean is a standardized, modular computational pipeline for discovering tumor-specific neoantigens from RNA sequencing and mass spectrometry data. It integrates multiple analytical approaches including de novo transcript assembly, transposable element expression, intron retention, gene fusion detection, and somatic mutation analysis to comprehensively profile potential neoantigen sources.

## Features

- **Modular Architecture**: Each analytical step is encapsulated as an independent module with standardized input/output interfaces
- **Configuration-Driven**: All paths and parameters managed via YAML configuration files — zero hard-coded paths
- **Checkpoint & Resume**: Built-in checkpoint system enables resuming interrupted runs without recomputation
- **Batch Processing**: Supports single-sample or multi-sample execution via sample-specific configuration files
- **Publication-Ready**: Designed for reproducibility with comprehensive logging, conda environment definitions, and clear dependency specifications

## Pipeline Overview

```
┌─────────────────────────────────────────────────────────────────────────────┐
│                         NeOcean Pipeline v2.0                               │
├─────────────────────────────────────────────────────────────────────────────┤
│                                                                             │
│  ┌──────────────┐                                                           │
│  │  RNA-seq QC  │  (External: fastqc, trimmomatic)                          │
│  └──────┬───────┘                                                           │
│         ▼                                                                   │
│  ┌────────────────────────────────────────────┐                             │
│  │  Module 01: Preprocessing                  │  STAR 1-pass → 2-pass      │
│  │  - STAR alignment (1-pass SJ detection)    │  → StringTie assembly      │
│  │  - Sample-specific genome index            │                             │
│  │  - BAM sort/index & transcript assembly    │                             │
│  └────────────────────────────────────────────┘                             │
│         │                                                                   │
│         ▼ (Parallel execution)                                              │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐   │
│  │ Module 02    │  │ Module 03    │  │ Module 04    │  │ Module 05    │   │
│  │ Fusion       │  │ Intron       │  │ TElocal      │  │ De Novo      │   │
│  │ Detection    │  │ Retention    │  │              │  │ Assembly     │   │
│  │ (SOAPfuse)   │  │ (iREAD)      │  │              │  │              │   │
│  └──────────────┘  └──────────────┘  └──────────────┘  └──────────────┘   │
│  ┌──────────────┐  ┌──────────────┐                                         │
│  │ Module 06    │  │ Module 07    │                                         │
│  │ TEProf2      │  │ Mutation     │                                         │
│  │ (12 steps)   │  │ (pVACtools)  │                                         │
│  └──────────────┘  └──────────────┘                                         │
│         │                                                                   │
│         ▼                                                                   │
│  ┌────────────────────────────────────────────┐                             │
│  │  Module 09: Peptide Screening              │                             │
│  │  - MS data processing                      │                             │
│  │  - Peptide grep against 9 databases        │                             │
│  │  - Coordinate mapping (DN/TEprof2)         │                             │
│  │  - Target-decoy FDR calculation            │                             │
│  │  - Cis/Trans classification                │                             │
│  └────────────────────────────────────────────┘                             │
│                                                                             │
└─────────────────────────────────────────────────────────────────────────────┘
```

## Directory Structure

```
neocean/
├── README.md                          # This file
├── config/
│   ├── config.template.yaml           # Global configuration template
│   └── sample.template.yaml           # Per-sample configuration template
├── workflow/
│   ├── main.sh                        # Main controller script
│   ├── run_module.sh                  # Generic module runner (optional)
│   └── common.sh                      # Common function library
├── modules/
│   ├── 01_preprocess.sh               # STAR alignment & StringTie
│   ├── 02_fusion.sh                   # Fusion detection (SOAPfuse)
│   ├── 03_intron_retention.sh         # Intron retention (iREAD)
│   ├── 04_telocal.sh                  # TE expression (TElocal)
│   ├── 05_denovo.sh                   # De novo assembly
│   ├── 06_teprof2.sh                  # TEProf2 (12-step pipeline)
│   ├── 07_mutation.sh                 # Somatic mutation → neoantigen
│   ├── 08_build_bn.sh                 # Build background normal DB
│   └── 09_peptide_screening.sh        # MS integration & peptide screening
├── tools/
│   ├── python/                        # Custom Python scripts
│   └── r/                             # Custom R scripts
├── envs/                              # Conda environment definitions
└── docs/
    ├── INSTALL.md                     # Installation guide
    ├── PARAMETERS.md                  # Parameter documentation
    └── PIPELINE.md                    # Detailed pipeline description
```

## Quick Start

### 1. Installation

See [docs/INSTALL.md](docs/INSTALL.md) for detailed installation instructions.

```bash
# Clone repository
git clone <repository-url> neocean
cd neocean

# Install conda environments
conda env create -f envs/teprof2.yaml
conda env create -f envs/cufflinks.yaml
```

### 2. Configuration

```bash
# Copy templates
cp config/config.template.yaml config/config.yaml
cp config/sample.template.yaml config/samples/YourSample.yaml

# Edit paths for your server environment
nano config/config.yaml
nano config/samples/YourSample.yaml
```

### 3. Run Pipeline

```bash
# Dry-run to check configuration
bash workflow/main.sh \
    --global-config config/config.yaml \
    --sample-config config/samples/YourSample.yaml \
    --dry-run

# Execute full pipeline
bash workflow/main.sh \
    --global-config config/config.yaml \
    --sample-config config/samples/YourSample.yaml

# Resume from interruption
bash workflow/main.sh \
    --global-config config/config.yaml \
    --sample-config config/samples/YourSample.yaml \
    --resume

# Run specific modules only
bash workflow/main.sh \
    --global-config config/config.yaml \
    --sample-config config/samples/YourSample.yaml \
    --modules preprocess,teprof2,peptide_screening
```

### 4. Output Structure

Each sample generates the following workspace structure:

```
${WORK_DIR}/${SAMPLE_NAME}/
├── logs/                      # Module logs
├── status/                    # Checkpoint files (*.done)
├── star/                      # BAM, GTF from preprocessing
├── fusion/                    # Fusion detection results
├── intron_retention/          # iREAD results
├── telocal/                   # TElocal quantification
├── denovo/                    # De novo assembly
├── teprof2/                   # TEProf2 full outputs
├── mutation/                  # Mutation-derived peptides
└── grep/                      # Peptide screening results
    ├── data/                  # Linked ORFs and databases
    └── res/                   # Final screening outputs
```

## Requirements

### Software
| Tool | Version | Module(s) |
|------|---------|-----------|
| STAR | ≥ 2.7.11b | Preprocessing, Intron Retention |
| SAMtools | ≥ 1.21 | Preprocessing, TEProf2 |
| StringTie | ≥ 2.2.1 | Preprocessing, De Novo, TEProf2 |
| BEDtools | ≥ v2.30.0 | Intron Retention, TElocal |
| gffread | ≥ 0.12.7 | De Novo |
| SOAPfuse | ≥ 0.01 | Fusion Detection |
| TElocal | ≥ 1.1.1 | TElocal |
| TEProf2 | ≥ v0.1 | TEProf2 |
| Cufflinks | ≥ 2.2.1 | De Novo |
| EMBOSS (getorf) | ≥ 6.5.7 | Fusion, IR, TElocal, De Novo, TEProf2 |
| Bedops | ≥ 2.4.20 | Intron Retention |
| PERL | ≥ 5.32.1 | Fusion, Intron Retention |
| Python | 2.7 (TEProf2) / ≥ 3.8 (pipeline) | TEProf2 / Controller |
| R | ≥ 3.4.1 | TEProf2, Peptide Screening |
| seqkit | any | Build BN, Decoy merge |

> **Note:** This pipeline has only been tested with the specified versions. Newer versions may cause errors or unexpected behavior.

### Conda Environments (preserved from original)
- `teprof2` — TEProf2 analysis (Python 2.7 + R: ggplot2, BSgenome, Xmisc, reshape2)
- `cufflinks_env` — Cuffmerge, gffread
- `de_nove` — getorf (EMBOSS)
- `telocal-xly` — TElocal

### Reference Data

#### Genome & Annotation (hg38 example)
```bash
# Genome FASTA
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz

# GENCODE annotation
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.primary_assembly.annotation.gtf.gz
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/gencode.v42.annotation.gtf.gz

# Ensembl alternative
wget https://ftp.ensembl.org/pub/release-113/gtf/homo_sapiens/Homo_sapiens.GRCh38.113.gtf.gz
```

#### Fusion Detection (SOAPfuse)
```bash
wget https://sourceforge.net/projects/soapfuse/files/SOAPfuse_Package/SOAPfuse-v1.27.tar.gz
wget ftp://ftp.ensembl.org/pub/release-110/gtf/homo_sapiens/Homo_sapiens.GRCh38.110.chr.gtf.gz
wget ftp://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/cytoBand.txt.gz
```

#### Intron Retention (iREAD)
```bash
wget https://raw.githubusercontent.com/genemine/iread/refs/heads/master/iread.py
wget ftp://ftp.ensembl.org/pub/release-77/gtf/homo_sapiens/Homo_sapiens.GRCh38.77.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-77/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa.gz
wget https://raw.githubusercontent.com/genemine/iread/master/meta/intron_human_3877.bed
wget https://raw.githubusercontent.com/genemine/iread/refs/heads/master/meta/intron_mouse_3875.bed
```

#### TElocal
```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_42/GRCh38.p13.genome.fa.gz
```

#### TEProf2
See [TEProf2 official repository](https://github.com/twlab/TEProf2Paper) for setup.
Key files: rmsk BED (UCSC), gencode dictionary, intron annotations, and 12 R/Python scripts.

#### Additional Required Files
- STAR genome index (build with `STAR --runMode genomeGenerate`)
- TE local index (`GRCh38_GENCODE_rmsk_TE.gtf.locInd` + `.locations`)
- UniProt human proteome FASTA (`uniprotkb_proteome_UP000005640.fasta`)
- (Optional) Sample HLA typing for mutation module

See [docs/INSTALL.md](docs/INSTALL.md) for the complete dependency list and build instructions.

## Module Details

| Module | Input | Output | Key Tools |
|--------|-------|--------|-----------|
| 01_preprocess | Clean RNA-seq FASTQ | Sorted BAM, StringTie GTF | STAR, SAMtools, StringTie |
| 02_fusion | Clean FASTQ | `*_fusion.getorf.fa` | SOAPfuse, getorf |
| 03_intron_retention | Clean FASTQ | `*_IR.getorf.fa` | STAR, iREAD, bedtools, getorf |
| 04_telocal | Sorted BAM | `*_TElocal.getorf.fa` | TElocal, bedtools, getorf |
| 05_denovo | Sorted BAM | `*_de_novo.getorf.fa` | StringTie, Cuffmerge, gffread, getorf |
| 06_teprof2 | BAM + GTF | `*_TEprof2.getorf.fa` | TEProf2 (12-step), parallel |
| 07_mutation | Tumor BAM/VCF + HLA | `mutation.fa` | pVACtools, hla2seq (placeholder) |
| 08_build_bn | UniProt FASTA | `BN_processed.fa` | seqkit (project-level) |
| 09_peptide_screening | All ORFs + MS data | Cis/Trans peptide lists | R scripts, target-decoy FDR |
