# NeOcean Parameter Documentation

## Global Configuration (`config/config.yaml`)

### `neocean` — Installation Paths

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `neocean.root_dir` | string | Yes | Absolute path to the NeOcean git repository |
| `neocean.resource_dir` | string | Yes | Absolute path to external resources (references, databases, custom scripts) |

### `reference` — Reference Genomes

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `reference.genome_fasta` | string | Yes | GRCh38 primary assembly FASTA |
| `reference.star_index` | string | Yes | Pre-built STAR genome index directory |
| `reference.gtf` | string | Yes | GENCODE primary assembly annotation GTF |
| `reference.gtf_telocal` | string | Yes | GENCODE basic annotation GTF (for TElocal) |
| `reference.te_index` | string | Yes | TElocal index file (`.locInd`) |
| `reference.te_locations` | string | Yes | TElocal locations file (`.locInd.locations`) |
| `reference.iread_bed` | string | Yes | iREAD intron BED file |

### `resources` — Tool-Specific Resources

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `resources.teprof2.bin_dir` | string | Yes | Directory containing TEProf2 Python/R scripts |
| `resources.teprof2.doc_dir` | string | Yes | Directory containing TEProf2 R documents |
| `resources.teprof2.arguments` | string | Yes | TEProf2 arguments configuration file |
| `resources.soapfuse.program_dir` | string | Yes | SOAPfuse installation directory |
| `resources.soapfuse.config` | string | Yes | SOAPfuse configuration file |
| `resources.iread.program_dir` | string | Yes | iREAD installation directory |
| `resources.grep.code_dir` | string | Yes | Peptide screening R/Python scripts directory |
| `resources.grep.bn_fasta` | string | No | Background normal protein FASTA |
| `resources.grep.uniprot_fasta` | string | Yes | UniProt human proteome FASTA |
| `resources.grep.ncbi_fasta` | string | No | NCBI RefSeq protein FASTA |

### `conda_envs` — Environment Names

| Parameter | Default | Description |
|-----------|---------|-------------|
| `conda_envs.teprof2` | `teprof2` | TEProf2 analysis environment |
| `conda_envs.cufflinks` | `cufflinks_env` | Cuffmerge/gffread environment |
| `conda_envs.denovo` | `de_nove` | EMBOSS getorf environment |
| `conda_envs.telocal` | `telocal-xly` | TElocal environment |

### `compute` — Computational Resources

| Parameter | Default | Description |
|-----------|---------|-------------|
| `compute.threads` | 40 | Maximum threads for parallel operations |
| `compute.parallel_jobs` | 5 | Number of concurrent module jobs |
| `compute.memory_gb` | 128 | Approximate memory limit (GB) |
| `compute.star_threads` | 20 | Threads for STAR alignment |
| `compute.stringtie_threads` | 20 | Threads for StringTie assembly |

### `modules` — Module Switches

Set to `false` to globally disable a module for all samples.

| Parameter | Default | Description |
|-----------|---------|-------------|
| `modules.preprocess` | true | RNA-seq alignment and assembly |
| `modules.fusion` | true | Fusion gene detection |
| `modules.intron_retention` | true | Intron retention analysis |
| `modules.telocal` | true | Transposable element expression |
| `modules.denovo` | true | De novo transcript assembly |
| `modules.teprof2` | true | TEProf2 profiling |
| `modules.mutation` | true | Somatic mutation analysis |
| `modules.build_bn` | false | Build background normal DB (project-level) |
| `modules.peptide_screening` | true | Peptide screening and MS integration |

---

## Sample Configuration (`config/samples/*.yaml`)

### `sample` — Sample Information

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sample.name` | string | Yes | Unique sample identifier |
| `sample.work_dir` | string | Yes | Sample-specific output directory |
| `sample.sequencing_type` | string | Yes | `RNA-seq` / `WGS` / `WES` |

### `sample.rna` — RNA-seq Data

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sample.rna.fq1` | string | Yes | Path to R1 FASTQ (gzipped) |
| `sample.rna.fq2` | string | Yes | Path to R2 FASTQ (gzipped) |
| `sample.rna.read_length` | int | Yes | Read length (bp) |
| `sample.rna.strandness` | string | Yes | Library strandness: `RF` / `FR` / `unstranded` |

### `sample.dna` — DNA Sequencing Data (Optional)

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sample.dna.tumor_bam` | string | No | Tumor WGS/WES BAM |
| `sample.dna.normal_bam` | string | No | Normal WGS/WES BAM (paired) |
| `sample.dna.vcf` | string | No | External somatic VCF (bypasses mutation calling) |

### `sample.hla` — HLA Typing (Optional)

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sample.hla.alleles` | list | No | List of HLA alleles, e.g. `["HLA-A*02:01"]` |
| `sample.hla.hla_file` | string | No | Path to HLA typing result file |

### `sample.pon` — Panel of Normal (Optional)

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `sample.pon.pon_vcf` | string | No | Panel of Normal VCF for tumor-only calling |

### `module_overrides` — Per-Sample Overrides

Override global module switches for individual samples.

```yaml
module_overrides:
  mutation: false    # Skip mutation for this sample
  teprof2: true      # Explicitly enable TEProf2
```

---

## Command-Line Arguments (`main.sh`)

| Argument | Short | Required | Description |
|----------|-------|----------|-------------|
| `--global-config` | `-g` | Yes | Path to global config YAML |
| `--sample-config` | `-s` | Yes | Path to sample config YAML |
| `--modules` | `-m` | No | Comma-separated module list. Default: all |
| `--resume` | `-r` | No | Skip modules with existing checkpoints |
| `--dry-run` | `-n` | No | Print execution plan without running |
| `--help` | `-h` | No | Show help message |

### Module Names for `--modules`

- `preprocess`
- `fusion`
- `intron_retention`
- `telocal`
- `denovo`
- `teprof2`
- `mutation`
- `peptide_screening`

Example:
```bash
bash workflow/main.sh -g config/config.yaml -s config/samples/MCJ.yaml -m preprocess,teprof2
```
