<p align="center">
  <strong>A comprehensive sequencing analysis pipeline for small RNA and other short-read data</strong>
</p>

<p align="center">
  <strong>Version 4.0</strong>
</p>

<p align="center">
  <a href="#features">Features</a> •
  <a href="#supported-organisms">Organisms</a> •
  <a href="#installation">Installation</a> •
  <a href="#quick-start">Quick Start</a> •
  <a href="#documentation">Documentation</a> •
  <a href="#output">Output</a> •
  <a href="#contributing">Contributing</a>
</p>

---

## Overview

AnnotationPipeline is a bioinformatics tool designed for processing and analyzing deep-sequencing libraries. It provides a comprehensive workflow from raw reads to publication-ready outputs including annotation counts, genome browser tracks, size profiles, and differential expression analysis.

Originally developed at the [Brennecke Lab](https://www.imba.oeaw.ac.at/research/julius-brennecke/) at IMBA (Institute of Molecular Biotechnology), this pipeline has been optimized for small RNA analysis in *Drosophila melanogaster* but supports multiple organisms and sequencing types.

## Features

### Multi-Platform Sequencing Support
Default values for Normalization, Mismatches, Length-Range and Trimming
| Type | Description | Normalization | Mismatches | Length Range | Trimming |
|------|-------------|----------------------|------------|--------------|----------|
| **sRNAseq** | Small RNA sequencing | 1M miRNAs | 0 | 18-35 | - |
| **sRNAseqIP** | Small RNA immunoprecipitation | 10M reads in fasta | 0 | 18-35 | - |
| **RNAseq** | Standard RNA sequencing | 10M unique-mapping reads | 0 | 18-200 | 6-45 |
| **RNAseq_QuantSeq** | QuantSeq libraries | 10M unique-mapping reads | 0 | 18-200 | - |
| **CLIPseq** | Cross-linking immunoprecipitation | 10M unique-mapping reads | 3 | 18-200 | - |
| **CHIPseq** | Chromatin immunoprecipitation | 10M reads in fasta | 0 | 18-200 | - |
| **GROseq** | Global run-on sequencing | 10M unique-mapping reads | 0 | 18-200 | 6-45 |
| **RIPseq** | RNA immunoprecipitation | - | 0 | 18-200 | - |
| **DNAseq** | DNA sequencing | 10M reads in fasta | 0 | 18-200 | - |
| **CapSeq** | Cap analysis gene expression | 10M unique-mapping reads | 0 | 35-200 | - |

###  Supported Organisms
Drosophila and ASM are fully supported, whereas Cel has limited functionality.
| Organism | Genome Versions | Annotation Source | Default Version |
|----------|-----------------|-------------------|-----------------|
| *Drosophila melanogaster* | dm3, dm6 | FlyBase | r6.40 (dm6), r5.57 (dm3) |
| *Caenorhabditis elegans* | Cel | WormBase | WS292 |
| Custom assemblies | ASM | User-provided | - |

### Analysis Pipeline

1. **Pre-processing**
   - BAM file merging (use `;` delimiter for multiple files)
   - Adaptor trimming (customizable)
   - Random nucleotide (N) trimming
   - Quality filtering
   - SRA download and demultiplexing support

2. **Read Processing**
   - N-containing read filtering
   - Length filtering (configurable min/max)
   - Read trimming (position-based)
   - PolyA trimming support
   - Sequence collapsing for efficient processing

3. **Mapping**
   - rRNA and mitochondrial DNA filtering
   - Genome mapping with configurable mismatch tolerance (0-3)
   - Splice-aware mapping for RNAseq (using STAR)
   - Unique and multi-mapper handling
   - Random multi-mapper assignment option

4. **Annotation**
   - Intersection with comprehensive genome annotations
   - Hierarchical annotation resolution using customizable rulesets
   - Majority-vote collapsing for multi-mappers

5. **Quantification**
   - Gene-level quantification (Salmon)
   - Transposon element quantification (sense/antisense)
   - Spike-in normalization support
   - SLAM-seq support for metabolic labeling

### Output

- **Annotation counts** - Per-category read counts with detailed breakdowns
- **UCSC Browser tracks** - Ready-to-use genome browser visualization with auto-scaling
- **Size profiles** - Small RNA length distributions
- **Differential expression** - Built-in DESeq2 analysis
- **Ping-pong analysis** - piRNA phasing analysis for small RNAs
- **Quality metrics** - Sequencing depth and PCR duplication estimates
- **GEO-ready exports** - Formatted for public data submission
- **BAM exports** - Collapsed or uncollapsed alignment files

## Requirements

### System Requirements

- **Operating System**: Linux (tested on CentOS)
- **Cluster**: SLURM workload manager
- **Storage**: Sufficient disk space for indices and results
- **Memory**: Varies by genome size (typically 32GB+ recommended)

### Dependencies

The pipeline uses **Singularity containers** to manage software dependencies automatically:

- `APmaster.simg` - Core analysis tools (mapping, trimming, etc.)
- `R.simg` - R environment for statistical analysis and plotting

**Included tools** (via containers):
- STAR (RNA-seq alignment)
- Bowtie/Bowtie2 (short read alignment)
- Salmon (transcript quantification)
- Samtools, BEDtools
- R with DESeq2, ggplot2, tidyverse

## Installation

### 1. Clone the Repository

```bash
git clone -b Public https://github.com/BrenneckeLab/AnnotationPipeline.git
cd AnnotationPipeline
```

### 2. Configure Settings

Edit `settings.txt` with paths appropriate for your environment:

```bash
# Default annotation versions
default_VERSION_dm3="r5.57"
default_VERSION_dm6="r6.40"
default_VERSION_Cel="WS292"

# Path for genome indices and annotation files (permanent storage)
# Do NOT place within the AnnotationPipeline directory
BASE_UTILITY_LOCATION="/path/to/permanent/storage/indices"

# Path for storing demultiplexed library files (permanent storage)
# Do NOT place within the AnnotationPipeline directory
BASE_LIBRARY_STORAGE_LOCATION="/path/to/permanent/storage/libraries"

# Path for results directory (web-accessible for UCSC hub)
# Do NOT place within the AnnotationPipeline directory
BASE_FOLDER="/path/to/results/folder"

# HTTP address for the results directory (for UCSC track hub)
HTTP_PATH="https://your-server.com/results/"

# Path for temporary files (fast storage recommended)
# Do NOT place within the AnnotationPipeline directory
BASE_FOLDER_TMP="/path/to/temp/storage"

# Set to N to skip automatic index creation during setup
RUN_INSTALLATION=Y

# Enable if your cluster uses hyperthreading
HYPER=N
```

#### Custom SLURM Settings

If your cluster requires additional SBATCH specifications, add them after the `#@!@#` marker in `settings.txt`:

```bash
#@!@#
#SBATCH --qos=short
#SBATCH --partition=your_partition
```

### 3. Run Setup

```bash
./setup.sh
```

The setup script will:
- Configure paths throughout the pipeline
- Download required Singularity containers (~2GB)
- Create necessary directory structures
- Copy utility files to the designated locations

### 4. Verify Installation

After setup, run a test with a small dataset to verify the installation is working correctly.

## Quick Start

### Basic Usage

```bash
./annotate_reads.sh -i library_file.txt -t sRNAseq -v dm6 -F "my_experiment"
```

Or using long options:

```bash
./annotate_reads.sh --input library_file.txt --type sRNAseq --genome-version dm6 --folder-name "my_experiment"
```

### Input File Format

Create a tab-separated file listing your libraries:

```
/path/to/library1.fa.gz    sample1_name
/path/to/library2.fa.gz    sample2_name    IP
/path/to/library3.bam      sample3_name    BCi7=ACGTAC BCi5=CGTACT
```

**Columns:**
1. Path to the sequencing file (bam, fa, fq, fa.gz, fq.gz)
2. Sample name (keep short but descriptive)
3. Optional Tags
   - "IP" tag for immunoprecipitation samples in mixed sRNAseq runs - changes normalization for this library to 10M uniquely aligned reads
   - "OX" tag for oxidized small RNA samples - changes normalization for this library to 10M uniquely aligned reads
   - "INVERT" tag to invert orientation of sequences by reverse complementing
   - "BCi7" Barcode sequence for i7 barcodes. [BCi7=######] - triggers demultiplexing
   - "BCi5" Barcode sequence for i5 barcodes. [BCi5=######] - triggers demultiplexing
   - "sRBC" Barcode sequence for small RNA barcodes contained in the 3' adaptor [sRBC=#####] - triggers demultiplexing
   - "SUB" subset library to n-reads [SUB=1000000]
   - "NORM" provide external normalization factor that overwrites other normalization options [NORM=1.5]
   - "COLOR" provide RGB color values for individual libaries [COLOR=100,200,150]
  

### Supported Input Formats

- Local files: `.bam`, `.fa`, `.fq`, `.fa.gz`, `.fq.gz`
- SRA accessions: `SRRxxxxxxx`, `ERRxxxxxxx`, `DRRxxxxxxx`
- Multiple files can be merged using `;` as delimiter

**Example for bam merging:**
```
/groups/lab/data/1a.bam;/groups/lab/data/1b.bam    merged_sample
```

## Command Line Options

The pipeline supports both short (`-x`) and long (`--option`) format options.

### Required Options

| Short | Long | Description |
|-------|------|-------------|
| `-i` | `--input` | Input file containing library information |

### Interactive Options (prompted if not set)

| Short | Long | Description |
|-------|------|-------------|
| `-F` | `--folder-name` | Folder name suffix for results |
| `-t` | `--type` | Library type (sRNAseq, RNAseq, etc.) |
| `-V` | `--version` | FlyBase/WormBase annotation version |

### Genome and Mapping Options

| Short | Long | Description | Default |
|-------|------|-------------|---------|
| `-v` | `--genome-version` | Genome version (dm3, dm6, Cel, ASM) | dm6 |
| `-s` | `--mismatches` | Allowed mismatches for mapping (0-3) | 0 (3 for CLIP) |
| `-Y` | `--y-chrom` | Include Y-chromosome in analysis | Off |
| `-E` | `--random-multi` | Randomly assign multi-mappers | Off |

### Pre-processing Options

| Short | Long | Description | Default |
|-------|------|-------------|---------|
| `-x` | `--force-preprocess` | Force pre-processing even for non-bam input | Off |
| `-N` | `--n-trimm` | Trim N random nucleotides from each end (0-8) | 0 |
| `-A` | `--custom-adaptors` | Use custom adaptor sequences for clipping | Off |
| `-S` | `--single-end` | Treat as single-end data | Off |
| `-2` | `--second-mate` | Use second mate of paired-end data | Off |

### Read Processing Options

| Short | Long | Description | Default |
|-------|------|-------------|---------|
| `-m` | `--min-length` | Minimum read length | 18 (35 for CapSeq) |
| `-M` | `--max-length` | Maximum read length | 35 (sRNA) / 200 (others) |
| `-T` | `--trimm` | Enable read trimming | Off |
| `-f` | `--first` | First base to keep (trimming) | 1 |
| `-l` | `--last` | Last base to keep (trimming) | 1000 |
| | `--polya` | Enable polyA trimming | Off |
| `-I` | `--invert` | Invert read orientation | Off |
| `-J` | `--max-count` | Maximum read count threshold | 1000000 |

### Analysis Options

| Short | Long | Description |
|-------|------|-------------|
| `-D` | `--dge` | Perform differential gene expression analysis |
| `-j` | `--slam` | Enable SLAM-seq analysis mode |
| `-1` | `--ping-pong` | Perform ping-pong and phasing analysis (sRNAseq) |
| `-5` | `--only-5end` | Count only 5' end of reads |
| `-n` | `--no-stranded` | Unstranded analysis mode |
| `-e` | `--extend` | Extend reads to fragment length (CHIPseq) - provide fragment length |
| | `--force-quant-unstranded` | Force unstranded quantification |
| | `--extra-seq` | Path to file containing sequences in fasta-format that will be added to the TE-consensus sequences for histogram generation and bowtie based quantification|

### Normalization Options

| Short | Long | Description |
|-------|------|-------------|
| `-W` | `--wig-norm` | Normalize to 10M uniquely mapping reads |
| `-w` | `--wig-fasta-norm` | Normalize to 10M reads in input fasta |
| `-L` | `--spike-in-norm` | Enable spike-in normalization |
| `-z` | `--no-norm` | Disable track normalization |

### Output Options

| Short | Long | Description |
|-------|------|-------------|
| |`--demux-only` | Only perform demultiplexing |
| `-P` | `--raw-paired` | Generate paired-read output file |
| `-p` | `--only-paired` | Output only paired reads |
| `-Q` | `--fastq-out` | Output adaptor-clipped FASTQ |
| `-q` | `--fastq-out-raw` | Output raw FASTQ |
| `-r` | `--raw` | Output raw (uncollapsed) data |
| `-b` | `--export-bam` | Export collapsed BAM files |
| `-B` | `--export-bam-uncollapsed` | Export uncollapsed BAM files |
| `-3` | `--export-salmon` | Export Salmon quantification files |
| `-G` | `--geo` | Create GEO submission export |
| `-c` | `--color` | Custom track color (R,G,B format) |
| `-a` | `--auto-scale` | Enable auto-scaling for track visualization |

### Advanced Options

| Short | Long | Description |
|-------|------|-------------|
| `-U` | `--user` | Run under different username |
| `-k` | `--keep-tmp` | Keep temporary files |
| `-d` | `--debug` | Enable debug mode |
| `-y` | `--demux-fasta` | Demultiplex FASTA files |

## Annotation Ruleset

When reads map to multiple annotation categories, the pipeline resolves ties  using a priority-based ruleset:

| Priority | Category | Description |
|----------|----------|-------------|
| 1-4 | rRNA / rRNA_AS | Ribosomal RNA (sense/antisense) |
| 5-6 | mito / mito_AS | Mitochondrial sequences |
| 7-8 | tRNA / tRNA_AS | Transfer RNA |
| 9-12 | snoRNA / snRNA | Small nucleolar/nuclear RNA |
| 13-16 | miRNA | MicroRNA and precursors |
| 17-18 | hpRNA | Hairpin RNA |
| 19-21 | exon | mRNA exons and UTRs |
| 22-23 | TE / TE_AS | Transposable elements |
| 24-31 | ncRNA / pseudo | Non-coding RNA and pseudogenes |
| 32-44 | intron | Intronic sequences |
| 99 | none | Unannotated |

The complete ruleset can be found in `utility-files/ruleset.txt`.

## Output Structure

```
results/
├── annotation_counts.txt          # Summary counts per annotation
├── annotation_TE_splitup.txt      # TE category breakdown (DNA, LINE, LTR, Satellite)
├── annotation_mRNA_splitup.txt    # mRNA category breakdown
├── TPM-table_salmon.txt           # Gene quantification (RNAseq)
├── TPM-table_salmon_TE_sense.txt  # TE quantification (sense)
├── TPM-table_salmon_TE_antisense.txt # TE quantification (antisense)
├── sequencing-depth.txt           # Depth estimates
├── read-duplication.txt           # PCR duplicate estimates (sRNAseq)
├── plots/                         # Generated visualizations
│   ├── annotation_counts.pdf      # Stacked bar charts
│   ├── sequencing-depth.pdf       # Depth visualization
│   ├── size-profiles/             # Length distribution plots
|   └── TE-histograms/              # Histogram showing read-coverage for Transposons
├── hub/                           # UCSC track hub
├── log/                           # Processing logs
│   ├── log.txt                    # Main run log with settings
│   ├── chrom.sizes                # Chromosome sizes
│   └── *.log                      # Per-library logs
└── individual-libraries/          # Per-sample outputs
    └── sample_name/
        ├── *_annotated.bed.gz     # Annotated mappings
        ├── *_annotated.fa.gz      # Annotated sequences
        ├── normalization.txt      # Normalization factors
        ├── *_annotation_counts.txt
        ├── *_TE_splitup.txt
        ├── *_mRNA_splitup.txt
        ├── *_sequencing-depth.txt
        ├── *_read-duplication.txt
        ├── size-profiles/         # Length distributions
        ├── wig-files/             # Browser tracks
        ├── cluster-logs/          # Job logs
        └── time-log.txt           # Timing information
```

### Output File Formats

#### Annotated FASTA/BED Name Tag Format

```
NR_1:ATTGCCTCTCATTTTCTCTCCCATATTA:count=1:mapping=u:ann=TE:final_ann=exon:filtered=N:mapcount=1:fine_ann=CDS
```

| Field | Description |
|-------|-------------|
| NR_# | Unique read identifier |
| Sequence | Processed sequence |
| count | Number of identical reads |
| mapping | u=unique, m=multi-mapper, n=not-mapped |
| ann | Raw annotation for this mapping |
| final_ann | Resolved final annotation |
| filtered | Whether filtered from size profiles |
| mapcount | Total mapping locations |
| fine_ann | Detailed annotation category |

#### Normalization File

The `normalization.txt` file contains all computed normalization factors. Allways the last line in the file was used for the normalization of UCSC-tracks, TE-histograms and size-profiles. e.g.:
- Factor for normalization to 10M reads in input fasta (after filtering)
- Factor for normalization to 10M uniquely mapping reads
- Factor for normalization to 1M miRNAs for small RNA libraries

## Examples

### Small RNA-seq Analysis

```bash
# Basic sRNAseq run
./annotate_reads.sh -i samples.txt -t sRNAseq -v dm6 -F "my_experiment"

# With ping-pong analysis
./annotate_reads.sh -i samples.txt -t sRNAseq -v dm6 --ping-pong -F "piRNA_analysis"

# Mixed sRNAseq with IP samples
./annotate_reads.sh -i samples_with_ip.txt -t sRNAseq -v dm6 -F "ip_experiment"

# Custom length range
./annotate_reads.sh -i samples.txt -t sRNAseq -v dm6 -m 20 -M 30 -F "custom_length"
```


## Troubleshooting

### Common Issues

**1. Setup fails with permission errors**
- Ensure write access to all paths specified in `settings.txt`
- Check that the Singularity containers can be downloaded
- Verify the HTTP path is accessible

**2. Jobs fail on cluster**
- Check SLURM logs in the `cluster-logs` directory
- Verify sufficient memory allocation
- Check disk space in temporary directory
- Review `time-log.txt` for timing information

**3. Annotation version not found**
- Run the pipeline without `-V` to interactively select or create annotations
- Check `versions.txt` in the utility location for available versions
- New annotations can be generated automatically

**4. UCSC hub not loading**
- Verify the HTTP path is correctly configured and accessible
- Check that `hub.txt` exists in the results directory
- Ensure proper permissions on the results folder

**5. Script started in background**
- The pipeline must be run in the foreground
- Do not use `&` when starting the script

### Getting Help

- Check the log files in your results directory
- Main log: `log/log.txt` contains run settings and replication command
- Individual sample logs are in `log/*.log`
- Cluster submission logs are in `cluster-logs/`
- RAM usage logs are available for debugging memory issues

## Updating the Pipeline

After installation, use the `update.sh` script to update to newer versions:

```bash
./update.sh
```


## 1kb Tile Analysis

For dm6 genome, the pipeline includes support for 1kb tile analysis. Example R scripts and data files are available in:

```
utility-files/dmel/dm6/1kb_tile/
├── kbWindows_mainchr_genes_dm6.txt
├── kbWindows_mainchr_genes_tss_dm6.txt
├── dm6_1kb_windows_DNAseq_mappability_mainchr_sorted.txt
└── sRNAseq_tile_analyses_example_Nxf3_paper.R
```

## Citation

If you use AnnotationPipeline in your research, please cite:

> Brennecke Lab. AnnotationPipeline: A comprehensive sequencing analysis pipeline.
> https://github.com/BrenneckeLab/AnnotationPipeline


## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- **Lab Website**: [Brennecke Lab at IMBA](https://www.imba.oeaw.ac.at/research/julius-brennecke/)
- **Issues**: [GitHub Issues](https://github.com/BrenneckeLab/AnnotationPipeline/issues)

---
