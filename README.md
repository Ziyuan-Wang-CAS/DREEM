# DREEM: Detection of RNA folding Ensembles using Expectation-Maximization clustering

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

DREEM (Detection of RNA folding Ensembles using Expectation-Maximization clustering) is a computational tool designed to detect alternative structures formed by RNA molecules using single-molecule chemical probing data. The software leverages Expectation-Maximization (EM) clustering algorithms to identify diverse RNA folding conformations.

**This is an enhanced version based on the original DREEM from [Code Ocean Capsule 6175523](https://codeocean.com/capsule/6175523/tree/v1)**, with significant improvements for analyzing live-cell RNA structures using proximity labeling techniques.

### Key Features

#### Core Capabilities
- **RNA Structural Diversity Detection**: Identifies multiple folding conformations that RNA molecules can adopt
- **Single-Molecule Analysis**: High-precision analysis based on single-molecule chemical probing data
- **Clustering Analysis**: Uses EM algorithms to cluster similar RNA structures
- **Visualization Output**: Generates comprehensive interactive HTML plots and statistical reports

#### Enhanced Features (This Version)
- **Customizable Base Type Analysis**: Configurable selection of nucleotide types for mutation rate analysis (e.g., excluding A and C bases)
- **Control-Experimental Comparison Pipeline**: Advanced denoising algorithm using Fisher's exact test and Benjamini-Hochberg correction
- **Low Signal-to-Noise Ratio Optimization**: Specifically optimized for proximity labeling sequencing data with low mutation rates
- **Live-Cell RNA Structure Analysis**: Designed for analyzing RNA secondary structures and clustering in living cells
- **End-to-End Automation**: Complete pipeline from FASTQ files to final clustering results with robust error handling

### Enhanced Analysis Pipeline

This version implements an improved four-step analysis workflow:

1. **Mapping**: Aligns sequencing reads to a reference genome using Bowtie2
2. **Bit Vector Conversion**: Converts sequencing reads into binary bit vector representations with customizable nucleotide selection
3. **Statistical Denoising**: Compares experimental and control samples using Fisher's exact test and multiple testing correction
4. **EM Clustering**: Applies Expectation-Maximization algorithms to cluster filtered, high-confidence bit vectors

#### Statistical Denoising Algorithm
- **Fisher's Exact Test**: Statistical comparison between experimental and control mutation rates at each position
- **Multiple Testing Correction**: Benjamini-Hochberg procedure to control false discovery rate
- **Fold-Change Filtering**: Additional filtering based on mutation rate fold-changes
- **Signal Enhancement**: Removes background noise to improve clustering accuracy for low SNR data

## Installation Requirements

### System Requirements
- Linux/macOS operating system
- Python 3.6 or higher
- At least 8GB RAM
- Multi-core CPU (recommended)

### External Dependencies
Ensure the following external tools are installed and added to your PATH environment variable:
- **Bowtie2**: Sequence alignment tool
- **FastQC**: Sequencing quality control
- **TrimGalore**: Adapter and quality trimming
- **Picard**: BAM file quality assessment

### Python Dependencies
Key Python packages required:
- numpy
- scipy
- pandas
- biopython
- plotly
- statsmodels
- matplotlib

## Environment Setup

### Creating Conda Environment

We provide a complete conda environment configuration with all exact package versions:

```bash
# Create conda environment
conda env create -f environment.yml

# Activate environment
conda activate dreem
```

### Installation Verification
```bash
# Check Python packages
python -c "import numpy, scipy, pandas, Bio; print('Python packages OK')"

# Check external tools
bowtie2 --version
fastqc --version
trim_galore --version
java -jar code/picard.jar -h
```

## Usage Guide

### Data Preparation

1. **Reference Sequence**: Prepare FASTA format reference sequence file
2. **Sequencing Data**: FASTQ format sequencing files (supports single-end and paired-end)
3. **Aligned Data**: Or provide pre-aligned BAM files

### Directory Structure
Recommended data organization:
```
project/
├── data/
│   ├── sample_name/
│   │   ├── sample_R1.fastq.gz  # Forward reads
│   │   ├── sample_R2.fastq.gz  # Reverse reads (paired-end)
│   │   └── reference.fasta     # Reference sequence
├── code/                       # DREEM program files
└── results/                    # Output results
```

### Basic Usage

#### Method 1: Main Program Script (Recommended)

```bash
cd DREEM/code
python Run_DREEM.py [input_dir] [output_dir] [sample_name] [ref_name] [start_pos] [end_pos] [options]
```

**Parameter Description:**
- `input_dir`: Directory path containing input files
- `output_dir`: Results output directory path
- `sample_name`: Name of the analysis sample
- `ref_name`: Name of the reference sequence
- `start_pos`: Start position of analysis region (1-based)
- `end_pos`: End position of analysis region (1-based)

**Optional Parameters:**
- `--single`: Single-end sequencing data
- `--struct`: Enable secondary structure prediction
- `--fastq`: Start analysis from FASTQ files (no BAM file)
- `--ctrl`: Control sample

#### Method 2: Using Shell Scripts (Recommended for Routine Analysis)

We provide several shell scripts for different analysis scenarios:

##### 2.1 Single Sample Analysis (`DREEM_Run.sh`)
This script is designed for analyzing a single sample with predefined parameters.

```bash
# Edit configuration parameters
vim DREEM_Run.sh

# Modify the following parameters in the script:
# REF_NAME=your_reference_name
# POS_START=start_position  
# POS_END=end_position
# DATA_NAME=your_sample_name
# FASTQ=yes/no (whether to start from FASTQ files)
# CTRL=yes/no (whether this is a control sample)
# STRUCT=yes/no (whether to perform structure prediction)

# Run analysis
bash DREEM_Run.sh
```

##### 2.2 Parallel Analysis for Large Regions (`DREEM_Submit.sh` + `DREEM_Parallel_Run.sh`)
For analyzing large genomic regions, these scripts enable parallel processing by splitting the region into chunks.

**Step 1**: Configure the submission script
```bash
# Edit DREEM_Submit.sh
vim DREEM_Submit.sh

# Set these parameters:
REF_NAME=your_reference_name       # Reference sequence name
POS_START=start_position           # Analysis start position
POS_END=end_position              # Analysis end position  
DATA_NAME=your_sample_name        # Sample name
FASTQ=yes/no                      # Start from FASTQ files
CTRL=yes/no                       # Control sample flag
STRUCT=yes/no                     # Structure prediction flag
CHUNK_SIZE=100                    # Size of each parallel chunk
```

**Step 2**: Submit parallel jobs (requires SGE/Grid Engine)
```bash
# Submit jobs to cluster
bash DREEM_Submit.sh
```

This will:
- Calculate the number of tasks needed based on region size and chunk size
- Submit an array job to SGE cluster system
- Each task analyzes a specific chunk of the region
- Results are saved in separate directories for each chunk

##### 2.3 Copying Statistical Results (`Copy_stats_json.sh`)
This script copies statistical JSON files between experiment and control results, useful for comparative analysis.

```bash
# Edit Copy_stats_json.sh
vim Copy_stats_json.sh

# Set the result directory names:
EXP_RESULTS=results_YYYYMMDD_experiment_name
CONTROL_RESULTS=results_YYYYMMDD_control_name

# Run the script
bash Copy_stats_json.sh
```

### Script Workflow and Features

#### Shell Script Advantages
The provided shell scripts offer several advantages over manual execution:

1. **Automated Environment Setup**: Scripts automatically activate the conda environment
2. **Robust Error Handling**: Built-in error checking and logging
3. **Parallel Processing**: Support for large-scale analysis via cluster systems
4. **Automatic Directory Management**: Creates output directories automatically
5. **Real-time Monitoring**: Provides process IDs and log file locations for monitoring
6. **Date-stamped Results**: Automatically includes dates in result directory names

#### Cluster System Integration
The parallel scripts are designed for SGE (Sun Grid Engine) cluster systems:
- `DREEM_Submit.sh`: Configures and submits array jobs
- `DREEM_Parallel_Run.sh`: Executes individual tasks in the array
- Supports automatic task distribution and load balancing
- Each parallel task processes a specific genomic region chunk

#### Script Configuration Details

**DREEM_Run.sh Configuration Parameters:**
```bash
REF_NAME=human_mitochondrial_genome-NCBI_NC_012920.1  # Reference name
POS_START=12337                                        # Start position  
POS_END=14148                                          # End position
DATA_NAME=MTS_s4U_APDF20_0906_rep1                   # Sample name
FASTQ=no                                              # Process FASTQ files
CTRL=no                                               # Control sample
STRUCT=yes                                            # Structure prediction
```

**Parallel Processing Parameters:**
```bash
CHUNK_SIZE=100                    # Nucleotides per parallel task
HOSTNAME=fnode006                 # Cluster node specification  
QUEUE=b1.q                       # SGE queue name
```

### Usage Examples

#### Example 1: Complete Pipeline Analysis (from FASTQ)
```bash
python Run_DREEM.py \
    ../data/sample1 \
    ../results/sample1_analysis \
    sample1 \
    mt_genome \
    12337 \
    14148 \
    --fastq \
    --struct
```

#### Example 2: Analysis from BAM Files
```bash
python Run_DREEM.py \
    ../data/sample1 \
    ../results/sample1_analysis \
    sample1 \
    mt_genome \
    12337 \
    14148 \
    --struct
```

#### Example 3: Single-End Sequencing Analysis
```bash
python Run_DREEM.py \
    ../data/sample1 \
    ../results/sample1_analysis \
    sample1 \
    mt_genome \
    12337 \
    14148 \
    --single \
    --fastq
```

#### Example 4: Using Shell Scripts for Routine Analysis

**Single Sample Analysis:**
```bash
# Edit parameters in DREEM_Run.sh
vim DREEM_Run.sh

# Run the analysis
bash DREEM_Run.sh

# Monitor progress
tail -f logs/YYYYMMDD_sample_name_position.log
```

**Parallel Analysis for Large Regions:**
```bash
# Configure parallel analysis
vim DREEM_Submit.sh
# Set: REF_NAME, POS_START, POS_END, DATA_NAME, CHUNK_SIZE

# Submit parallel jobs
bash DREEM_Submit.sh

# Monitor job status (SGE commands)
qstat -u $USER                    # Check job status
qdel JOB_ID                       # Cancel jobs if needed

# Check individual task logs
ls qsub_logs/DREEM/YYYYMMDD_sample_position/
tail -f qsub_logs/DREEM/YYYYMMDD_sample_position/DREEM_Parallel_Run.sh.o*
```

**Control-Experimental Analysis Workflow (Enhanced Feature):**
```bash
# 1. Process control sample first (generates background statistics)
vim DREEM_Run.sh  # Set CTRL=yes, DATA_NAME=control_sample
bash DREEM_Run.sh

# 2. Process experimental sample (applies denoising algorithm)  
vim DREEM_Run.sh  # Set CTRL=no, DATA_NAME=experimental_sample
bash DREEM_Run.sh

# 3. Copy control statistics for comparative analysis
vim Copy_stats_json.sh
# Set: EXP_RESULTS and CONTROL_RESULTS directory names
bash Copy_stats_json.sh

# The pipeline automatically:
# - Generates control_group_stats.json from control sample
# - Applies Fisher's exact test for experimental vs control comparison
# - Performs multiple testing correction (Benjamini-Hochberg)
# - Filters positions based on statistical significance and fold-change
# - Produces denoised bit vectors for improved clustering
```

### Step-by-Step Execution

If you need to run the analysis step by step, you can execute each step individually:

#### Step 1: Sequence Mapping
```bash
python Mapping.py [sample_name] [ref_name] [paired] [cpus] [seed_length] [max_fragment_length] [input_dir] [output_dir] [picard_path]
```

#### Step 2: Bit Vector Generation
```bash
python BitVector.py [sample_name] [ref_name] [start_pos] [end_pos] [surrounding_bases] [qscore_file] [qscore_cutoff] [input_dir] [output_dir] [paired] [picard_path] [from_fastq]
```

#### Step 3: EM Clustering
```bash
python EM_Clustering.py [sample_name] [ref_name] [start_pos] [end_pos] [min_iterations] [info_threshold] [conv_cutoff] [num_runs] [max_k] [cpus] [norm_perc_bases] [exclude_ac] [signal_threshold] [struct_prediction] [input_dir] [output_dir] [control]
```

## Configuration Parameters

### Main Parameter Descriptions

#### Mapping Parameters
- `CPUS = 24`: Number of threads to use
- `L = 12`: Bowtie2 seed length
- `X = 1000`: Maximum fragment length for paired-end sequencing

#### Bit Vector Parameters
- `SUR_BASES = 10`: Number of bases surrounding deletions
- `QSCORE_CUTOFF = 20`: Quality score threshold for valid bases

#### Clustering Parameters
- `MIN_ITS = 300`: Minimum iterations per EM run
- `INFO_THRESH = 1.0`: Threshold for informative sites
- `CONV_CUTOFF = 0.5`: Log-likelihood difference for convergence
- `NUM_RUNS = 10`: Number of independent EM runs per K
- `MAX_K = 3`: Maximum number of clusters K
- `SIG_THRESH = 0.005`: Threshold for signal-to-noise distinction
- `NORM_PERC_BASES = 10`: Percentage of bases for normalization
- `exc_AC = True`: Whether to exclude A and C bases from analysis (customizable base type selection)

#### Enhanced Denoising Parameters
- **Fisher's Exact Test**: Statistical test for experimental vs. control comparison
- **q-value threshold**: Multiple testing corrected p-value cutoff (typically 0.05)
- **Fold-change threshold**: Minimum fold-change between experimental and control mutation rates
- **Background noise filtering**: Automated removal of positions with non-significant differences

## Output Results

### File Structure
After analysis completion, the output directory contains:

```
results/
├── mapping/                    # Mapping results
│   ├── *.bam                  # Alignment files
│   ├── *_stats.txt            # Mapping statistics
│   └── plots/                 # Quality control plots
├── bitvectors/                # Bit vectors
│   ├── *.txt                  # Bit vector files
│   ├── *_read_coverage.html   # Interactive coverage plots
│   ├── *_pop_avg.html         # Population average plots
│   ├── *_DMS_mutations.html   # DMS mutation plots
│   └── *_mutation_histogram.html # Mutation histogram plots
└── clustering/                # Clustering results
    ├── *_clusters.json        # Clustering results
    ├── LogLikes_Iterations.html # EM convergence plots
    ├── DMSModRate.html         # Modification rate plots
    ├── DMSModRate_Clusters.html # Cluster-specific mod rates
    ├── control_mut_pop_avg.html # Control sample plots (if available)
    ├── experimental_mut_pop_avg.html # Experimental sample plots
    └── scatter_plots/          # Cluster scatter plots
        ├── *_normmus.html      # Normalized cluster plots
        └── *_mus.html          # Raw cluster plots
```

### Main Output File Descriptions

#### 1. Bit Vector Files (*.txt)
- Each line represents a bit vector for one read pair
- Binary format representing modification status at each position

#### 2. Clustering Result Files (*.json)
- Contains detailed information for each cluster
- Cluster centers, member reads, statistical parameters, etc.

#### 3. Interactive Visualization Plots (HTML Format)
All plots are generated as interactive HTML files using Plotly, allowing for zooming, hovering, and data exploration:

- **Coverage Plots** (`*_read_coverage.html`): Display read coverage at each position
- **Population Average Plots** (`*_pop_avg.html`): Show population-level modification rates
- **DMS Mutation Plots** (`*_DMS_mutations.html`): Visualize DMS modification patterns
- **Mutation Histograms** (`*_mutation_histogram.html`): Distribution of mutations per read
- **EM Convergence Plots** (`LogLikes_Iterations.html`): Show log-likelihood convergence during EM clustering
- **Modification Rate Plots** (`DMSModRate.html`, `DMSModRate_Clusters.html`): Modification rates overall and by cluster
- **Cluster Scatter Plots** (`*_normmus.html`, `*_mus.html`): 2D projections of different clusters
- **Control/Experimental Comparison** (`control_mut_pop_avg.html`, `experimental_mut_pop_avg.html`): Side-by-side comparison plots

**Viewing Interactive Plots:**
All HTML plots can be opened in any web browser and offer interactive features:
```bash
# Open plots in browser
firefox results/bitvectors/sample_pop_avg.html
# or
open results/clustering/DMSModRate.html  # macOS
# or double-click the HTML files in file explorer
```

## Performance Optimization

### Hardware Recommendations
- **CPU**: Multi-core processor, 16+ cores recommended
- **Memory**: At least 16GB, 32GB+ recommended
- **Storage**: SSD hard drive, sufficient temporary space

### Software Optimization
- Adjust thread count based on available memory
- Use parallel analysis for multiple samples
- Regularly clean temporary files

## Citation

If you use this enhanced version of DREEM in your research, please cite:

**Original DREEM method:**
> Tomezsko, P.J., Corbin, V.D.A., Gupta, P. et al. Determination of RNA structural diversity and its role in HIV-1 RNA splicing. Nature 582, 438–442 (2020).

**Enhanced version source:**
> This implementation is based on the original DREEM code from Code Ocean Capsule 6175523 (https://codeocean.com/capsule/6175523/tree/v1) with significant enhancements for live-cell proximity labeling data analysis.

**Enhanced version development:**
> Wang, Z. Enhanced DREEM implementation with statistical denoising and live-cell RNA structure analysis capabilities. ShuoHan Lab, CEMCS (2025).

## License

This enhanced version maintains compatibility with the original MIT License from the base DREEM implementation, with additional enhancements developed by:

**Copyright (c) 2025 ShuoHan Lab, CEMCS**  
Enhanced features including statistical denoising algorithms and live-cell analysis capabilities.

The original DREEM components remain under their respective MIT License. See LICENSE file for complete details.

## Contact

### Original DREEM
- **Corresponding Author**: Silvi Rouskin (srouskin@wi.mit.edu)
- **Original Code**: [Code Ocean Capsule 6175523](https://codeocean.com/capsule/6175523/tree/v1)

### Enhanced Version
- **Developer**: Ziyuan Wang (wangziyuan2024@sibcb.ac.cn)
- **Institution**: ShuoHan Lab, Center for Excellence in Molecular Cell Science (CEMCS)
- **GitHub**: [https://github.com/Ziyuan-Wang-CAS/DREEM]
- **For enhancement-related issues**: Please contact the developer directly or through GitHub issues

## Changelog

### v1.0.0
- Initial release
- Complete DREEM analysis pipeline support

---
