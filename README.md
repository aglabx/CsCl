# CsCl - Virtual Centrifuge for Metagenomic Data Analysis

CsCl is a computational tool that emulates DNA ultracentrifugation in cesium chloride (CsCl) gradients, providing a modern conceptual analogue to this classical technique from molecular biology. By calculating GC content and k-mer frequency abundance for each sequencing read, CsCl generates visualizations reminiscent of ultracentrifugation banding patterns, enabling exploration of genome composition and structure.

## Overview

Named after cesium chloride used in physical density gradient centrifugation, CsCl offers a computational alternative to traditional ultracentrifugation techniques, bridging historically established conceptual frameworks with modern high-throughput sequencing analysis. By plotting GC content against k-mer frequency, it creates a visualization similar to a density gradient that can help distinguish different genomic fractions or microbial populations in a sample.

### Key Features

- Fast k-mer counting using the KMC library for efficient processing of large datasets
- Calculation of GC content and abundance score (median k-mer frequency) for each read
- Visualization of genomic fractions based on their compositional and abundance characteristics
- Support for direct k-mer counting or using pre-computed KMC databases
- Automatic generation of Python plotting scripts for 2D representation

## Theoretical Background

CsCl builds upon two fundamental genomic properties:

1. **GC Content**: The percentage of guanine and cytosine bases in DNA directly correlates with buoyant density in traditional cesium chloride gradients. Different genomic regions and organellar DNA often have characteristic GC content.

2. **K-mer Abundance**: The frequency of short nucleotide sequences of length k serves as a robust proxy for copy number variation within a genome. Highly repetitive genomic regions exhibit higher k-mer frequencies compared to unique, single-copy regions.

The computational 2D space defined by these two properties aims to resolve distinct DNA fractions similar to physical banding patterns observed in traditional ultracentrifugation.

## Requirements

- C++17 compatible compiler (e.g., GCC 7+ or Clang 5+)
- KMC (k-mer counter) installed and available in PATH
- Python 3 with pandas, matplotlib, and numpy (for visualization)

## Installation

```bash
# Clone the repository
git clone https://github.com/yourusername/CsCl.git
cd CsCl

# Compile
make

# Optional: Install system-wide
sudo make install
```

### Debug Build

For debugging purposes:

```bash
make clean
make DEBUG=1
```

## Usage

Basic usage:

```bash
./CsCl <output_prefix> -k <kmer_size> <input_fastq1> [input_fastq2 ...]
```

Options:

- `<output_prefix>`: Prefix for output files (_gc_abundance.tsv, _kmc_db*, etc.)
- `-k <kmer_size>`: K-mer size (required if KMC needs to run, ignored if --use-kmc-db is used)
- `--min-count <min_count>`: Minimum k-mer count threshold for KMC (default: 2)
- `--use-kmc-db <prefix>`: Use existing KMC database files (<prefix>.kmc_pre, <prefix>.kmc_suf)
- `<input_fastq...>`: One or more input FASTQ files

## Examples

### Basic Analysis

```bash
./CsCl my_sample -k 31 sample1.fastq sample2.fastq
```

### Use Existing KMC Database

```bash
./CsCl my_sample --use-kmc-db existing_db sample1.fastq sample2.fastq
```

### Visualization

After running CsCl, you can generate plots using the auto-generated Python script:

```bash
python my_sample_plot.py
# or
./my_sample_plot.py
```

## Output

CsCl produces several output files:

1. `<prefix>_gc_abundance.tsv`: Tab-separated file containing ReadID, GC_Content, and Median_Kmer_Frequency
2. `<prefix>_kmc_db.kmc_pre` and `<prefix>_kmc_db.kmc_suf`: KMC database files
3. `<prefix>_plot.py`: Python script for visualization
4. Generated plots when running the Python script:
   - `<prefix>_hexbin_plot.png`: Hexbin plot for large datasets
   - `<prefix>_scatter_plot.png`: Scatter plot for smaller datasets

## How It Works

1. **Phase 1**: Counts k-mers in the input FASTQ files using KMC
2. **Phase 2**: Processes each read to calculate GC content and median k-mer frequency
3. **Phase 3**: Outputs results and generates plotting scripts

## Applications

### Genome Characterization
CsCl provides a valuable approach for exploring genome organization and identifying genomic fractions with distinct compositional and abundance characteristics, particularly useful for initial exploration of large genomic datasets.

### Metagenomics
In complex metagenomic datasets containing DNA from numerous microbial species, CsCl can provide an overview of the community's compositional landscape. Different microbial species often exhibit distinct GC content ranges and varying abundance levels, potentially forming separate clusters in the visualization.

### Cancer Genomics
CsCl can help visualize significant alterations such as aneuploidy and copy number variations in cancer genomes. Regions with copy number gains would likely be represented by higher k-mer abundance scores and may appear as distinct clusters in the 2D plot.

### Study of Repetitive Elements
Different families of repetitive elements (satellite DNA, LINEs, SINEs, LTR retrotransposons) often possess distinct ranges of GC content and copy numbers, potentially forming distinguishable clusters in the CsCl plot.

## Comparison with Related Tools

While tools like KAT (K-mer Analysis Toolkit) also utilize k-mer frequencies and GC composition, CsCl is distinguished by:

1. Operating at the read level rather than the k-mer level
2. Specifically emulating the classic DNA ultracentrifugation technique
3. Focusing on visualization that connects to historically established conceptual frameworks
4. Combining both GC content and abundance information in an intuitive 2D representation

## Citation

If you use CsCl in your research, please cite:

[Citation information: "Enhancing Computational Analysis of DNA Ultracentrifugation through Integration with Existing Bioinformatics Approaches"]

## License

[License information]

## Contact

[Your contact information]