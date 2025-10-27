# Circe: Peak Coaccessibility Analysis

A Python package for analyzing peak coaccessibility in ATAC-seq data using two different methods: cosine similarity and the circe algorithm.

## Overview

This project provides tools to identify coaccessible peaks in ATAC-seq data, which can help understand chromatin accessibility patterns and regulatory element interactions. The analysis can be performed using either:

1. **Cosine Similarity**: A statistical approach that computes similarity between peak accessibility profiles
2. **Circe Method**: An advanced algorithm that uses the `circe-py` package for more sophisticated coaccessibility analysis

## Features

- **Dual Analysis Methods**: Choose between cosine similarity or circe-based coaccessibility analysis
- **Parallel Processing**: Efficient computation using multiple CPU cores
- **Flexible Filtering**: Filter results based on quantile thresholds or absolute values
- **Window-based Analysis**: Configurable genomic window sizes for local coaccessibility
- **Command-line Interface**: Easy-to-use CLI for batch processing
- **H5AD Support**: Native support for AnnData format commonly used in single-cell genomics

## Installation

### Using Conda (Recommended)

The project includes a conda environment specification. To set up the environment:

```bash
# Create environment from the specification
conda env create -f conda_env_list.txt

# Activate the environment
conda activate circe
```

### Using pip

Alternatively, install dependencies using pip:

```bash
pip install -r requirements.txt
```

### Key Dependencies

- **Python 3.12+**
- **anndata**: For handling single-cell genomics data
- **scanpy**: Single-cell analysis toolkit
- **circe-py**: Advanced coaccessibility analysis
- **scikit-learn**: Machine learning utilities
- **pandas**: Data manipulation
- **numpy**: Numerical computing
- **joblib**: Parallel processing

## Usage

### Command Line Interface

The main script `run_peak_coaccess.py` provides a comprehensive CLI for peak coaccessibility analysis:

#### Basic Usage

```bash
# Basic circe analysis (default method)
python run_peak_coaccess.py input.h5ad -o output.tsv

# Cosine similarity analysis
python run_peak_coaccess.py input.h5ad -m cosine -o output.tsv
```

#### Advanced Options

```bash
# Custom parameters for cosine similarity
python run_peak_coaccess.py input.h5ad -m cosine -t 0.2 -w 500000 -o output.tsv

# Circe analysis with custom window size
python run_peak_coaccess.py input.h5ad -m circe -w 100000 -o output.tsv

# Filter results and save intermediate files
python run_peak_coaccess.py input.h5ad -m cosine -o output.tsv --filter --quantile 0.95 --save-intermediate

# Use specific number of CPU cores
python run_peak_coaccess.py input.h5ad -m circe -j 8 -o output.tsv
```

### Command Line Arguments

| Argument | Short | Description | Default |
|----------|-------|-------------|---------|
| `--input_atac` | `-i` | Input h5ad file containing ATAC-seq data | Required |
| `--output` | `-o` | Output TSV file for coaccessibility results | Required |
| `--method` | `-m` | Analysis method: `cosine` or `circe` | `circe` |
| `--window` | `-w` | Window extension size in base pairs | `250000` |
| `--jobs` | `-j` | Number of parallel jobs (-1 for all cores) | `-1` |
| `--filter` | | Filter results based on threshold | False |
| `--quantile` | | Quantile threshold for filtering (0.0-1.0) | `0.9` |
| `--filter-threshold` | | Absolute threshold for filtering | None |
| `--verbose` | | Verbosity level (0-10) | `5` |
| `--save-intermediate` | | Save unfiltered results | False |

### Input Data Requirements

Your input h5ad file should contain:

- **ATAC-seq accessibility matrix**: Cells × Peaks
- **Peak annotations**: Must include chromosome (`chr`), start, and end coordinates
- **Standard AnnData format**: Compatible with scanpy/anndata

Example peak annotation format:
```
peak_id    chr    start    end
chr1_start_end     chr1   1000     2000
chr2:start-end     chr1   5000     6000
```

### Output Format

The analysis produces a TSV file with the following columns:

- **Peak1**: First peak identifier
- **Peak2**: Second peak identifier  
- **coaccess**: Coaccessibility score

Additionally, a `peak_names.txt` file is created containing all peak identifiers.

## Methods

### Cosine Similarity Method

The cosine similarity approach computes the similarity between peak accessibility profiles across cells. Key features:

- **Local Analysis**: Only compares peaks within a specified genomic window
- **Parallel Processing**: Processes each chromosome independently
- **Threshold Filtering**: Filters connections based on similarity scores
- **Overlap Detection**: Accounts for peak overlaps in the analysis

### Circe Method

The circe method uses the `circe-py` [package](https://github.com/cantinilab/Circe) for advanced coaccessibility analysis:

- **Advanced Algorithm**: More sophisticated statistical modeling
- **Configurable Parameters**: Multiple tuning parameters for optimization
- **Distance Constraints**: Optional distance-based filtering
- **Sampling-based**: Uses statistical sampling for robust results

## Examples

### Example 1: Basic Analysis

```bash
# Run circe analysis on ATAC-seq data
python run_peak_coaccess.py atac_data.h5ad -o coaccess_results.tsv
```

### Example 2: Cosine Similarity with Custom Parameters

```bash
# Use cosine similarity with higher threshold and larger window
python run_peak_coaccess.py atac_data.h5ad -m cosine -t 0.3 -w 500000 -o cosine_results.tsv
```

### Example 3: Filtered Results

```bash
# Filter to top 5% of connections
python run_peak_coaccess.py atac_data.h5ad -m circe --filter --quantile 0.95 -o filtered_results.tsv
```

## Performance Considerations

- **Memory Usage**: Large datasets may require significant RAM
- **CPU Cores**: Use `-j` parameter to control parallel processing
- **Window Size**: Larger windows increase computation time but may capture more distant interactions
- **Filtering**: Early filtering can reduce memory usage and processing time

## Troubleshooting

### Common Issues

1. **Missing chromosome annotations**: Ensure your h5ad file has a `chr` column in `adata.var`
2. **Memory errors**: Reduce the number of parallel jobs or use smaller window sizes
3. **File format errors**: Verify your input file is a valid h5ad format

### Error Messages

- `'chr' column not found`: Add chromosome annotations to your peak metadata
- `Input file does not exist`: Check the file path and permissions
- `Error loading data`: Verify the h5ad file is not corrupted

## Contact

For questions or issues, please contact the project maintainers or create an issue in the project repository.
