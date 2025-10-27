"""
CLI for Peak Coaccessibility Analysis

This script provides a command-line interface for analyzing peak coaccessibility
using either cosine similarity or the circe method.
"""

import argparse
import sys
import os
from pathlib import Path
import pandas as pd
import anndata as ad
import scanpy as sc

from process_peak_coaccessibility import (
    compute_cosine_sim_peak_coaccess, 
    compute_circe_peak_coaccess, 
    filter_coaccess_peaks
)


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Analyze peak coaccessibility using cosine similarity or circe method",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic cosine similarity analysis
  python run_peak_coaccess.py input.h5ad -o output.tsv

  # Use circe method with custom parameters
  python run_peak_coaccess.py input.h5ad -m circe -o output.tsv

  # Cosine similarity with custom threshold and window
  python run_peak_coaccess.py input.h5ad -m cosine -t 0.2 -w 500000 -o output.tsv

  # Filter results and save
  python run_peak_coaccess.py input.h5ad -m cosine -o output.tsv --filter --quantile 0.95
        """
    )
    
    # Required arguments
    parser.add_argument(
        '--input_atac', '-i',
        help='Input h5ad file containing ATAC-seq data'
    )
    parser.add_argument(
        '-o', '--output',
        required=True,
        help='Output TSV file for coaccessibility results'
    )
    
    # Method selection
    parser.add_argument(
        '-m', '--method',
        choices=['cosine', 'circe'],
        default='circe',
        help='Method to use for coaccessibility analysis (default: circe)'
    )
    

    parser.add_argument(
        '-w', '--window',
        type=int,
        default=250000,
        help='Window extension size in base pairs (default: 250000)'
    )
    parser.add_argument(
        '-j', '--jobs',
        type=int,
        default=2,
        help='Number of parallel jobs default: 2'
    )
    
    # Filtering options
    parser.add_argument(
        '--filter',
        action='store_true',
        help='Filter results based on threshold'
    )
    parser.add_argument(
        '--quantile',
        type=float,
        required=False,
        default=0.9,
        help='Quantile threshold for filtering (default: 0.9)'
    )
    parser.add_argument(
        '--filter-threshold',
        type=float,
        required=False,
        default=None,
        help='Absolute threshold for filtering (overrides quantile)'
    )
    
    # Output options
    parser.add_argument(
        '--verbose',
        type=int,
        default=5,
        help='Verbosity level (default: 5)'
    )
    parser.add_argument(
        '--save-intermediate',
        action='store_true',
        help='Save intermediate results before filtering'
    )
    
    return parser.parse_args()


def validate_input_file(input_file):
    """Validate that the input file exists and is readable."""
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist.", file=sys.stderr)
        sys.exit(1)
    
    if not input_file.endswith(('.h5ad', '.h5')):
        print(f"Warning: Input file '{input_file}' may not be in h5ad format.", file=sys.stderr)


def load_atac_data(input_file):
    """Load ATAC-seq data from h5ad file."""
    try:
        print(f"Loading data from {input_file}...")
        adata = sc.read_h5ad(input_file)
        print(f"Loaded data with {adata.n_obs} cells and {adata.n_vars} peaks")
        
        # Check for required columns
        if 'chr' not in adata.var.columns:
            print("Error: 'chr' column not found in var. Please ensure your data has chromosome annotations.", file=sys.stderr)
            sys.exit(1)
            
        return adata
    except Exception as e:
        print(f"Error loading data: {e}", file=sys.stderr)
        sys.exit(1)


def run_cosine_analysis(adata, args):
    """Run cosine similarity analysis."""
    print("Running cosine similarity analysis...")
    print(f"Parameters: threshold={args.filter_threshold}, window={args.window}, jobs={args.jobs}")
    
    coaccess_df = compute_cosine_sim_peak_coaccess(
        atac_data=adata,
        threshold=args.filter_threshold,
        window_ext=args.window,
        n_jobs=args.jobs,
        verbose=args.verbose
    )
    
    return coaccess_df


def run_circe_analysis(adata, args):
    """Run circe analysis."""
    print("Running circe analysis...")
    
    coaccess_df = compute_circe_peak_coaccess(
        atac_data=adata,
        njobs=args.jobs,
        window_size=args.window
        )
    
    return coaccess_df


def filter_results(coaccess_df, filter_threshold, quantile):
    """Filter coaccessibility results."""
    
    filtered_df = filter_coaccess_peaks(coaccess_df, threshold=filter_threshold, quantile_threshold=quantile)
    print(f"Filtered from {len(coaccess_df)} to {len(filtered_df)} connections")
    
    return filtered_df


def save_results(coaccess_df, output_file, args):
    """Save results to output file."""
    print(f"Saving results to {output_file}...")
    
    # Create output directory if it doesn't exist
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    # Save results
    coaccess_df.to_csv(output_file, index=False, sep='\t')
    print(f"Saved {len(coaccess_df)} coaccessibility connections")
    
    # Save intermediate results if requested
    if args.save_intermediate and args.filter:
        intermediate_file = output_file.replace('.tsv', '_unfiltered.tsv')
        print(f"Saving unfiltered results to {intermediate_file}...")
        # Note: This would require storing the unfiltered results


def main():
    """Main function."""
    args = parse_arguments()
    
    # Validate input
    validate_input_file(args.input_atac)
    
    # Load data
    adata = load_atac_data(args.input_atac)
    
    # Run analysis based on method
    if args.method == 'cosine':
        coaccess_df = run_cosine_analysis(adata, args)
    elif args.method == 'circe':
        coaccess_df = run_circe_analysis(adata, args)
    
    print(f"Found {len(coaccess_df)} coaccessibility connections")
    
    # Filter results if requested
    if args.filter:
        coaccess_df = filter_results(coaccess_df, args.filter_threshold, args.quantile)
    
    # Save results
    save_results(coaccess_df, args.output, args)

    peak_names_file_out = Path(args.output).parent / "peak_names.txt"

    with open(peak_names_file_out, "w") as fh:
        fh.write('\n'.join(list(adata.var_names)))
    
    
    print("Analysis complete!")


if __name__ == "__main__":
    main()




