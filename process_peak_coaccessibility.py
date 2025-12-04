"""
Peak Coaccessibility Analysis using Cosine Similarity

This module provides functionality to analyze peak coaccessibility using cosine similarity
with parallel processing for efficiency.
"""

import numpy as np
import pandas as pd
import anndata as ad
import circe as ci
import os
import matplotlib.pyplot as plt
from sklearn.metrics.pairwise import cosine_similarity
from joblib import Parallel, delayed
from typing import Tuple, List, Union, Optional


def find_local_connections(X_chr: np.ndarray, peaks_chr: pd.DataFrame, 
                          threshold: float, window_ext: int) -> List[Tuple[str, str, float]]:
    """
    Find local connections between peaks using cosine similarity.
    
    Args:
        X_chr: ATAC-seq data matrix for a chromosome (peaks x cells)
        peaks_chr: DataFrame with peak information (must have 'start', 'end' columns)
        threshold: Minimum cosine similarity threshold for connections
        window_ext: Window extension size for local analysis
        
    Returns:
        List of tuples containing (peak1, peak2, similarity_score)
    """
    starts = peaks_chr["start"].to_numpy()
    ends = peaks_chr["end"].to_numpy()
    names = peaks_chr.index.to_numpy()
    out = []

    for i in range(len(peaks_chr)):
        # find local window
        left = np.searchsorted(starts, starts[i] - window_ext, side='left')
        right = np.searchsorted(starts, ends[i] + window_ext, side='right')

        # compute similarities only in that slice
        if right - left <= 1: 
            continue
        sim_block = cosine_similarity(X_chr[i:i+1], X_chr[left:right]).ravel()

        for j, s in zip(range(left, right), sim_block):
            if j <= i: 
                continue
            # check overlap and threshold
            if (starts[i] - window_ext < ends[j]) and (starts[j] - window_ext < ends[i]):
                if abs(s) >= threshold:
                    out.append((names[i], names[j], s))
    return out


def compute_cosine_sim_peak_coaccess(atac_data: ad.AnnData, 
                                threshold: float = 0.1, 
                                window_ext: int = 250_000,
                                n_jobs: int = -1,
                                verbose: int = 5) -> pd.DataFrame:
    """
    Compute peak coaccessibility using cosine similarity with parallel processing.
    
    Args:
        atac_data: AnnData object containing ATAC-seq data with peak annotations
        threshold: Minimum cosine similarity threshold for connections (default: 0.1)
        window_ext: Window extension size for local analysis in base pairs (default: 250,000)
        n_jobs: Number of parallel jobs (-1 for all cores, default: -1)
        verbose: Verbosity level for parallel processing (default: 5)
        
    Returns:
        DataFrame with columns ['Peak1', 'Peak2', 'coaccess'] containing peak pairs
        and their coaccessibility scores
        
    Example:
        >>> import scanpy as sc
        >>> adata = sc.read_h5ad('atac_data.h5ad')
        >>> coaccess_df = compute_peak_coaccessibility(adata, threshold=0.2, window_ext=500000)
    """
    # Process each chromosome in parallel
    connections = Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(find_local_connections)(
            atac_data[:, atac_data.var["chr"] == chrom].X.T,
            atac_data.var[atac_data.var["chr"] == chrom],
            threshold,
            window_ext
        )
        for chrom in atac_data.var["chr"].unique()
    )
    
    # Flatten results and create DataFrame
    coaccess_network = pd.DataFrame([x for sub in connections for x in sub],
                      columns=["Peak1", "Peak2", "coaccess"])
    
    return coaccess_network


def compute_circe_peak_coaccess(
    atac_data: ad.AnnData,
    njobs=-1,
    window_size = 250000,
    organism = 'human'
    ) -> pd.DataFrame:
    """
    Compute co-accessibility network using circe.
    
    Args:
        atac_data: Preprocessed AnnData object
        
    Returns:
        DataFrame containing the co-accessibility network with renamed 'coaccess' column
    """

    atac_data.var.index = list(map(lambda x:x.replace(":", "_").replace("-", "_"), atac_data.var.index))
    # Add region information using circe
    atac_data = ci.add_region_infos(atac_data)

    # Compute the co-accessibility network
    ci.compute_atac_network(
        atac_data,
        window_size=window_size,
        unit_distance=1000,
        distance_constraint=None,
        s=None,
        organism=organism,
        max_alpha_iteration=100,
        distance_parameter_convergence=1e-22,
        max_elements=200,
        n_samples=100,
        n_samples_maxtry=500,
        key="atac_network",
        seed=42,
        njobs=njobs,
        threads_per_worker=1,
        verbose=0,
        chromosomes_sizes=None
        )
    
    # Extract the network
    circe_network = ci.extract_atac_links(atac_data)
    
    # Rename score column to coaccess
    circe_network = circe_network.rename(columns={"score": "coaccess"})
    
    return circe_network


def filter_coaccess_peaks(coaccess_df, threshold:float=None, quantile_threshold=0.9):

    if threshold is None:
        threshold = coaccess_df["coaccess"].quantile(quantile_threshold)
    
    return coaccess_df[coaccess_df["coaccess"] >= threshold]
    






