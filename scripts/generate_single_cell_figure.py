#!/usr/bin/env python3
"""
Generate single-cell analysis example figure for the README.

This script creates a demonstration of the clean_umap function using the
PBMC3k dataset.

Usage:
    python scripts/generate_single_cell_figure.py

Output:
    - figures/single_cell_clean_umap.png - Clean UMAP visualization example
"""

import matplotlib.pyplot as plt
import warnings
from pathlib import Path

warnings.filterwarnings('ignore')

# Import buencolors
import buencolors as bc

# Check for scanpy/anndata availability
try:
    import scanpy as sc
    import anndata as ad
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False
    print("Error: scanpy and anndata are required to generate single-cell figures.")
    print("Install with: pip install scanpy anndata")
    exit(1)


def create_output_dir():
    """Create figures directory if it doesn't exist."""
    figures_dir = Path("figures")
    figures_dir.mkdir(exist_ok=True)
    return figures_dir


def generate_clean_umap_example(output_dir):
    """Generate clean_umap() example showing cell type clustering."""
    print("Generating clean UMAP example...")

    # Load the preprocessed PBMC3k dataset
    adata = sc.datasets.pbmc3k_processed()

    # Create figure with clean UMAP
    with plt.style.context('pretty-plot'):
        fig = plt.figure(figsize=(10, 8))
        ax = bc.clean_umap(adata, color='louvain', palette='lawhoops')

        # Add a title
        fig.suptitle('PBMC3k Single-Cell UMAP - Cell Type Clusters',
                     fontsize=14, fontweight='bold', y=0.98)

        output_path = output_dir / 'single_cell_clean_umap.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"✓ Generated: {output_path}")
        plt.close()


def main():
    """Main function to generate all single-cell figures."""
    print("=" * 60)
    print("BuenColors Single-Cell Figure Generator")
    print("=" * 60)

    if not SCANPY_AVAILABLE:
        return

    # Create output directory
    output_dir = create_output_dir()

    # Generate examples
    generate_clean_umap_example(output_dir)

    print("\n" + "=" * 60)
    print("All single-cell figures generated successfully!")
    print("=" * 60)


if __name__ == "__main__":
    main()

