# BuenColors Examples

This directory contains Jupyter notebooks demonstrating all features of the BuenColors package.

## Notebooks

### 1. Helper Functions Demo (`helpers_demo.ipynb`)

Demonstrates all helper functions including:
- `eject_legend()` - Move legends outside plots
- `rotate_discrete_xticks()` - Rotate axis labels
- `grab_legend()` - Extract legends to separate figures (with optional `remove` parameter)
- `get_density()` - Compute point density for scatter plots
- `shuffle()` - Randomize data order
- `number_to_color()` - Map values to colors

### 2. Palettes Demo (`palettes_demo.ipynb`)

Showcases the 117+ color palettes including:
- `list_palettes()` - Browse available palettes
- `display_palette()` - Visualize palettes
- `get_palette()` - Extract colors
- Wes Anderson collection
- Scientific visualization palettes
- Custom gradients

### 3. Single-Cell Analysis Demo (`single_cell_demo.ipynb`)

Demonstrates single-cell RNA-seq visualization functions:
- `clean_umap()` - Publication-quality UMAP plots
- Gene expression visualization
- Custom axis indicators
- Multi-panel figures for publications
- Integration with Scanpy workflows
- Using BuenColors palettes with single-cell data

**Note:** Requires `scanpy` and `anndata` to be installed.

## Running the Notebooks

```bash
# Install jupyter if needed
pip install jupyter

# Launch Jupyter
jupyter notebook

# Or use JupyterLab
pip install jupyterlab
jupyter lab
```

Then open the desired notebook and run all cells to see the examples.

## Requirements

All notebooks require:
- `buencolors` (this package)
- `numpy`
- `pandas`
- `matplotlib`

These are automatically installed with BuenColors.

The single-cell demo additionally requires:
- `scanpy`
- `anndata`

Install with: `pip install scanpy anndata`
