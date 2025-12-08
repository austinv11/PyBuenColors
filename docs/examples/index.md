# Examples

Interactive Jupyter notebooks demonstrating the features and capabilities of PyBuenColors.

## Available Notebooks

### [Helper Functions](helpers_demo.ipynb)

Comprehensive demonstration of utility functions for creating publication-ready plots:

- **`eject_legend()`** - Move legends outside the plot area to avoid obscuring data
- **`rotate_discrete_xticks()`** - Rotate x-axis labels for better readability
- **`grab_legend()`** - Extract legends to separate figures for independent saving
- **`get_density()`** - Compute point density for colored scatter plots
- **`shuffle()`** - Randomize data order to prevent plotting artifacts
- **`number_to_color()`** - Map numerical values to colors from any palette

### [Color Palettes](palettes_demo.ipynb)

Explore the 117+ scientific color palettes available in PyBuenColors:

- **`list_palettes()`** - Browse and filter available palettes
- **`display_palette()`** - Visualize individual palettes
- **`get_palette()`** - Extract colors for use in plots
- Wes Anderson-inspired palettes
- Scientific visualization color schemes
- Sequential and diverging gradients

### [Single-Cell Analysis](single_cell_demo.ipynb)

Specialized functions for single-cell RNA-seq visualization with Scanpy integration:

- **`clean_umap()`** - Create publication-quality UMAP plots with minimal decorations
- Gene expression visualization with custom colormaps
- L-shaped axis indicators for dimensional reduction plots
- Multi-panel figures for publications
- Integration with Scanpy workflows

!!! note "Optional Dependencies"
    The single-cell analysis notebook requires `scanpy` and `anndata`:
    ```bash
    pip install buencolors scanpy anndata
    ```

## Running Locally

To run these notebooks on your own machine:

```bash
# Install Jupyter
pip install jupyter

# Clone the repository
git clone https://github.com/austinv11/PyBuenColors.git
cd PyBuenColors/docs/examples

# Launch Jupyter
jupyter notebook
```

All notebooks are designed to run independently and include detailed explanations of each feature.