import matplotlib.lines as mlines

from .helpers import eject_legend, shuffle

try:
    import anndata as ad
    ANNDATA_AVAILABLE = True
except ImportError:
    ANNDATA_AVAILABLE = False

if ANNDATA_AVAILABLE:
    try:
        import scanpy as sc
        SCANPY_AVAILABLE = True
    except ImportError:
        SCANPY_AVAILABLE = False
else:
    SCANPY_AVAILABLE = False


if SCANPY_AVAILABLE:
    def clean_umap(adata, color, ax=None, axis_len=0.2, thickness=3.0, **kwargs):
        """Plot a clean UMAP with minimal decorations and custom axis indicators.

        Plots a Scanpy UMAP with no borders/ticks, but adds a small 'L' shaped
        axis indicator in the bottom left corner with arrowheads. The legend is
        ejected to the right side of the plot. Cells are automatically shuffled
        to avoid non-random ordering artifacts.

        Parameters
        ----------
        adata : AnnData
            AnnData object containing the UMAP coordinates
        color : str
            Column in adata.obs or gene name to color cells by
        ax : plt.Axes, optional
            Existing matplotlib axis. If None, creates a new axis.
        axis_len : float, optional
            Length of the custom axis arrows in relative axes coordinates (0-1). Default is 0.2.
        thickness : float, optional
            Line width of the custom axes. Default is 3.0.
        **kwargs
            Additional keyword arguments passed to sc.pl.umap

        Returns
        -------
        plt.Axes
            The matplotlib axes containing the UMAP plot

        Examples
        --------
        >>> import scanpy as sc
        >>> import buencolors
        >>>
        >>> # Load example dataset
        >>> adata = sc.datasets.pbmc3k_processed()
        >>>
        >>> # Plot clean UMAP colored by cell type
        >>> buencolors.clean_umap(adata, color='louvain')
        >>> plt.show()
        >>>
        >>> # Customize axis length and thickness
        >>> buencolors.clean_umap(adata, color='louvain', axis_len=0.3, thickness=4.0)
        >>> plt.show()
        >>>
        >>> # Use with existing axis and pass additional arguments
        >>> fig, ax = plt.subplots(figsize=(8, 6))
        >>> buencolors.clean_umap(adata, color='CST3', ax=ax, cmap='viridis')
        >>> plt.tight_layout()
        >>> plt.show()
        """

        # 1. Plot base UMAP with frame disabled
        # We will first shuffle the data to avoid non-random cell ordering issues
        adata = shuffle(adata.copy())
        ax = sc.pl.umap(adata, color=color, show=False, frameon=False, ax=ax, **kwargs)

        # 2. Ensure all standard decorations are gone (in case frameon misses something)
        ax.set_xlabel("")
        ax.set_ylabel("")
        ax.set_xticks([])
        ax.set_yticks([])
        # Remove spines
        for spine in ax.spines.values():
            spine.set_visible(False)
        # Set aspect to equal
        ax.set_aspect('equal')
        # Set colors of background to transparent
        ax.set_facecolor('none')
        ax.figure.patch.set_alpha(0)

        # 3. Eject the legend to the right side of the plot
        eject_legend(ax)

        # 4. Add Custom Axes (The "L" shape) with arrowheads
        # We use ax.transAxes so coordinates are relative (0,0 = bottom-left)
        # Use small offset to prevent arrows from being cut off
        offset_x = 0.03
        offset_y = 0.03

        # Draw the L-shape using Line2D with markers for arrowheads
        # Vertical line (UMAP 2) with arrowhead pointing up at the end only
        line_v = mlines.Line2D([offset_x, offset_x], [offset_y, axis_len + offset_y],
                               transform=ax.transAxes, color='black', linewidth=thickness,
                               marker='^', markersize=10, markevery=[1],
                               markerfacecolor='black', markeredgecolor='black',
                               solid_capstyle='projecting', solid_joinstyle='miter')
        # Horizontal line (UMAP 1) with arrowhead pointing right at the end only
        line_h = mlines.Line2D([offset_x, axis_len + offset_x], [offset_y, offset_y],
                               transform=ax.transAxes, color='black', linewidth=thickness,
                               marker='>', markersize=10, markevery=[1],
                               markerfacecolor='black', markeredgecolor='black',
                               solid_capstyle='projecting', solid_joinstyle='miter')

        ax.add_line(line_v)
        ax.add_line(line_h)


        # 5. Add Labels
        # Position labels to emulate xlabel/ylabel: centered along each axis, anchored at origin
        # Horizontal label (UMAP 1) - centered below the horizontal axis
        ax.text(axis_len/2 + offset_x, offset_y - 0.03, "UMAP 1", transform=ax.transAxes,
                ha='center', va='top', fontsize=10)
        # Vertical label (UMAP 2) - centered along vertical axis, rotated 90 degrees
        ax.text(offset_x - 0.03, axis_len/2 + offset_y, "UMAP 2", transform=ax.transAxes,
                ha='right', va='center', fontsize=10, rotation=90)

        # 6. Add buffer region to prevent arrows and legend from being cut off
        # Use tight_layout with rect parameter to preserve space while keeping aspect ratio
        # Wrap in try/except because tight_layout conflicts with colorbars (gene expression plots)
        try:
            ax.figure.tight_layout(rect=(0.05, 0.05, 0.95, 0.95))
        except RuntimeError:
            # Colorbar present, skip tight_layout
            pass
        # Adjust subplot parameters to ensure nothing is clipped
        ax.figure.subplots_adjust(left=0.1, bottom=0.1, right=0.85)

        return ax