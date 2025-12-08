import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import gaussian_kde
import pandas as pd
from typing import TypeVar, overload, Optional, Union
from numpy.typing import ArrayLike
import matplotlib.colors as mcolors
import matplotlib.cm as cm

try:
    import anndata  # noqa: F401
    ANNDATA_AVAILABLE = True
except ImportError:
    ANNDATA_AVAILABLE = False


def eject_legend(ax: Optional[plt.Axes] = None) -> None:
    """Eject the legend from the plot area to the right side of the axes.

    This is useful for creating publication-quality plots where the legend
    should not overlap with the data.

    Parameters
    ----------
    ax : plt.Axes, optional
        The axes to modify. If None, uses the current axes.

    Returns
    -------
    None

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import buencolors
    >>>
    >>> # Create a plot with multiple lines
    >>> x = np.linspace(0, 10, 100)
    >>> plt.plot(x, np.sin(x), label='sin(x)')
    >>> plt.plot(x, np.cos(x), label='cos(x)')
    >>> plt.plot(x, np.sin(x) * np.cos(x), label='sin(x)cos(x)')
    >>>
    >>> # Eject legend to the right
    >>> buencolors.eject_legend()
    >>> plt.show()
    >>>
    >>> # Or use with a specific axes
    >>> fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))
    >>> ax1.plot(x, np.sin(x), label='sin(x)')
    >>> ax1.plot(x, np.cos(x), label='cos(x)')
    >>> ax2.plot(x, x**2, label='x²')
    >>> ax2.plot(x, x**3, label='x³')
    >>>
    >>> buencolors.eject_legend(ax1)
    >>> buencolors.eject_legend(ax2)
    >>> plt.show()
    """
    if ax is None:
        ax = plt.gca()
    ax.legend(bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)


def rotate_discrete_xticks(ax: Optional[plt.Axes] = None, rotation: float = 45) -> None:
    """Rotate the x-tick labels for discrete axes with proper alignment.

    This will rotate the x-tick labels and ensure that the ending of the labels
    are directly under the ticks for better readability.

    Parameters
    ----------
    ax : plt.Axes, optional
        The axes to modify. If None, uses the current axes.
    rotation : float, optional
        The angle of rotation for the x-tick labels. Default is 45 degrees.

    Returns
    -------
    None

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import buencolors
    >>>
    >>> # Bar plot with long category names
    >>> categories = ['First Category', 'Second Category', 'Third Category',
    ...               'Fourth Category', 'Fifth Category']
    >>> values = [23, 45, 56, 78, 32]
    >>>
    >>> plt.bar(categories, values)
    >>> buencolors.rotate_discrete_xticks
    >>> plt.show()
    >>>
    >>> # Custom rotation angle
    >>> plt.bar(categories, values)
    >>> buencolors.rotate_discrete_xticks(rotation=60)
    >>> plt.show()
    >>>
    >>> # Use with specific axes
    >>> fig, ax = plt.subplots()
    >>> ax.bar(categories, values)
    >>> buencolors.rotate_discrete_xticks(ax, rotation=30)
    >>> plt.show()
    """
    if ax is None:
        ax = plt.gca()
    for label in ax.get_xticklabels():
        label.set_rotation(rotation)
        label.set_ha('right')


def grab_legend(ax: Optional[plt.Axes] = None, remove: bool = True) -> plt.Figure:
    """Grab the legend from the axes and return it in a new figure for external saving or modification.

    This is useful for creating separate legend files for publications or presentations,
    or for combining legends from multiple plots.

    Parameters
    ----------
    ax : plt.Axes, optional
        The axes to grab the legend from. If None, uses the current axes.
    remove : bool, optional
        If True (default), remove the legend from the original axes after extraction.
        If False, keep the legend on the original axes.

    Returns
    -------
    plt.Figure
        A new figure containing only the legend.

    Examples
    --------
    >>> import matplotlib.pyplot as plt
    >>> import numpy as np
    >>> import buencolors
    >>>
    >>> # Create a plot and extract legend (removing it from original)
    >>> x = np.linspace(0, 10, 100)
    >>> plt.plot(x, np.sin(x), label='sin(x)', color='blue')
    >>> plt.plot(x, np.cos(x), label='cos(x)', color='red')
    >>> plt.plot(x, np.tan(x), label='tan(x)', color='green')
    >>> plt.legend()
    >>> plt.ylim(-2, 2)
    >>>
    >>> # Extract and save legend separately (legend removed from plot)
    >>> legend_fig = buencolors.grab_legend()
    >>> legend_fig.savefig('legend.pdf', bbox_inches='tight')
    >>> plt.savefig('plot.pdf')  # Plot saved without legend
    >>> plt.show()
    >>>
    >>> # Extract legend while keeping it on the original axes
    >>> fig, ax = plt.subplots()
    >>> ax.scatter([1, 2, 3], [1, 4, 9], label='Data A', color='blue')
    >>> ax.scatter([1, 2, 3], [2, 3, 5], label='Data B', color='red')
    >>> ax.legend()
    >>>
    >>> # Keep legend on original plot
    >>> legend_fig = buencolors.grab_legend(ax, remove=False)
    >>> legend_fig.savefig('my_legend.png', dpi=300, bbox_inches='tight')
    >>> plt.show()  # Original plot still has legend
    >>>
    >>> # Remove legend from original (default behavior)
    >>> fig, ax = plt.subplots()
    >>> ax.plot(x, np.exp(x), label='exp(x)')
    >>> ax.legend()
    >>> legend_fig = buencolors.grab_legend(ax, remove=True)
    >>> plt.show()  # Original plot has no legend
    """
    if ax is None:
        ax = plt.gca()

    # Get the legend from the axes
    legend = ax.get_legend()

    if legend is None:
        raise ValueError("No legend found on the provided axes")

    # Get handles and labels before removing the legend
    handles, labels = ax.get_legend_handles_labels()
    frameon = legend.get_frame_on()

    # Draw the canvas to get accurate legend dimensions
    ax.figure.canvas.draw()

    # Get legend bounding box in inches
    bbox = legend.get_window_extent().transformed(ax.figure.dpi_scale_trans.inverted())

    # Remove the legend from the original axes if requested
    if remove:
        legend.remove()

    # Create a new figure for the legend with appropriate size
    fig, ax_legend = plt.subplots(figsize=(bbox.width, bbox.height))

    # Hide the axis
    ax_legend.axis('off')

    # Ensure no title on the legend figure
    ax_legend.set_title('')
    fig.suptitle('')

    # Add the legend to the axis
    ax_legend.legend(handles=handles,
                     labels=labels,
                     loc='center',
                     frameon=frameon)

    return fig


def get_density(x: ArrayLike, y: ArrayLike, n: int = 200) -> np.ndarray:
    """Compute the density of points in a grid.

    Based on the R implementation. This is useful for coloring scatter plot points
    by their local density to visualize data concentration.

    Credit to Kamil Slowikowski
    See post: http://slowkow.com/notes/ggplot2-color-by-density/

    Parameters
    ----------
    x : array-like
        X coordinates
    y : array-like
        Y coordinates
    n : int, optional
        Number of bins to divide grid. Default is 200. Higher values give more detail but are slower.

    Returns
    -------
    np.ndarray
        A vector of densities for plotting

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import buencolors
    >>>
    >>> # Create sample data with varying density
    >>> x = np.concatenate([
    ...     np.random.normal(0, 0.1, 10000),
    ...     np.random.normal(0, 0.1, 1000)
    ... ])
    >>> y = np.concatenate([
    ...     np.random.normal(0, 0.1, 10000),
    ...     np.random.normal(0.1, 0.2, 1000)
    ... ])
    >>>
    >>> # Compute density for each point
    >>> density = buencolors.get_density(x, y)
    >>>
    >>> # Plot with color mapped to density
    >>> plt.scatter(x, y, c=density, cmap='viridis', s=1)
    >>> plt.colorbar(label='Point Density')
    >>> plt.show()
    >>>
    >>> # Use with a BuenColors palette
    >>> density = buencolors.get_density(x, y, n=300)  # Higher resolution
    >>> plt.scatter(x, y, c=density, cmap='solar_flare', s=1, alpha=0.5)
    >>> plt.colorbar()
    >>> plt.show()
    """
    # Convert to numpy arrays
    x = np.asarray(x)
    y = np.asarray(y)

    # Create 2D kernel density estimation
    kde = gaussian_kde(np.vstack([x, y]))

    # Create grid
    x_grid = np.linspace(x.min(), x.max(), n)
    y_grid = np.linspace(y.min(), y.max(), n)

    # Evaluate KDE on grid
    xx, yy = np.meshgrid(x_grid, y_grid)
    positions = np.vstack([xx.ravel(), yy.ravel()])
    z = np.reshape(kde(positions), xx.shape)

    # Find which bin each point falls into
    # digitize returns 1-based bin indices, so subtract 1 for 0-based indexing
    ix = np.digitize(x, x_grid) - 1
    iy = np.digitize(y, y_grid) - 1

    # Clip indices to valid range [0, n-1]
    ix = np.clip(ix, 0, n - 1)
    iy = np.clip(iy, 0, n - 1)

    # Return density values at each point
    return z[iy, ix]


@overload
def shuffle(x: pd.DataFrame) -> pd.DataFrame: ...

@overload
def shuffle(x: pd.Series) -> pd.Series: ...

@overload
def shuffle(x: np.ndarray) -> np.ndarray: ...

@overload
def shuffle(x: ArrayLike) -> list: ...

@overload
def shuffle(x: list) -> list: ...

@overload
def shuffle(x: tuple) -> tuple: ...

if ANNDATA_AVAILABLE:
    @overload
    def shuffle(x: anndata.AnnData) -> anndata.AnnData: ...

def shuffle(x):
    """Shuffle the order of rows/elements to make plots independent of point ordering.

    This function accepts various data types and returns a shuffled version while
    preserving the original type when possible. Useful for preventing visual artifacts
    caused by data ordering in scatter plots.

    Parameters
    ----------
    x : array-like, DataFrame, or Series
        Data to shuffle. Can be a list, tuple, numpy array, pandas DataFrame, or pandas Series.

    Returns
    -------
    same type as input
        Shuffled version of the input data. For DataFrames and Series, the index is reset.
        For numpy arrays, a shuffled copy is returned. For lists, a shuffled list is returned.

    Examples
    --------
    >>> import numpy as np
    >>> import pandas as pd
    >>> import matplotlib.pyplot as plt
    >>> import buencolors
    >>>
    >>> # Shuffle a DataFrame for plotting
    >>> df = pd.DataFrame({
    ...     'x': np.random.randn(1000),
    ...     'y': np.random.randn(1000),
    ...     'category': np.random.choice(['A', 'B', 'C'], 1000)
    ... })
    >>>
    >>> # Without shuffling, later categories may cover earlier ones
    >>> df_shuffled = buencolors.shuffle(df)
    >>> colors = {'A': 'red', 'B': 'blue', 'C': 'green'}
    >>> df_shuffled['color'] = df_shuffled['category'].map(colors)
    >>> plt.scatter(df_shuffled['x'], df_shuffled['y'], c=df_shuffled['color'], alpha=0.5)
    >>> plt.show()
    >>>
    >>> # Shuffle a numpy array
    >>> arr = np.random.randn(100, 2)
    >>> arr_shuffled = buencolors.shuffle(arr)
    >>>
    >>> # Shuffle a list
    >>> lst = [1, 2, 3, 4, 5]
    >>> lst_shuffled = buencolors.shuffle(lst)

    Notes
    -----
    This also supports scanpy AnnData objects by shuffling the obs index when the anndata package is installed.
    """
    if isinstance(x, pd.DataFrame):
        # Shuffle DataFrame rows and reset index
        return x.sample(frac=1).reset_index(drop=False)
    elif isinstance(x, pd.Series):
        # Shuffle Series and reset index
        return x.sample(frac=1).reset_index(drop=False)
    elif isinstance(x, np.ndarray):
        # Shuffle numpy array rows (axis 0)
        shuffled = x.copy()
        np.random.shuffle(shuffled)
        return shuffled
    elif isinstance(x, (list, tuple)):
        # Shuffle list or tuple
        x_list = list(x)
        np.random.shuffle(x_list)
        # Return same type as input
        return type(x)(x_list) if isinstance(x, tuple) else x_list
    elif ANNDATA_AVAILABLE and 'anndata' in str(type(x)).lower():
        # Shuffle AnnData obs index
        idx = np.random.permutation(x.n_obs)
        return x[idx, :]
    else:
        # Handle other iterables by converting to list
        x_list = list(x)
        np.random.shuffle(x_list)
        return x_list


def number_to_color(
    values: ArrayLike,
    palette: Union[str, mcolors.Colormap],
    value_range: Optional[tuple[float, float]] = None,
    n_bins: Optional[int] = None,
    return_rgb: bool = False
) -> Union[list[str], np.ndarray]:
    """Map numeric values to colors from a palette.

    This is the pythonic equivalent of R's numberToColorVec function, with support
    for both continuous and discrete (binned) color mapping.

    Parameters
    ----------
    values : array-like
        Numeric values to map to colors
    palette : str or Colormap
        Name of a registered colormap or a Colormap object. Can be any matplotlib
        colormap or a BuenColors palette name.
    value_range : tuple[float, float], optional
        (min, max) values to clip and normalize. If None, uses the min/max of values.
        Values outside this range will be clipped to the range bounds.
    n_bins : int, optional
        Number of discrete color bins to use. If None (default), uses continuous mapping.
        If specified, values will be binned into n_bins discrete colors, similar to
        the R implementation (which uses 100 bins). Common values: 100, 256.
    return_rgb : bool, optional
        If True, return RGB arrays instead of hex strings. Default is False.

    Returns
    -------
    list[str] or np.ndarray
        If return_rgb is False: list of hex color strings (e.g., '#FF5733')
        If return_rgb is True: numpy array of shape (n, 4) with RGBA values in [0, 1]

    Examples
    --------
    >>> import numpy as np
    >>> import matplotlib.pyplot as plt
    >>> import buencolors
    >>>
    >>> # Basic usage: map values to colors
    >>> values = np.random.randn(1000) ** 2
    >>> colors = buencolors.number_to_color(values, "viridis")
    >>> plt.scatter(np.arange(len(values)), values, c=colors)
    >>> plt.show()
    >>>
    >>> # Use a BuenColors palette
    >>> colors = buencolors.number_to_color(values, "solar_flare")
    >>> plt.scatter(np.arange(len(values)), values, c=colors)
    >>> plt.show()
    >>>
    >>> # Discrete binned colors (like R implementation)
    >>> colors = buencolors.number_to_color(values, "brewer_spectra", n_bins=100)
    >>> plt.scatter(np.arange(len(values)), values, c=colors, s=20)
    >>> plt.show()
    >>>
    >>> # With custom value range (capping)
    >>> colors = buencolors.number_to_color(values, "plasma", value_range=(0, 2))
    >>> plt.scatter(np.arange(len(values)), values, c=colors)
    >>> plt.show()
    >>>
    >>> # Get RGB arrays instead of hex
    >>> colors_rgb = buencolors.number_to_color(values[:10], "magma", return_rgb=True)
    >>> print(colors_rgb.shape)  # (10, 4) RGBA values
    (10, 4)
    >>>
    >>> # Combined with shuffle for better visualization
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     'x': np.random.randn(500),
    ...     'y': np.random.randn(500),
    ...     'value': np.random.randn(500) ** 2
    ... })
    >>> df['color'] = buencolors.number_to_color(df['value'], "flame_flame", n_bins=50)
    >>> df = buencolors.shuffle(df)
    >>> plt.scatter(df['x'], df['y'], c=df['color'], s=30, alpha=0.6)
    >>> plt.show()

    Notes
    -----
    The default behavior (n_bins=None) uses continuous color mapping for smoother
    gradients. Set n_bins=100 to match the R implementation's discrete binning.

    For more advanced discrete colormapping, consider using matplotlib's BoundaryNorm.
    """
    values = np.asarray(values, dtype=float)

    # Get the colormap
    if isinstance(palette, str):
        cmap = cm.get_cmap(palette)
    else:
        cmap = palette

    # Handle range clipping and normalization
    if value_range is not None:
        vmin, vmax = value_range
        # Clip values to the specified range
        clipped_values = np.clip(values, vmin, vmax)
    else:
        clipped_values = values.copy()
        vmin, vmax = np.nanmin(clipped_values), np.nanmax(clipped_values)

    # Handle edge case where all values are the same
    if vmin == vmax:
        # Return the middle color of the colormap for all values
        if return_rgb:
            return np.tile(cmap(0.5), (len(values), 1))
        else:
            return [mcolors.to_hex(cmap(0.5))] * len(values)

    if n_bins is not None:
        # Discrete binning mode (like R implementation)
        # Create bins
        breaks = np.linspace(vmin, vmax, n_bins)

        # Sample colors from the colormap at evenly spaced intervals
        color_indices = np.linspace(0, 1, n_bins + 1)
        colors_array = cmap(color_indices)

        # Bin each value into the appropriate color
        # digitize returns bin indices (1-based), we want 0-based for color indexing
        bin_indices = np.digitize(clipped_values, breaks)
        # Clip to valid range [0, n_bins]
        bin_indices = np.clip(bin_indices, 0, n_bins)

        if return_rgb:
            # Return RGBA arrays for each binned value
            return colors_array[bin_indices]
        else:
            # Return hex strings for each binned value
            return [mcolors.to_hex(colors_array[idx]) if not np.isnan(clipped_values[i]) else '#000000'
                    for i, idx in enumerate(bin_indices)]
    else:
        # Continuous mapping mode (default)
        # Normalize values to [0, 1]
        norm = mcolors.Normalize(vmin=vmin, vmax=vmax)
        normalized = norm(clipped_values)

        # Map to colors
        if return_rgb:
            # Return RGBA arrays
            return cmap(normalized)
        else:
            # Return hex strings
            return [mcolors.to_hex(cmap(v)) if not np.isnan(v) else '#000000'
                    for v in normalized]