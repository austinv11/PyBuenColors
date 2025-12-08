import matplotlib.style
import matplotlib as mpl
from pathlib import Path
from .palettes import register_cmaps, get_palette, display_palette, list_palettes, get_registered_cmaps
from .helpers import eject_legend, rotate_discrete_xticks, grab_legend, get_density, shuffle, number_to_color
from .single_cell import clean_umap

# -----------------------------------------------------------------------------
# 1. Style Sheet Registration
# -----------------------------------------------------------------------------
# Determine the absolute path to the styles directory
# Using.parent ensures this works regardless of where the package is installed
PACKAGE_DIR = Path(__file__).parent
STYLE_DIR = PACKAGE_DIR / 'styles'

# Dynamic manipulation of USER_LIBRARY_PATHS
# This allows 'plt.style.use("pretty-plot")' to work without moving files.
if str(STYLE_DIR) not in matplotlib.style.core.USER_LIBRARY_PATHS:
    matplotlib.style.core.USER_LIBRARY_PATHS.append(str(STYLE_DIR))

# Re-index styles
matplotlib.style.reload_library()

# -----------------------------------------------------------------------------
# 2. Colormap Registration
# -----------------------------------------------------------------------------
# Execute the registration function from palettes.py
register_cmaps()

# -----------------------------------------------------------------------------
# 3. Metadata
# -----------------------------------------------------------------------------
try:
    from ._version import __version__
except ImportError:
    # Fallback version if _version.py doesn't exist yet
    __version__ = "unknown"


__all__ = (
    "eject_legend",
    "rotate_discrete_xticks",
    "grab_legend",
    "get_density",
    "shuffle",
    "number_to_color",
    "get_palette",
    "display_palette",
    "list_palettes",
    "get_registered_cmaps",
    "clean_umap",
)