import matplotlib
import matplotlib.style
from pathlib import Path
from .palettes import register_cmaps, get_palette, display_palette, list_palettes, get_registered_cmaps
from .helpers import eject_legend, rotate_discrete_xticks, grab_legend, get_density, shuffle, number_to_color, eject_colorbar, grab_colorbar
from .single_cell import clean_umap

# -----------------------------------------------------------------------------
# 1. Style Sheet Registration
# -----------------------------------------------------------------------------
# Determine the absolute path to the styles directory
# Using.parent ensures this works regardless of where the package is installed
PACKAGE_DIR = Path(__file__).parent
STYLE_DIR = PACKAGE_DIR / 'styles'



def _register_style_dir(style_dir: Path) -> None:
    """Make the bundled ``.mplstyle`` files visible to ``plt.style.use`` and
    ``plt.style.available`` without copying them anywhere.

    matplotlib >= 3.11 exposes ``USER_LIBRARY_PATHS`` on ``matplotlib.style``
    and deprecates ``matplotlib.style.core``; older versions only have it on
    ``matplotlib.style.core``. Try the public location first so no deprecation
    warning is triggered on new releases.
    """
    style_dir = str(style_dir)

    user_paths = getattr(matplotlib.style, "USER_LIBRARY_PATHS", None)  # matplotlib >= 3.11
    if user_paths is None:
        try:
            from matplotlib.style.core import USER_LIBRARY_PATHS as user_paths  # matplotlib < 3.11
        except ImportError:
            user_paths = None

    if user_paths is not None:
        if style_dir not in user_paths:
            user_paths.append(style_dir)
        # Re-index styles
        matplotlib.style.reload_library()
        return

    # Last-resort fallback using only public API: load each style sheet and
    # inject it into the style library directly.
    for path in sorted(Path(style_dir).glob("*.mplstyle")):
        matplotlib.style.library[path.stem] = matplotlib.rc_params_from_file(
            path, use_default_template=False
        )
    matplotlib.style.available[:] = sorted(
        name for name in matplotlib.style.library if not name.startswith("_")
    )


# This allows 'plt.style.use("pretty-plot")' to work without moving files.
_register_style_dir(STYLE_DIR)

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
    "eject_colorbar",
    "grab_colorbar",
    "clean_umap",
)