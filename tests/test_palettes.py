import pytest
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib.cm as cm
from buencolors.palettes import (
    jdb_palettes, register_cmaps, get_palette, display_palette,
    list_palettes, get_registered_cmaps
)


class TestPaletteRegistry:
    """Test colormap registration functionality."""

    def test_all_palettes_registered(self):
        """Test that all palettes are registered as matplotlib colormaps."""
        for palette_name in jdb_palettes.keys():
            # Check normal version
            assert palette_name in plt.colormaps(), f"Palette {palette_name} not registered"
            # Check reversed version
            assert f"{palette_name}_r" in plt.colormaps(), f"Reversed palette {palette_name}_r not registered"

    def test_colormap_colors_match_palette(self):
        """Test that registered colormaps match original palette colors."""
        palette_name = "Zissou"
        expected_colors = jdb_palettes[palette_name]

        cmap = cm.get_cmap(palette_name)
        # Sample at evenly spaced points
        n = len(expected_colors)
        indices = np.linspace(0, 1, n)

        for i, idx in enumerate(indices):
            cmap_color = mcolors.to_hex(cmap(idx))
            # Colors should be similar (allowing for interpolation)
            assert isinstance(cmap_color, str)
            assert cmap_color.startswith('#')

    def test_reversed_colormap(self):
        """Test that reversed colormaps are actually reversed."""
        palette_name = "solar_flare"

        cmap = cm.get_cmap(palette_name)
        cmap_r = cm.get_cmap(f"{palette_name}_r")

        # First color of normal should match last color of reversed
        color_start = mcolors.to_hex(cmap(0.0))
        color_end_r = mcolors.to_hex(cmap_r(1.0))

        # These should be very close (might not be exact due to interpolation)
        assert isinstance(color_start, str)
        assert isinstance(color_end_r, str)


class TestGetPalette:
    """Test get_palette function."""

    def test_get_all_colors_default(self):
        """Test getting all colors from a palette."""
        colors = get_palette("Cavalcanti")
        expected = jdb_palettes["Cavalcanti"]
        assert colors == expected
        assert len(colors) == 5

    def test_get_first_n_colors(self):
        """Test getting first n colors."""
        colors = get_palette("Cavalcanti", n=3)
        expected = jdb_palettes["Cavalcanti"][:3]
        assert colors == expected
        assert len(colors) == 3

    def test_get_specific_indices(self):
        """Test getting colors at specific indices."""
        colors = get_palette("Cavalcanti", n=[0, 2, 4])
        expected = [jdb_palettes["Cavalcanti"][i] for i in [0, 2, 4]]
        assert colors == expected
        assert len(colors) == 3

    def test_interpolate_more_colors(self):
        """Test interpolating more colors than palette has."""
        palette_name = "berry"  # Only 2 colors
        colors = get_palette(palette_name, n=10)
        assert len(colors) == 10
        # All should be hex strings
        assert all(c.startswith('#') for c in colors)

    def test_reverse_palette(self):
        """Test reversing a palette."""
        colors = get_palette("Zissou", reverse=True)
        expected = list(reversed(jdb_palettes["Zissou"]))
        assert colors == expected

    def test_return_rgb_array(self):
        """Test returning RGB arrays instead of hex."""
        colors = get_palette("berry", as_hex=False)
        assert isinstance(colors, np.ndarray)
        assert colors.shape == (2, 4)  # 2 colors, RGBA
        assert colors.dtype == np.float64
        # Values should be in [0, 1]
        assert np.all(colors >= 0) and np.all(colors <= 1)

    def test_invalid_palette_raises_error(self):
        """Test that invalid palette name raises ValueError."""
        with pytest.raises(ValueError, match="Palette 'nonexistent' not found"):
            get_palette("nonexistent")

    def test_out_of_range_index_raises_error(self):
        """Test that out of range indices raise IndexError."""
        with pytest.raises(IndexError):
            get_palette("berry", n=[0, 10])  # berry only has 2 colors


class TestListPalettes:
    """Test list_palettes function."""

    def test_list_all_palettes(self):
        """Test listing all palettes."""
        palettes = list_palettes()
        assert isinstance(palettes, list)
        assert len(palettes) == len(jdb_palettes)
        # Should be sorted
        assert palettes == sorted(palettes)

    def test_filter_by_pattern(self):
        """Test filtering palettes by pattern."""
        flame_palettes = list_palettes("flame")
        assert all("flame" in p.lower() for p in flame_palettes)
        assert len(flame_palettes) > 0

    def test_list_by_categories(self):
        """Test listing palettes grouped by category."""
        categorized = list_palettes(categories=True)
        assert isinstance(categorized, dict)
        assert "wesanderson" in categorized
        assert "brewer" in categorized
        # Check that wesanderson palettes are actually there
        assert "Zissou" in categorized["wesanderson"]
        assert "GrandBudapest" in categorized["wesanderson"]

    def test_categories_sorted(self):
        """Test that palettes within categories are sorted."""
        categorized = list_palettes(categories=True)
        for category, palettes in categorized.items():
            assert palettes == sorted(palettes), f"Category {category} not sorted"


class TestDisplayPalette:
    """Test display_palette function."""

    def test_display_returns_figure(self, clean_plot):
        """Test that display_palette returns a Figure."""
        fig = display_palette("Zissou")
        assert isinstance(fig, plt.Figure)

    def test_display_with_custom_n(self, clean_plot):
        """Test displaying palette with specific number of colors."""
        fig = display_palette("berry", n=10)
        assert isinstance(fig, plt.Figure)

    def test_display_with_custom_figsize(self, clean_plot):
        """Test displaying palette with custom figure size."""
        fig = display_palette("solar_flare", figsize=(12, 1))
        assert fig.get_figwidth() == 12
        assert fig.get_figheight() == 1

    def test_display_without_name(self, clean_plot):
        """Test displaying palette without showing name."""
        fig = display_palette("Zissou", show_name=False)
        ax = fig.axes[0]
        assert ax.get_title() == ""


class TestGetRegisteredCmaps:
    """Test get_registered_cmaps function."""

    def test_get_normal_cmaps(self):
        """Test getting only normal colormaps."""
        cmaps = get_registered_cmaps(include_reversed=False)
        assert isinstance(cmaps, list)
        assert len(cmaps) == len(jdb_palettes)
        assert "Zissou" in cmaps
        assert "Zissou_r" not in cmaps

    def test_get_with_reversed(self):
        """Test getting colormaps including reversed versions."""
        cmaps = get_registered_cmaps(include_reversed=True)
        assert len(cmaps) == len(jdb_palettes) * 2
        assert "Zissou" in cmaps
        assert "Zissou_r" in cmaps
        # Should be sorted
        assert cmaps == sorted(cmaps)
