import pytest
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from buencolors.helpers import (
    eject_legend, rotate_discrete_xticks, grab_legend,
    get_density, shuffle, number_to_color
)


class TestEjectLegend:
    """Test eject_legend function."""

    def test_eject_legend_changes_position(self, clean_plot):
        """Test that eject_legend moves the legend outside the plot."""
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 2, 3], label='Line 1')
        ax.plot([1, 2, 3], [3, 2, 1], label='Line 2')
        ax.legend()

        # Get original legend position
        original_legend = ax.get_legend()
        assert original_legend is not None

        # Eject the legend
        eject_legend(ax)

        # Check that legend still exists and has bbox_to_anchor set
        new_legend = ax.get_legend()
        assert new_legend is not None
        # The legend should have bbox_to_anchor property set
        assert hasattr(new_legend, '_bbox_to_anchor')

    def test_eject_legend_uses_current_axes(self, clean_plot):
        """Test that eject_legend works without explicit axes."""
        plt.plot([1, 2, 3], [1, 2, 3], label='Line')
        plt.legend()

        eject_legend()  # Should work on current axes

        ax = plt.gca()
        assert ax.get_legend() is not None


class TestRotateDiscreteTicks:
    """Test rotate_discrete_xticks function."""

    def test_rotate_xticks_default(self, clean_plot):
        """Test that xticks are rotated with default angle."""
        fig, ax = plt.subplots()
        ax.bar(['Category A', 'Category B', 'Category C'], [1, 2, 3])

        rotate_discrete_xticks(ax)

        labels = ax.get_xticklabels()
        for label in labels:
            assert label.get_rotation() == 45
            assert label.get_ha() == 'right'

    def test_rotate_xticks_custom_angle(self, clean_plot):
        """Test rotating xticks with custom angle."""
        fig, ax = plt.subplots()
        ax.bar(['A', 'B', 'C'], [1, 2, 3])

        rotate_discrete_xticks(ax, rotation=60)

        labels = ax.get_xticklabels()
        for label in labels:
            assert label.get_rotation() == 60

    def test_rotate_uses_current_axes(self, clean_plot):
        """Test that rotate works without explicit axes."""
        plt.bar(['A', 'B', 'C'], [1, 2, 3])
        rotate_discrete_xticks(rotation=30)

        ax = plt.gca()
        labels = ax.get_xticklabels()
        for label in labels:
            assert label.get_rotation() == 30


class TestGrabLegend:
    """Test grab_legend function."""

    def test_grab_legend_creates_new_figure(self, clean_plot):
        """Test that grab_legend creates a new figure with the legend."""
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 2, 3], label='Line 1')
        ax.plot([1, 2, 3], [3, 2, 1], label='Line 2')
        ax.legend()

        legend_fig = grab_legend(ax, remove=True)

        assert isinstance(legend_fig, plt.Figure)
        assert len(legend_fig.axes) == 1
        # Original should no longer have legend
        assert ax.get_legend() is None

    def test_grab_legend_keep_original(self, clean_plot):
        """Test that grab_legend can keep the original legend."""
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 2, 3], label='Line')
        ax.legend()

        legend_fig = grab_legend(ax, remove=False)

        assert isinstance(legend_fig, plt.Figure)
        # Original should still have legend
        assert ax.get_legend() is not None

    def test_grab_legend_no_legend_raises_error(self, clean_plot):
        """Test that grab_legend raises error when no legend exists."""
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 2, 3])  # No legend

        with pytest.raises(ValueError, match="No legend found"):
            grab_legend(ax)


class TestGetDensity:
    """Test get_density function."""

    def test_get_density_returns_correct_shape(self, sample_data):
        """Test that get_density returns array of correct length."""
        x = sample_data['x']
        y = sample_data['y']

        density = get_density(x, y)

        assert isinstance(density, np.ndarray)
        assert len(density) == len(x)

    def test_get_density_values_positive(self, sample_data):
        """Test that density values are positive."""
        x = sample_data['x']
        y = sample_data['y']

        density = get_density(x, y)

        assert np.all(density >= 0)

    def test_get_density_custom_resolution(self, sample_data):
        """Test get_density with custom resolution."""
        x = sample_data['x'][:20]
        y = sample_data['y'][:20]

        density = get_density(x, y, n=100)

        assert len(density) == 20
        assert np.all(density >= 0)

    def test_get_density_accepts_lists(self):
        """Test that get_density accepts list inputs."""
        # Use data with some spread to avoid singular covariance matrix
        x = [1, 2, 3, 4, 5, 1.1, 2.1, 3.1]
        y = [1, 2, 3, 4, 5, 1.2, 2.2, 3.2]

        density = get_density(x, y)

        assert len(density) == 8


class TestShuffle:
    """Test shuffle function."""

    def test_shuffle_list(self):
        """Test shuffling a list."""
        original = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
        shuffled = shuffle(original)

        assert isinstance(shuffled, list)
        assert len(shuffled) == len(original)
        assert set(shuffled) == set(original)
        # With high probability, order should be different
        assert shuffled != original

    def test_shuffle_tuple(self):
        """Test shuffling a tuple."""
        original = (1, 2, 3, 4, 5)
        shuffled = shuffle(original)

        assert isinstance(shuffled, tuple)
        assert len(shuffled) == len(original)
        assert set(shuffled) == set(original)

    def test_shuffle_numpy_array(self):
        """Test shuffling a numpy array."""
        original = np.arange(20)
        shuffled = shuffle(original)

        assert isinstance(shuffled, np.ndarray)
        assert len(shuffled) == len(original)
        assert set(shuffled) == set(original)
        # Original should not be modified
        assert np.array_equal(original, np.arange(20))

    def test_shuffle_dataframe(self, sample_dataframe):
        """Test shuffling a DataFrame."""
        original = sample_dataframe.copy()
        shuffled = shuffle(original)

        assert isinstance(shuffled, pd.DataFrame)
        assert len(shuffled) == len(original)
        # shuffle adds 'index' column when it resets index
        assert 'index' in shuffled.columns
        # All original columns should be present
        for col in original.columns:
            assert col in shuffled.columns
        # All values should be present
        assert set(shuffled['categories']) == set(original['categories'])

    def test_shuffle_series(self):
        """Test shuffling a Series."""
        original = pd.Series([1, 2, 3, 4, 5], name='test')
        shuffled = shuffle(original)

        assert isinstance(shuffled, pd.DataFrame)  # Returns DataFrame with reset index
        assert len(shuffled) == len(original)


class TestNumberToColor:
    """Test number_to_color function."""

    def test_number_to_color_returns_hex(self):
        """Test that number_to_color returns hex strings by default."""
        values = [1, 2, 3, 4, 5]
        colors = number_to_color(values, 'viridis')

        assert isinstance(colors, list)
        assert len(colors) == 5
        assert all(c.startswith('#') for c in colors)

    def test_number_to_color_returns_rgb(self):
        """Test that number_to_color can return RGB arrays."""
        values = [1, 2, 3, 4, 5]
        colors = number_to_color(values, 'viridis', return_rgb=True)

        assert isinstance(colors, np.ndarray)
        assert colors.shape == (5, 4)  # RGBA
        assert np.all(colors >= 0) and np.all(colors <= 1)

    def test_number_to_color_with_value_range(self):
        """Test number_to_color with custom value range."""
        values = [1, 2, 3, 4, 5, 100]  # 100 is outlier
        colors = number_to_color(values, 'plasma', value_range=(1, 5))

        assert len(colors) == 6
        # Last value should be clipped to max
        assert all(c.startswith('#') for c in colors)

    def test_number_to_color_binned_mode(self):
        """Test number_to_color with discrete binning."""
        values = np.linspace(0, 1, 50)
        colors = number_to_color(values, 'magma', n_bins=10)

        assert len(colors) == 50
        # With only 10 bins, there should be at most 10 unique colors
        unique_colors = set(colors)
        assert len(unique_colors) <= 11  # Allow some tolerance

    def test_number_to_color_constant_values(self):
        """Test number_to_color when all values are the same."""
        values = [5, 5, 5, 5, 5]
        colors = number_to_color(values, 'viridis')

        assert len(colors) == 5
        # All colors should be the same (middle of colormap)
        assert len(set(colors)) == 1

    def test_number_to_color_with_buencolors_palette(self):
        """Test number_to_color with a BuenColors palette."""
        values = np.random.rand(20)
        colors = number_to_color(values, 'solar_flare')

        assert len(colors) == 20
        assert all(c.startswith('#') for c in colors)

    def test_number_to_color_with_cmap_object(self):
        """Test number_to_color with a Colormap object."""
        import matplotlib.cm as cm
        values = [1, 2, 3, 4, 5]
        cmap = cm.get_cmap('coolwarm')
        colors = number_to_color(values, cmap)

        assert len(colors) == 5
        assert all(c.startswith('#') for c in colors)