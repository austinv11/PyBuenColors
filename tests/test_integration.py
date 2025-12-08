import pytest
import matplotlib.pyplot as plt
import matplotlib.style as mplstyle
import buencolors as bc


class TestPackageImport:
    """Test package-level imports and initialization."""

    def test_package_has_version(self):
        """Test that package has a version attribute."""
        assert hasattr(bc, '__version__')
        assert isinstance(bc.__version__, str)
        # Version should not be 'unknown' in installed package
        # (it might be 'unknown' during development)

    def test_all_functions_exported(self):
        """Test that all expected functions are exported."""
        expected_exports = [
            'eject_legend',
            'rotate_discrete_xticks',
            'grab_legend',
            'get_density',
            'shuffle',
            'number_to_color',
            'get_palette',
            'display_palette',
            'list_palettes',
            'get_registered_cmaps',
            'clean_umap'
        ]

        for func_name in expected_exports:
            assert hasattr(bc, func_name), f"Function {func_name} not exported"

    def test_import_registers_colormaps(self):
        """Test that importing buencolors registers all colormaps."""
        # Import should register colormaps
        from buencolors.palettes import jdb_palettes

        for palette_name in jdb_palettes.keys():
            assert palette_name in plt.colormaps()
            assert f"{palette_name}_r" in plt.colormaps()


class TestStyleRegistration:
    """Test matplotlib style registration."""

    def test_pretty_plot_style_available(self):
        """Test that 'pretty-plot' style is available after import."""
        available_styles = plt.style.available
        assert 'pretty-plot' in available_styles

    def test_pretty_plot_style_can_be_used(self, clean_plot):
        """Test that 'pretty-plot' style can be applied."""
        # This should not raise an error
        plt.style.use('pretty-plot')

        # Create a simple plot to verify style is applied
        fig, ax = plt.subplots()
        ax.plot([1, 2, 3], [1, 2, 3])

        # Check some style attributes
        # pretty-plot removes top and right spines
        assert not ax.spines['top'].get_visible()
        assert not ax.spines['right'].get_visible()
        assert ax.spines['bottom'].get_visible()
        assert ax.spines['left'].get_visible()

    def test_style_properties(self, clean_plot):
        """Test that pretty-plot style sets expected properties."""
        with plt.style.context('pretty-plot'):
            fig, ax = plt.subplots()

            # Check that background is transparent (none)
            assert ax.get_facecolor() == (0.0, 0.0, 0.0, 0.0) or ax.get_fc() == 'none'


class TestEndToEndWorkflow:
    """Test common end-to-end workflows."""

    def test_basic_plot_workflow(self, clean_plot):
        """Test a basic plotting workflow with buencolors."""
        import numpy as np

        # Use pretty-plot style
        plt.style.use('pretty-plot')

        # Create plot with BuenColors palette
        x = np.linspace(0, 10, 100)
        colors = bc.get_palette('Zissou', n=3)

        fig, ax = plt.subplots()
        ax.plot(x, np.sin(x), label='sin', color=colors[0])
        ax.plot(x, np.cos(x), label='cos', color=colors[1])
        ax.plot(x, np.sin(x) * np.cos(x), label='sin*cos', color=colors[2])

        ax.legend()
        bc.eject_legend(ax)

        # Verify plot was created successfully
        assert len(ax.lines) == 3
        assert ax.get_legend() is not None

    def test_colormap_workflow(self, clean_plot):
        """Test using BuenColors colormaps with matplotlib."""
        import numpy as np

        plt.style.use('pretty-plot')

        # Create heatmap with BuenColors colormap
        data = np.random.randn(10, 10)
        fig, ax = plt.subplots()
        im = ax.imshow(data, cmap='solar_flare')

        # Verify colormap is applied
        assert im.get_cmap().name == 'solar_flare'

    def test_number_to_color_workflow(self, clean_plot):
        """Test mapping numbers to colors workflow."""
        import numpy as np

        np.random.seed(42)
        values = np.random.randn(50) ** 2
        x = np.arange(50)

        colors = bc.number_to_color(values, 'brewer_spectra')

        fig, ax = plt.subplots()
        ax.scatter(x, values, c=colors)

        # Verify scatter plot created
        assert len(ax.collections) == 1

    def test_density_coloring_workflow(self, clean_plot):
        """Test density-based coloring workflow."""
        import numpy as np

        np.random.seed(42)
        x = np.random.randn(1000)
        y = np.random.randn(1000)

        density = bc.get_density(x, y)
        colors = bc.number_to_color(density, 'viridis')

        # Shuffle for better visualization
        import pandas as pd
        df = pd.DataFrame({'x': x, 'y': y, 'color': colors})
        df = bc.shuffle(df)

        fig, ax = plt.subplots()
        ax.scatter(df['x'], df['y'], c=df['color'], s=5, alpha=0.5)

        # Verify plot created
        assert len(ax.collections) == 1

    def test_legend_extraction_workflow(self, clean_plot):
        """Test extracting legend to separate figure."""
        import numpy as np

        plt.style.use('pretty-plot')

        x = np.linspace(0, 10, 100)
        fig, ax = plt.subplots()
        ax.plot(x, np.sin(x), label='sin(x)')
        ax.plot(x, np.cos(x), label='cos(x)')
        ax.legend()

        # Extract legend
        legend_fig = bc.grab_legend(ax, remove=True)

        # Verify legend was extracted
        assert isinstance(legend_fig, plt.Figure)
        assert ax.get_legend() is None  # Removed from original

    def test_palette_display_workflow(self, clean_plot):
        """Test displaying and exploring palettes."""
        # List palettes
        palettes = bc.list_palettes()
        assert len(palettes) > 100

        # Filter palettes
        flame_palettes = bc.list_palettes('flame')
        assert all('flame' in p.lower() for p in flame_palettes)

        # Display a palette
        fig = bc.display_palette('solar_flare')
        assert isinstance(fig, plt.Figure)

        # Get palette colors
        colors = bc.get_palette('Zissou')
        assert len(colors) == 5