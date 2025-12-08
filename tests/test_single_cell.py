import pytest
import numpy as np
import matplotlib.pyplot as plt

# Try to import single-cell dependencies
try:
    import anndata as ad
    import scanpy as sc
    from buencolors.single_cell import clean_umap
    SCANPY_AVAILABLE = True
except ImportError:
    SCANPY_AVAILABLE = False


@pytest.mark.skipif(not SCANPY_AVAILABLE, reason="scanpy or anndata not installed")
class TestCleanUMAP:
    """Test clean_umap function (requires scanpy and anndata)."""

    @pytest.fixture
    def mock_adata(self):
        """Create a mock AnnData object with UMAP coordinates."""
        # Create synthetic data
        n_obs = 100
        n_vars = 50

        # Random count matrix
        X = np.random.negative_binomial(5, 0.3, (n_obs, n_vars))

        # Create AnnData
        adata = ad.AnnData(X)
        adata.obs['cell_type'] = np.random.choice(['TypeA', 'TypeB', 'TypeC'], n_obs)
        adata.var_names = [f'Gene_{i}' for i in range(n_vars)]

        # Add UMAP coordinates
        np.random.seed(42)
        adata.obsm['X_umap'] = np.random.randn(n_obs, 2)

        return adata

    def test_clean_umap_returns_axes(self, mock_adata, clean_plot):
        """Test that clean_umap returns matplotlib axes."""
        # Note: clean_umap already sets show=False internally
        ax = clean_umap(mock_adata, color='cell_type')

        assert ax is not None
        assert isinstance(ax, plt.Axes)

    def test_clean_umap_removes_decorations(self, mock_adata, clean_plot):
        """Test that clean_umap removes standard plot decorations."""
        ax = clean_umap(mock_adata, color='cell_type')

        # Check that standard decorations are removed
        assert ax.get_xlabel() == ""
        assert ax.get_ylabel() == ""
        assert len(ax.get_xticks()) == 0 or all(not label.get_visible() for label in ax.get_xticklabels())
        assert len(ax.get_yticks()) == 0 or all(not label.get_visible() for label in ax.get_yticklabels())

        # Check that spines are hidden
        for spine in ax.spines.values():
            assert not spine.get_visible()

    def test_clean_umap_has_legend(self, mock_adata, clean_plot):
        """Test that clean_umap creates a legend."""
        ax = clean_umap(mock_adata, color='cell_type')

        legend = ax.get_legend()
        assert legend is not None

    def test_clean_umap_custom_axis_length(self, mock_adata, clean_plot):
        """Test clean_umap with custom axis length."""
        ax = clean_umap(mock_adata, color='cell_type', axis_len=0.3)

        assert ax is not None
        assert isinstance(ax, plt.Axes)

    def test_clean_umap_custom_thickness(self, mock_adata, clean_plot):
        """Test clean_umap with custom line thickness."""
        ax = clean_umap(mock_adata, color='cell_type', thickness=5.0)

        assert ax is not None

    def test_clean_umap_with_existing_axis(self, mock_adata, clean_plot):
        """Test clean_umap with pre-existing axis."""
        fig, ax = plt.subplots()
        result_ax = clean_umap(mock_adata, color='cell_type', ax=ax)

        assert result_ax is ax

    def test_clean_umap_with_gene_color(self, mock_adata, clean_plot):
        """Test clean_umap colored by gene expression."""
        # Add some gene expression
        mock_adata.X = np.random.rand(mock_adata.n_obs, mock_adata.n_vars)

        ax = clean_umap(mock_adata, color='Gene_0')

        assert ax is not None

    def test_clean_umap_with_kwargs(self, mock_adata, clean_plot):
        """Test clean_umap with additional kwargs passed to sc.pl.umap."""
        ax = clean_umap(mock_adata, color='cell_type', size=50)

        assert ax is not None


@pytest.mark.skipif(SCANPY_AVAILABLE, reason="Test only when scanpy is not available")
def test_clean_umap_not_available():
    """Test that clean_umap is not available when dependencies are missing."""
    try:
        from buencolors.single_cell import clean_umap
        # If we get here, the import succeeded when it shouldn't have
        pytest.fail("clean_umap should not be available without scanpy")
    except (ImportError, AttributeError):
        # This is expected
        pass