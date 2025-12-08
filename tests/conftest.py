import pytest
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend for testing
import matplotlib.pyplot as plt


@pytest.fixture
def sample_data():
    """Generate sample data for testing."""
    np.random.seed(42)
    return {
        'x': np.random.randn(100),
        'y': np.random.randn(100),
        'values': np.random.rand(100),
        'categories': np.random.choice(['A', 'B', 'C'], 100)
    }


@pytest.fixture
def sample_dataframe(sample_data):
    """Create a sample DataFrame for testing."""
    return pd.DataFrame(sample_data)


@pytest.fixture
def clean_plot():
    """Fixture to clean up matplotlib figures after each test."""
    yield
    plt.close('all')


@pytest.fixture(autouse=True)
def reset_random_state():
    """Reset random state before each test for reproducibility."""
    np.random.seed(42)