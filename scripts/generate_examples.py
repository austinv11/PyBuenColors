#!/usr/bin/env python3
"""
Generate example figures showcasing BuenColors style and utilities.

This script creates minimal, high-quality figures demonstrating:
- The pretty-plot matplotlib style
- Helper function utilities

Usage:
    python scripts/generate_examples.py

Output:
    - figures/default_style.png - Example with default matplotlib style
    - figures/pretty_style.png - Example with pretty-plot style
    - figures/util_*.png - Individual utility function examples
"""

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import pandas as pd
from pathlib import Path
import buencolors
from matplotlib.patches import Patch


def create_output_dir():
    """Create figures directory if it doesn't exist."""
    figures_dir = Path("figures")
    figures_dir.mkdir(exist_ok=True)
    return figures_dir


def generate_default_style_example(output_dir):
    """Generate example with default matplotlib style."""
    print("Generating default style example...")

    # Prepare data
    np.random.seed(42)
    x = np.linspace(0, 10, 100)
    y1 = np.sin(x) + np.random.normal(0, 0.1, 100)
    y2 = np.cos(x) + np.random.normal(0, 0.1, 100)
    y3 = np.sin(x) * np.cos(x) + np.random.normal(0, 0.1, 100)

    # Use default matplotlib style
    with mpl.style.context('default'):
        fig, ax = plt.subplots(figsize=(10, 6))

        ax.plot(x, y1, label='sin(x)', linewidth=2)
        ax.plot(x, y2, label='cos(x)', linewidth=2)
        ax.plot(x, y3, label='sin(x)·cos(x)', linewidth=2)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('Default matplotlib Style', fontsize=14, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)

        output_path = output_dir / 'default_style.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"✓ Generated: {output_path}")
        plt.close()


def generate_pretty_style_example(output_dir):
    """Generate example with pretty-plot style."""
    print("Generating pretty-plot style example...")

    # Prepare data (same as default for fair comparison)
    np.random.seed(42)
    x = np.linspace(0, 10, 100)
    y1 = np.sin(x) + np.random.normal(0, 0.1, 100)
    y2 = np.cos(x) + np.random.normal(0, 0.1, 100)
    y3 = np.sin(x) * np.cos(x) + np.random.normal(0, 0.1, 100)

    # Use pretty-plot style
    with plt.style.context('pretty-plot'):
        fig, ax = plt.subplots(figsize=(10, 6))

        ax.plot(x, y1, label='sin(x)', linewidth=2)
        ax.plot(x, y2, label='cos(x)', linewidth=2)
        ax.plot(x, y3, label='sin(x)·cos(x)', linewidth=2)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('BuenColors pretty-plot Style', fontsize=14, fontweight='bold')
        ax.legend()

        output_path = output_dir / 'pretty_style.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"✓ Generated: {output_path}")
        plt.close()


def generate_eject_legend_example(output_dir):
    """Generate eject_legend() utility example."""
    print("Generating eject_legend example...")

    np.random.seed(42)
    with plt.style.context('pretty-plot'):
        fig, ax = plt.subplots(figsize=(10, 6))
        x = np.linspace(0, 10, 100)
        ax.plot(x, np.sin(x), label='sin(x)', linewidth=2.5)
        ax.plot(x, np.cos(x), label='cos(x)', linewidth=2.5)
        ax.plot(x, np.sin(x) * np.cos(x), label='sin(x)·cos(x)', linewidth=2.5)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('eject_legend() - Legend Moved Outside Plot', fontweight='bold')
        buencolors.eject_legend(ax)

        output_path = output_dir / 'util_eject_legend.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"✓ Generated: {output_path}")
        plt.close()


def generate_rotate_xticks_example(output_dir):
    """Generate rotate_discrete_xticks() utility example."""
    print("Generating rotate_discrete_xticks example...")

    np.random.seed(42)
    with plt.style.context('pretty-plot'):
        fig, ax = plt.subplots(figsize=(10, 6))
        categories = ['Machine Learning', 'Deep Learning', 'Computer Vision',
                      'Natural Language Processing', 'Reinforcement Learning']
        values = [85, 92, 78, 88, 75]
        colors = buencolors.get_palette('Royal2', len(categories))
        bars = ax.bar(categories, values, color=colors, edgecolor='black', linewidth=1.2)
        ax.set_ylabel('Performance Score')
        ax.set_title('rotate_discrete_xticks() - Rotated Labels for Readability', fontweight='bold')
        buencolors.rotate_discrete_xticks(ax, rotation=45)

        output_path = output_dir / 'util_rotate_xticks.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"✓ Generated: {output_path}")
        plt.close()


def generate_grab_legend_example(output_dir):
    """Generate grab_legend() utility example."""
    print("Generating grab_legend example...")

    np.random.seed(42)
    with plt.style.context('pretty-plot'):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
        x = np.linspace(0, 10, 100)

        # Left plot: remove=True (default)
        ax1.plot(x, np.sin(x), label='sin(x)', linewidth=2.5)
        ax1.plot(x, np.cos(x), label='cos(x)', linewidth=2.5)
        ax1.plot(x, np.sin(x) * np.cos(x), label='sin(x)·cos(x)', linewidth=2.5)
        ax1.legend(loc='upper right')
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_title('Before grab_legend(remove=True)', fontweight='bold')

        # Extract legend (removes it from ax1)
        legend_fig1 = buencolors.grab_legend(ax1, remove=True)

        # Right plot: remove=False
        ax2.plot(x, np.sin(x), label='sin(x)', linewidth=2.5)
        ax2.plot(x, np.cos(x), label='cos(x)', linewidth=2.5)
        ax2.plot(x, np.sin(x) * np.cos(x), label='sin(x)·cos(x)', linewidth=2.5)
        ax2.legend(loc='upper right')
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_title('grab_legend(remove=False)', fontweight='bold')

        # Extract legend (keeps it on ax2)
        legend_fig2 = buencolors.grab_legend(ax2, remove=False)

        # plt.suptitle('grab_legend() - Extract Legend to Separate Figure',
        #             fontsize=14, fontweight='bold', y=1.00)

        output_path = output_dir / 'util_grab_legend.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"✓ Generated: {output_path}")
        plt.close()

        # Clean up the legend figures
        plt.close(legend_fig1)
        plt.close(legend_fig2)


def generate_density_example(output_dir):
    """Generate get_density() utility example."""
    print("Generating get_density example...")

    np.random.seed(42)
    with plt.style.context('pretty-plot'):
        fig, ax = plt.subplots(figsize=(10, 8))
        n = 2000
        x = np.concatenate([
            np.random.normal(0, 0.5, n // 2),
            np.random.normal(3, 0.8, n // 2)
        ])
        y = np.concatenate([
            np.random.normal(0, 0.5, n // 2),
            np.random.normal(3, 0.8, n // 2)
        ])
        density = buencolors.get_density(x, y, n=200)
        scatter = ax.scatter(x, y, c=density, cmap='ocean_earth',
                           s=12, alpha=0.6, edgecolors='none')
        plt.colorbar(scatter, ax=ax, label='Point Density')
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('get_density() - Density-Based Scatter Plot Coloring', fontweight='bold')

        output_path = output_dir / 'util_density.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"✓ Generated: {output_path}")
        plt.close()


def generate_shuffle_example(output_dir):
    """Generate shuffle() utility example."""
    print("Generating shuffle example...")

    np.random.seed(42)
    with plt.style.context('pretty-plot'):
        fig, ax = plt.subplots(figsize=(10, 8))
        df = pd.DataFrame({
            'x': np.random.randn(3000),
            'y': np.random.randn(3000),
            'category': np.repeat(['Group A', 'Group B', 'Group C'], 1000)
        })
        color_map = {'Group A': '#E41A1C', 'Group B': '#377EB8', 'Group C': '#4DAF4A'}
        df['color'] = df['category'].map(color_map)
        df_shuffled = buencolors.shuffle(df)
        ax.scatter(df_shuffled['x'], df_shuffled['y'], c=df_shuffled['color'],
                 s=40, alpha=0.7, edgecolors='black', linewidth=0.3)
        ax.set_xlabel('x')
        ax.set_ylabel('y')
        ax.set_title('shuffle() - Randomized Point Order for Fair Visualization', fontweight='bold')

        # Add legend
        legend_elements = [Patch(facecolor=color_map[cat], label=cat, edgecolor='black')
                          for cat in ['Group A', 'Group B', 'Group C']]
        ax.legend(handles=legend_elements, loc='upper right')

        output_path = output_dir / 'util_shuffle.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"✓ Generated: {output_path}")
        plt.close()


def generate_number_to_color_example(output_dir):
    """Generate number_to_color() utility example."""
    print("Generating number_to_color example...")

    np.random.seed(42)
    with plt.style.context('pretty-plot'):
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))

        n_points = 500
        values = np.random.randn(n_points) ** 2
        x = np.random.randn(n_points)
        y = np.random.randn(n_points)

        # Continuous
        colors_continuous = buencolors.number_to_color(values, 'solar_flare')
        ax1.scatter(x, y, c=colors_continuous, s=50, alpha=0.7, edgecolors='black', linewidth=0.3)
        ax1.set_xlabel('x')
        ax1.set_ylabel('y')
        ax1.set_title('Continuous Color Mapping', fontweight='bold')

        # Discrete
        colors_discrete = buencolors.number_to_color(values, 'brewer_spectra', n_bins=8)
        ax2.scatter(x, y, c=colors_discrete, s=50, alpha=0.7, edgecolors='black', linewidth=0.3)
        ax2.set_xlabel('x')
        ax2.set_ylabel('y')
        ax2.set_title('Discrete Binning (8 bins)', fontweight='bold')

        plt.suptitle('number_to_color() - Map Numeric Values to Colors',
                    fontsize=14, fontweight='bold', y=1.02)

        output_path = output_dir / 'util_number_to_color.png'
        plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
        print(f"✓ Generated: {output_path}")
        plt.close()


def main():
    """Main execution function."""
    print("BuenColors Example Figure Generator")
    print("=" * 50)

    # Create output directory
    output_dir = create_output_dir()
    print(f"\nOutput directory: {output_dir}/")
    print()

    # Generate figures
    generate_default_style_example(output_dir)
    generate_pretty_style_example(output_dir)
    generate_eject_legend_example(output_dir)
    generate_rotate_xticks_example(output_dir)
    generate_grab_legend_example(output_dir)
    generate_density_example(output_dir)
    generate_shuffle_example(output_dir)
    generate_number_to_color_example(output_dir)

    print()
    print("=" * 50)
    print(f"✓ All figures generated successfully!")
    print(f"✓ Saved to: {output_dir.absolute()}/")
    print()
    print("Usage in README.md:")
    print("  ![Default Style](figures/default_style.png)")
    print("  ![Pretty-Plot Style](figures/pretty_style.png)")
    print("  ![eject_legend](figures/util_eject_legend.png)")
    print("  ![rotate_discrete_xticks](figures/util_rotate_xticks.png)")
    print("  ![grab_legend](figures/util_grab_legend.png)")
    print("  ![get_density](figures/util_density.png)")
    print("  ![shuffle](figures/util_shuffle.png)")
    print("  ![number_to_color](figures/util_number_to_color.png)")
    print()


if __name__ == "__main__":
    main()