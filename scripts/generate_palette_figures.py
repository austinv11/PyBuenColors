#!/usr/bin/env python3
"""
Generate palette visualization figures for README.

This script creates a comprehensive visualization of all BuenColors palettes.

Usage:
    python scripts/generate_palette_figures.py

Output:
    - figures/all_palettes.png - All 106 palettes in one figure
"""
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
import buencolors


def create_output_dir():
    """Create figures directory if it doesn't exist."""
    figures_dir = Path("figures")
    figures_dir.mkdir(exist_ok=True)
    return figures_dir


def generate_all_palettes_figure(output_dir):
    """Generate a single figure showing all 106 palettes."""
    # Get all registered colormap names (without _r variants)
    all_palettes = buencolors.get_registered_cmaps(include_reversed=False)

    n_palettes = len(all_palettes)
    print(f"Generating figure for {n_palettes} palettes...")

    # Create figure with one row per palette
    fig, axes = plt.subplots(n_palettes, 1, figsize=(12, n_palettes * 0.4))

    # Handle single palette case
    if n_palettes == 1:
        axes = [axes]

    for ax, palette_name in zip(axes, all_palettes):
        # Get colors (interpolate to 256 for smooth gradient)
        colors = buencolors.get_palette(palette_name, 256)

        # Create color bar
        ax.imshow([list(range(len(colors)))], aspect='auto',
                 cmap=matplotlib.colors.ListedColormap(colors),
                 interpolation='bilinear')

        # Configure axes
        ax.set_yticks([])
        ax.set_xticks([])
        ax.set_ylabel(palette_name, fontsize=9, rotation=0, ha='right',
                     va='center', fontfamily='monospace')

        # Remove spines
        for spine in ax.spines.values():
            spine.set_visible(False)

    plt.suptitle(f'BuenColors: All {n_palettes} Palettes',
                fontsize=16, fontweight='bold', y=0.9995)
    plt.tight_layout()

    output_path = output_dir / 'all_palettes.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"✓ Generated: {output_path}")
    plt.close()


def main():
    """Main execution function."""
    print("BuenColors Palette Figure Generator")
    print("=" * 50)

    # Create output directory
    output_dir = create_output_dir()
    print(f"\nOutput directory: {output_dir}/")
    print()

    # Generate the comprehensive palette figure
    generate_all_palettes_figure(output_dir)

    print()
    print("=" * 50)
    print(f"✓ Figure generated successfully!")
    print(f"✓ Saved to: {output_dir.absolute()}/")
    print()
    print("Usage in README.md:")
    print("  ![All Palettes](figures/all_palettes.png)")
    print()


if __name__ == "__main__":
    main()