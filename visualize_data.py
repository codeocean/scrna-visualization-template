#!/usr/bin/env python3
"""
Single Cell Data Visualization Script
Usage: python visualize_data.py [threshold] [title] [cell_type]
Example: python visualize_data.py 2.5 "My Experiment" t_cells

If no arguments provided, uses defaults:
  threshold = 2.0
  title = "Single Cell Analysis"
  cell_type = "t_cells"
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Hardcoded data path
DATA_PATH = "/data/experiment/single_cell_data.npz"

# Default values
DEFAULT_THRESHOLD = 2.0
DEFAULT_TITLE = "Single Cell Analysis"
DEFAULT_CELL_TYPE = "t_cells"

# Set style for prettier plots
sns.set_style("whitegrid")
plt.rcParams['figure.facecolor'] = 'white'
plt.rcParams['axes.facecolor'] = '#f8f9fa'

def load_data(filepath):
    """Load single cell data from npz file"""
    try:
        data = np.load(filepath, allow_pickle=True)
        return {
            'x': data['x'],
            'y': data['y'],
            'cell_types': data['cell_types'],
            'gene1': data['gene1'],
            'gene2': data['gene2'],
            'colors': data['colors']
        }
    except FileNotFoundError:
        print(f"\n❌ Error: File '{filepath}' not found!")
        print("Please generate data first using: python generate_data.py")
        exit(1)
    except Exception as e:
        print(f"\n❌ Error loading data: {e}")
        exit(1)

def plot_single_cell_data(data, threshold, title, highlighted_type):
    """Create a beautiful visualization of single cell data"""
    fig, axes = plt.subplots(1, 2, figsize=(16, 7))
    fig.suptitle(title, fontsize=20, fontweight='bold', y=0.98)
    
    x = data['x']
    y = data['y']
    cell_types = data['cell_types']
    gene1 = data['gene1']
    gene2 = data['gene2']
    colors = data['colors']
    
    # Create highlight mask
    highlight_mask = cell_types == highlighted_type
    
    # Plot 1: UMAP-like cell type clustering
    ax1 = axes[0]
    
    # Plot non-highlighted cells first
    ax1.scatter(x[~highlight_mask], y[~highlight_mask], 
               c=colors[~highlight_mask], s=50, alpha=0.6, 
               edgecolors='white', linewidth=0.5)
    
    # Plot highlighted cells on top
    if np.any(highlight_mask):
        ax1.scatter(x[highlight_mask], y[highlight_mask], 
                   c=colors[highlight_mask], s=120, alpha=0.9, 
                   edgecolors='gold', linewidth=2.5, zorder=100)
    
    ax1.set_xlabel('UMAP 1', fontsize=14, fontweight='bold')
    ax1.set_ylabel('UMAP 2', fontsize=14, fontweight='bold')
    ax1.set_title('Cell Type Clustering', fontsize=16, fontweight='bold', pad=15)
    ax1.grid(alpha=0.3, linestyle='--')
    
    # Add legend
    cell_types_legend = {
        't_cells': ('T Cells', '#e74c3c'),
        'b_cells': ('B Cells', '#3498db'),
        'monocytes': ('Monocytes', '#2ecc71'),
        'nk_cells': ('NK Cells', '#f39c12')
    }
    
    legend_elements = []
    for ct_name, (label, color) in cell_types_legend.items():
        marker_size = 120 if ct_name == highlighted_type else 50
        edge_color = 'gold' if ct_name == highlighted_type else 'white'
        edge_width = 2.5 if ct_name == highlighted_type else 0.5
        legend_elements.append(
            plt.scatter([], [], c=color, s=marker_size, alpha=0.9,
                       edgecolors=edge_color, linewidth=edge_width,
                       label=label + (' (Selected)' if ct_name == highlighted_type else ''))
        )
    
    ax1.legend(handles=legend_elements, loc='upper right', frameon=True, 
              fancybox=True, shadow=True, fontsize=11)
    
    # Plot 2: Gene expression with threshold
    ax2 = axes[1]
    
    # Scatter plot of gene expression
    ax2.scatter(gene1, gene2, c=colors, s=80, alpha=0.6,
               edgecolors='white', linewidth=0.5)
    
    # Add threshold lines
    ax2.axhline(y=threshold, color='red', linestyle='--', linewidth=2.5, 
               label=f'Threshold = {threshold}', alpha=0.8)
    ax2.axvline(x=threshold, color='red', linestyle='--', linewidth=2.5, alpha=0.8)
    
    # Highlight cells above threshold
    above_threshold = (gene1 > threshold) & (gene2 > threshold)
    if np.any(above_threshold):
        ax2.scatter(gene1[above_threshold], gene2[above_threshold], 
                   s=200, facecolors='none', edgecolors='red', 
                   linewidth=2.5, alpha=0.7, zorder=99)
    
    ax2.set_xlabel('Gene 1 Expression', fontsize=14, fontweight='bold')
    ax2.set_ylabel('Gene 2 Expression', fontsize=14, fontweight='bold')
    ax2.set_title(f'Gene Expression (Threshold: {threshold})', 
                 fontsize=16, fontweight='bold', pad=15)
    ax2.grid(alpha=0.3, linestyle='--')
    ax2.legend(loc='upper right', frameon=True, fancybox=True, shadow=True, fontsize=11)
    
    # Count cells above threshold
    n_above = np.sum(above_threshold)
    ax2.text(0.02, 0.98, f'Cells above threshold: {n_above}/{len(gene1)}', 
            transform=ax2.transAxes, fontsize=11, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.7))
    
    plt.tight_layout()
    return fig

def main():
    # Parse command line arguments with defaults
    if len(sys.argv) > 1:
        threshold = float(sys.argv[1])
    else:
        threshold = DEFAULT_THRESHOLD
    
    if len(sys.argv) > 2:
        title = sys.argv[2]
    else:
        title = DEFAULT_TITLE
    
    if len(sys.argv) > 3:
        cell_type = sys.argv[3]
        valid_types = ['t_cells', 'b_cells', 'monocytes', 'nk_cells']
        if cell_type not in valid_types:
            print(f"\n❌ Error: Invalid cell type '{cell_type}'")
            print(f"Valid options: {', '.join(valid_types)}")
            exit(1)
    else:
        cell_type = DEFAULT_CELL_TYPE
    
    # Display parameters
    print("\n" + "="*60)
    print("Single Cell Data Visualization")
    print("="*60)
    print(f"Data file:        {DATA_PATH}")
    print(f"Threshold:        {threshold}")
    print(f"Title:            {title}")
    print(f"Highlighted Type: {cell_type.replace('_', ' ').title()}")
    print("="*60 + "\n")
    
    # Load data
    print("Loading data...")
    data = load_data(DATA_PATH)
    print(f"✓ Loaded {len(data['x'])} cells")
    
    # Create visualization
    print("Creating visualization...")
    fig = plot_single_cell_data(data, threshold, title, cell_type)
    
    # Save plot
    #output_file = f"/results/single_cell_viz_{cell_type}_threshold_{threshold}.png"
    output_file = f"/results/single_cell_viz.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"✓ Plot saved as: {output_file}")
    
    plt.show()
    print("\nDone!")

if __name__ == "__main__":
    main()
