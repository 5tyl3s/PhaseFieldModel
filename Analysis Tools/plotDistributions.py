"""
Analysis script for grain structure using pole figures and size distributions.
Analyzes grain map and orientations to:
1. Plot pole figure
2. Calculate grain sizes
3. Plot distribution of grain sizes per pixel (probability distribution)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.spatial.transform import Rotation
from scipy.ndimage import gaussian_filter
import os

# Set up paths
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(SCRIPT_DIR)
BUILD_DIR = os.path.join(PROJECT_ROOT, "build")

# File paths
GRAIN_GRID_FILE = os.path.join(BUILD_DIR, "GrainGrid.csv")
GRAIN_ORIENTATIONS_FILE = os.path.join(BUILD_DIR, "GrainOrientations.csv")


def load_data():
    """Load grain grid and orientation data."""
    print("Loading data...")
    
    # Load grain grid - it appears to be space-separated or comma-separated
    grain_grid = np.loadtxt(GRAIN_GRID_FILE, delimiter=',', dtype=int)
    print(f"Grain grid shape: {grain_grid.shape}")
    
    # Load grain orientations
    orientations_df = pd.read_csv(GRAIN_ORIENTATIONS_FILE)
    print(f"Loaded {len(orientations_df)} grain orientations")
    print(f"Columns: {list(orientations_df.columns)}")
    
    return grain_grid, orientations_df


def euler_to_pole_vectors(theta1, Phi, theta2):
    """
    Convert Euler angles (in degrees) to pole figure vectors.
    Uses ZXZ convention (standard in materials science).
    Returns unit vectors pointing along the <100> direction (pole figure coordinates).
    """
    # Convert to radians
    theta1_rad = np.radians(theta1)
    Phi_rad = np.radians(Phi)
    theta2_rad = np.radians(theta2)
    
    # Create rotation from Euler angles (ZXZ convention)
    rot = Rotation.from_euler('ZXZ', [theta1_rad, Phi_rad, theta2_rad])
    
    # The pole figure typically shows the <100> orientation
    # Initial vector along Z-axis
    initial_pole = np.array([0, 0, 1])
    
    # Rotate the pole
    rotated_pole = rot.apply(initial_pole)
    
    return rotated_pole


def create_pole_figure(grain_grid, orientations_df, output_dir=None):
    """
    Create a stereographic projection pole figure heatmap showing per-pixel 
    orientation distribution (not per-grain).
    """
    print("\nCreating pole figure heatmap...")
    
    # Create a mapping from grain_id to orientation
    grain_orientation_map = {}
    for idx, row in orientations_df.iterrows():
        grain_id = int(row['grain_id'])
        grain_orientation_map[grain_id] = (row['theta1_deg'], row['Phi_deg'], row['theta2_deg'])
    
    # Calculate pole vectors for each pixel
    poles = []
    for i in range(grain_grid.shape[0]):
        for j in range(grain_grid.shape[1]):
            grain_id = grain_grid[i, j]
            if grain_id in grain_orientation_map:
                theta1, Phi, theta2 = grain_orientation_map[grain_id]
                pole = euler_to_pole_vectors(theta1, Phi, theta2)
                poles.append(pole)
    
    poles = np.array(poles)
    
    # Convert Cartesian coordinates to stereographic projection coordinates
    x, y, z = poles[:, 0], poles[:, 1], poles[:, 2]
    
    # Stereographic projection
    rho = np.sqrt(x**2 + y**2) / (1 + np.abs(z))
    theta = np.arctan2(y, x)
    
    # Convert polar to Cartesian for 2D histogram
    proj_x = rho * np.cos(theta)
    proj_y = rho * np.sin(theta)
    
    # Create 2D histogram heatmap with fine binning
    hist, xedges, yedges = np.histogram2d(proj_x, proj_y, bins=300, range=[[-1.5, 1.5], [-1.5, 1.5]])
    
    # Apply Gaussian filter for smooth interpolation
    hist_smooth = gaussian_filter(hist.T, sigma=2.0)
    
    # Create figure
    fig, ax = plt.subplots(figsize=(12, 12))
    
    # Plot smoothed heatmap
    extent = [-1.5, 1.5, -1.5, 1.5]
    im = ax.imshow(hist_smooth, extent=extent, origin='lower', cmap='hot', aspect='auto', interpolation='bilinear')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Pixel Count (Smoothed)', fontsize=12, fontweight='bold')
    
    # Add reference circle (equal area projection circle)
    circle = plt.Circle((0, 0), 1, fill=False, edgecolor='cyan', linestyle='--', linewidth=2, alpha=0.7)
    ax.add_patch(circle)
    
    ax.set_xlim(-1.5, 1.5)
    ax.set_ylim(-1.5, 1.5)
    ax.set_aspect('equal')
    
    ax.set_xlabel('X (Stereographic)', fontsize=12, fontweight='bold')
    ax.set_ylabel('Y (Stereographic)', fontsize=12, fontweight='bold')
    ax.set_title('Pole Figure Heatmap\n<100> Orientation Distribution (Per-Pixel, Interpolated)', fontsize=14, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    if output_dir is None:
        output_dir = os.path.dirname(GRAIN_ORIENTATIONS_FILE)
    
    output_file = os.path.join(output_dir, "pole_figure.png")
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Pole figure heatmap with interpolation saved to {output_file}")
    
    plt.close(fig)


def calculate_grain_sizes(grain_grid):
    """
    Calculate the size (number of pixels) for each grain.
    Returns a dictionary mapping grain_id to pixel_count.
    """
    print("\nCalculating grain sizes...")
    
    unique_grains, counts = np.unique(grain_grid, return_counts=True)
    grain_sizes = dict(zip(unique_grains, counts))
    
    print(f"Found {len(grain_sizes)} unique grains")
    print(f"Grain size range: {min(counts)} to {max(counts)} pixels")
    print(f"Mean grain size: {np.mean(counts):.1f} pixels")
    print(f"Median grain size: {np.median(counts):.1f} pixels")
    
    return grain_sizes


def plot_grain_size_distribution_per_pixel(grain_grid, grain_sizes, output_dir=None):
    """
    Plot the probability distribution of grain sizes per pixel.
    For each pixel, determine the size of the grain it belongs to,
    then plot the distribution of these grain sizes across all pixels.
    """
    print("\nPlotting grain size distribution per pixel...")
    
    # Get the grain size for each pixel
    pixel_grain_sizes = np.array([[grain_sizes[gid] for gid in row] for row in grain_grid])
    
    # Flatten to get all grain sizes for all pixels
    all_pixel_grain_sizes = pixel_grain_sizes.flatten()
    
    # Create histogram
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Plot 1: Histogram of grain sizes per pixel
    ax = axes[0]
    counts, bins, patches = ax.hist(all_pixel_grain_sizes, bins=50, edgecolor='black', alpha=0.7, density=True)
    ax.set_xlabel('Grain Size (pixels)', fontsize=12)
    ax.set_ylabel('Probability Density', fontsize=12)
    ax.set_title('Distribution of Grain Sizes Per Pixel\n(Probability Density)', fontsize=12)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Add statistics text
    stats_text = f'Mean: {np.mean(all_pixel_grain_sizes):.1f}\n'
    stats_text += f'Median: {np.median(all_pixel_grain_sizes):.1f}\n'
    stats_text += f'Std Dev: {np.std(all_pixel_grain_sizes):.1f}\n'
    stats_text += f'Min: {np.min(all_pixel_grain_sizes)}\n'
    stats_text += f'Max: {np.max(all_pixel_grain_sizes)}'
    ax.text(0.98, 0.97, stats_text, transform=ax.transAxes, 
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
            fontsize=10, family='monospace')
    
    # Plot 2: Cumulative distribution
    ax = axes[1]
    sorted_sizes = np.sort(all_pixel_grain_sizes)
    cumulative = np.arange(1, len(sorted_sizes) + 1) / len(sorted_sizes)
    ax.plot(sorted_sizes, cumulative, linewidth=2, color='blue')
    ax.set_xlabel('Grain Size (pixels)', fontsize=12)
    ax.set_ylabel('Cumulative Probability', fontsize=12)
    ax.set_title('Cumulative Distribution of Grain Sizes Per Pixel', fontsize=12)
    ax.grid(True, alpha=0.3)
    
    if output_dir is None:
        output_dir = os.path.dirname(GRAIN_ORIENTATIONS_FILE)
    
    output_file = os.path.join(output_dir, "grain_size_distribution.png")
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Grain size distribution plot saved to {output_file}")
    
    plt.close(fig)
    
    return all_pixel_grain_sizes


def plot_grain_size_map(grain_grid, grain_sizes, output_dir=None):
    """
    Create a visual map showing grain sizes in the grid.
    """
    print("\nCreating grain size map...")
    
    # Create array showing grain size at each pixel
    grain_size_map = np.zeros_like(grain_grid, dtype=float)
    for i in range(grain_grid.shape[0]):
        for j in range(grain_grid.shape[1]):
            grain_id = grain_grid[i, j]
            grain_size_map[i, j] = grain_sizes[grain_id]
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    im = ax.imshow(grain_size_map, cmap='viridis', origin='lower', aspect='auto')
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('Grain Size (pixels)', fontsize=12)
    
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_title('Grain Size Map', fontsize=14)
    
    if output_dir is None:
        output_dir = os.path.dirname(GRAIN_ORIENTATIONS_FILE)
    
    output_file = os.path.join(output_dir, "grain_size_map.png")
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Grain size map saved to {output_file}")
    
    plt.close(fig)


def plot_grain_id_map(grain_grid, output_dir=None):
    """
    Create a visual map showing grain IDs with different colors.
    """
    print("\nCreating grain ID map...")
    
    fig, ax = plt.subplots(figsize=(12, 10))
    
    # Create random colormap for grains
    im = ax.imshow(grain_grid, cmap='tab20', origin='lower', aspect='auto', interpolation='nearest')
    cbar = plt.colorbar(im, ax=ax, label='Grain ID')
    
    ax.set_xlabel('X', fontsize=12)
    ax.set_ylabel('Y', fontsize=12)
    ax.set_title('Grain ID Map', fontsize=14)
    
    if output_dir is None:
        output_dir = os.path.dirname(GRAIN_ORIENTATIONS_FILE)
    
    output_file = os.path.join(output_dir, "grain_id_map.png")
    fig.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Grain ID map saved to {output_file}")
    
    plt.close(fig)


def calculate_two_point_correlation_2d(grain_grid, max_distance=50):
    """
    Calculate 2D two-point correlation function.
    For each displacement (dx, dy), finds the probability that two points 
    separated by (dx, dy) belong to the same grain.
    
    Returns: 2D correlation array [dy, dx]
    """
    print("\nCalculating 2D two-point correlation function...")
    
    ny, nx = grain_grid.shape
    max_dist = min(max_distance, min(ny, nx) // 4)
    
    # Initialize 2D correlation array
    correlation_2d = np.zeros((2*max_dist + 1, 2*max_dist + 1))
    
    # Sample the grid to speed up calculation
    sample_step = max(1, min(ny, nx) // 40)  # Sample ~40 points per dimension
    
    # Compute correlation for all (dx, dy) pairs
    for dy in range(-max_dist, max_dist + 1):
        for dx in range(-max_dist, max_dist + 1):
            same_grain_count = 0
            total_count = 0
            
            # Check all point pairs with this separation
            for i in range(0, ny, sample_step):
                for j in range(0, nx, sample_step):
                    i2 = i + dy
                    j2 = j + dx
                    
                    # Check if target point is within bounds
                    if 0 <= i2 < ny and 0 <= j2 < nx:
                        if grain_grid[i, j] == grain_grid[i2, j2]:
                            same_grain_count += 1
                        total_count += 1
            
            if total_count > 0:
                correlation_2d[dy + max_dist, dx + max_dist] = same_grain_count / total_count
            else:
                correlation_2d[dy + max_dist, dx + max_dist] = 0
    
    print(f"2D correlation function computed for displacements ±{max_dist} pixels")
    print(f"Correlation at (0, 0): {correlation_2d[max_dist, max_dist]:.3f}")
    print(f"Correlation at (±1, 0): {correlation_2d[max_dist, max_dist+1]:.3f}")
    print(f"Correlation at (0, ±1): {correlation_2d[max_dist+1, max_dist]:.3f}")
    
    return correlation_2d, max_dist


def create_combined_analysis_figure(grain_grid, orientations_df, grain_sizes, correlation_2d, max_dist, output_dir=None):
    """
    Create a comprehensive combined figure with all analysis results in high resolution.
    """
    print("\nCreating combined high-resolution analysis figure...")
    
    # Create a large figure with subplots
    fig = plt.figure(figsize=(24, 20), dpi=150)  # Large figure for high resolution
    gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.3)
    
    # 1. Pole Figure - Per-Pixel Orientation Heatmap (top left)
    ax1 = fig.add_subplot(gs[0, 0])
    
    # Create a mapping from grain_id to orientation
    grain_orientation_map = {}
    for idx, row in orientations_df.iterrows():
        grain_id = int(row['grain_id'])
        grain_orientation_map[grain_id] = (row['theta1_deg'], row['Phi_deg'], row['theta2_deg'])
    
    # Calculate pole vectors for each pixel
    poles = []
    for i in range(grain_grid.shape[0]):
        for j in range(grain_grid.shape[1]):
            grain_id = grain_grid[i, j]
            if grain_id in grain_orientation_map:
                theta1, Phi, theta2 = grain_orientation_map[grain_id]
                pole = euler_to_pole_vectors(theta1, Phi, theta2)
                poles.append(pole)
    
    poles = np.array(poles)
    
    # Convert Cartesian to stereographic projection
    x, y, z = poles[:, 0], poles[:, 1], poles[:, 2]
    rho = np.sqrt(x**2 + y**2) / (1 + np.abs(z))
    theta_ang = np.arctan2(y, x)
    
    # Convert polar to Cartesian for 2D histogram
    proj_x = rho * np.cos(theta_ang)
    proj_y = rho * np.sin(theta_ang)
    
    # Create 2D histogram with fine binning
    hist, xedges, yedges = np.histogram2d(proj_x, proj_y, bins=250, range=[[-1.5, 1.5], [-1.5, 1.5]])
    
    # Apply Gaussian filter for smooth interpolation
    hist_smooth = gaussian_filter(hist.T, sigma=1.5)
    
    # Plot smoothed heatmap
    extent = [-1.5, 1.5, -1.5, 1.5]
    im1 = ax1.imshow(hist_smooth, extent=extent, origin='lower', cmap='hot', aspect='auto', interpolation='bilinear')
    cbar1 = plt.colorbar(im1, ax=ax1)
    cbar1.set_label('Pixel Count', fontsize=9, fontweight='bold')
    
    # Add reference circle
    circle = plt.Circle((0, 0), 1, fill=False, edgecolor='cyan', linestyle='--', linewidth=1.5, alpha=0.6)
    ax1.add_patch(circle)
    
    ax1.set_xlim(-1.5, 1.5)
    ax1.set_ylim(-1.5, 1.5)
    ax1.set_aspect('equal')
    
    ax1.set_xlabel('X (Stereo.)', fontsize=11, fontweight='bold')
    ax1.set_ylabel('Y (Stereo.)', fontsize=11, fontweight='bold')
    ax1.set_title('Pole Figure Heatmap (Interpolated)\n<100> Per-Pixel Orientation', fontsize=13, fontweight='bold')
    
    # 2. Grain Size Distribution Histogram (top middle)
    ax2 = fig.add_subplot(gs[0, 1])
    pixel_grain_sizes = np.array([[grain_sizes[gid] for gid in row] for row in grain_grid])
    all_pixel_grain_sizes = pixel_grain_sizes.flatten()
    
    counts, bins, patches = ax2.hist(all_pixel_grain_sizes, bins=60, edgecolor='black', alpha=0.7, color='steelblue')
    ax2.set_xlabel('Grain Size (pixels)', fontsize=12, fontweight='bold')
    ax2.set_ylabel('Frequency', fontsize=12, fontweight='bold')
    ax2.set_title('Grain Size Distribution\nProbability Density', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='y')
    
    # Add statistics
    stats_text = f'Mean: {np.mean(all_pixel_grain_sizes):.1f}\n'
    stats_text += f'Median: {np.median(all_pixel_grain_sizes):.1f}\n'
    stats_text += f'Std: {np.std(all_pixel_grain_sizes):.1f}'
    ax2.text(0.98, 0.97, stats_text, transform=ax2.transAxes, 
            verticalalignment='top', horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
            fontsize=10, family='monospace')
    
    # 3. Cumulative Distribution (top right)
    ax3 = fig.add_subplot(gs[0, 2])
    sorted_sizes = np.sort(all_pixel_grain_sizes)
    cumulative = np.arange(1, len(sorted_sizes) + 1) / len(sorted_sizes)
    ax3.plot(sorted_sizes, cumulative, linewidth=2.5, color='darkblue')
    ax3.fill_between(sorted_sizes, cumulative, alpha=0.3, color='steelblue')
    ax3.set_xlabel('Grain Size (pixels)', fontsize=12, fontweight='bold')
    ax3.set_ylabel('Cumulative Probability', fontsize=12, fontweight='bold')
    ax3.set_title('Cumulative Distribution\nGrain Sizes Per Pixel', fontsize=14, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    
    # 4. Grain Size Map (middle left and middle center)
    ax4 = fig.add_subplot(gs[1, :2])
    grain_size_map = np.zeros_like(grain_grid, dtype=float)
    for i in range(grain_grid.shape[0]):
        for j in range(grain_grid.shape[1]):
            grain_id = grain_grid[i, j]
            grain_size_map[i, j] = grain_sizes[grain_id]
    
    im4 = ax4.imshow(grain_size_map, cmap='viridis', origin='lower', aspect='auto')
    cbar4 = plt.colorbar(im4, ax=ax4)
    cbar4.set_label('Grain Size (pixels)', fontsize=11, fontweight='bold')
    ax4.set_xlabel('X (pixels)', fontsize=12, fontweight='bold')
    ax4.set_ylabel('Y (pixels)', fontsize=12, fontweight='bold')
    ax4.set_title('Spatial Grain Size Distribution', fontsize=14, fontweight='bold')
    
    # 5. 2D Two-Point Correlation (middle right)
    ax5 = fig.add_subplot(gs[1, 2])
    
    im5 = ax5.imshow(correlation_2d, cmap='RdYlBu_r', origin='lower', aspect='auto', 
                     extent=[-max_dist, max_dist, -max_dist, max_dist], vmin=0, vmax=1)
    cbar5 = plt.colorbar(im5, ax=ax5)
    cbar5.set_label('P(same grain)', fontsize=11, fontweight='bold')
    
    ax5.set_xlabel('ΔX (pixels)', fontsize=12, fontweight='bold')
    ax5.set_ylabel('ΔY (pixels)', fontsize=12, fontweight='bold')
    ax5.set_title('2D Two-Point Correlation\nSpatial Autocorrelation', fontsize=14, fontweight='bold')
    ax5.grid(True, alpha=0.2, color='white', linewidth=0.5)
    
    # Add crosshairs at origin
    ax5.axhline(y=0, color='white', linestyle='--', linewidth=1, alpha=0.5)
    ax5.axvline(x=0, color='white', linestyle='--', linewidth=1, alpha=0.5)
    
    # 6. Grain ID Map (bottom left and middle)
    ax6 = fig.add_subplot(gs[2, :2])
    im6 = ax6.imshow(grain_grid, cmap='tab20', origin='lower', aspect='auto', interpolation='nearest')
    cbar6 = plt.colorbar(im6, ax=ax6, label='Grain ID')
    ax6.set_xlabel('X (pixels)', fontsize=12, fontweight='bold')
    ax6.set_ylabel('Y (pixels)', fontsize=12, fontweight='bold')
    ax6.set_title(f'Grain Map ({len(np.unique(grain_grid))} grains, {grain_grid.shape[0]}×{grain_grid.shape[1]} grid)', 
                  fontsize=14, fontweight='bold')
    
    # 7. Statistics Summary (bottom right)
    ax7 = fig.add_subplot(gs[2, 2])
    ax7.axis('off')
    
    # Extract key correlation values
    center = max_dist
    corr_center = correlation_2d[center, center]
    corr_dx1 = correlation_2d[center, center+1] if center+1 < correlation_2d.shape[1] else 0
    corr_dy1 = correlation_2d[center+1, center] if center+1 < correlation_2d.shape[0] else 0
    corr_diag = correlation_2d[center+1, center+1] if center+1 < correlation_2d.shape[0] else 0
    
    # Create statistics text
    stats_summary = f"""
ANALYSIS SUMMARY

Grid Size: {grain_grid.shape[0]} × {grain_grid.shape[1]} px
Total Pixels: {grain_grid.shape[0] * grain_grid.shape[1]:,}

Grain Statistics:
  • Number of grains: {len(grain_sizes):,}
  • Min size: {min(grain_sizes.values())} px
  • Max size: {max(grain_sizes.values())} px
  • Mean size: {np.mean(list(grain_sizes.values())):.1f} px
  • Median size: {np.median(list(grain_sizes.values())):.1f} px
  • Std Dev: {np.std(list(grain_sizes.values())):.1f} px

2D Correlation (P = same grain):
  • Center (0,0): {corr_center:.3f}
  • ΔX=±1, ΔY=0: {corr_dx1:.3f}
  • ΔX=0, ΔY=±1: {corr_dy1:.3f}
  • ΔX=±1, ΔY=±1: {corr_diag:.3f}
  • Range: ±{max_dist} px
"""
    
    ax7.text(0.05, 0.95, stats_summary, transform=ax7.transAxes,
            verticalalignment='top', horizontalalignment='left',
            fontsize=11, family='monospace',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.7, pad=1))
    
    # Main title
    fig.suptitle('Comprehensive Grain Structure Analysis', 
                fontsize=18, fontweight='bold', y=0.995)
    
    if output_dir is None:
        output_dir = os.path.dirname(GRAIN_ORIENTATIONS_FILE)
    
    output_file = os.path.join(output_dir, "combined_analysis_highres.png")
    fig.savefig(output_file, dpi=150, bbox_inches='tight', facecolor='white')
    print(f"Combined high-resolution analysis figure saved to {output_file}")
    print(f"Image resolution: {int(fig.get_size_inches()[0]*150)}×{int(fig.get_size_inches()[1]*150)} pixels")
    
    plt.close(fig)


def main():
    """Main analysis pipeline."""
    print("=" * 60)
    print("GRAIN STRUCTURE ANALYSIS")
    print("=" * 60)
    
    # Load data
    grain_grid, orientations_df = load_data()
    
    # Set output directory
    output_dir = os.path.dirname(GRAIN_ORIENTATIONS_FILE)
    
    # Create pole figure with per-pixel orientation distribution
    create_pole_figure(grain_grid, orientations_df, output_dir)
    
    # Calculate grain sizes
    grain_sizes = calculate_grain_sizes(grain_grid)
    
    # Plot grain size distributions
    pixel_grain_sizes = plot_grain_size_distribution_per_pixel(grain_grid, grain_sizes, output_dir)
    
    # Create additional visualizations
    plot_grain_size_map(grain_grid, grain_sizes, output_dir)
    plot_grain_id_map(grain_grid, output_dir)
    
    # Calculate 2D two-point correlation
    correlation_2d, max_dist = calculate_two_point_correlation_2d(grain_grid, max_distance=50)
    
    # Create combined high-resolution figure with all analyses
    create_combined_analysis_figure(grain_grid, orientations_df, grain_sizes, correlation_2d, max_dist, output_dir)
    
    print("\n" + "=" * 60)
    print("ANALYSIS COMPLETE")
    print("=" * 60)
    print(f"Output files saved to: {output_dir}")


if __name__ == "__main__":
    main()
