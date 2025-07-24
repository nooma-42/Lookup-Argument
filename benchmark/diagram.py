#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Lookup Argument Benchmark Visualization Tool

This script generates performance plots from benchmark data, showing proving time
vs lookup table size for different N:n ratios and proof systems.

Usage:
    python diagram.py [options]
    python diagram.py --protocols Lasso,LogupGKR --output-dir plots/
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.style as style
import argparse
import sys
import os
from pathlib import Path
from typing import List, Dict, Set
import numpy as np

class BenchmarkVisualizer:
    """Main class for processing benchmark data and generating visualizations."""
    
    def __init__(self, csv_file: str = "benchmark_results.csv"):
        """Initialize with CSV file path."""
        self.csv_file = csv_file
        self.data = None
        self.available_systems = set()
        self.available_ratios = set()
        
    def load_data(self) -> bool:
        """Load and validate CSV data."""
        try:
            self.data = pd.read_csv(self.csv_file)
            
            # Validate required columns
            required_cols = ['System', 'K', 'N_to_n_Ratio', 'ProveTime_ms', 'SetupTime_ms', 'VerifyTime_ms', 'TotalTime_ms', 'ProofSize_bytes']
            missing_cols = [col for col in required_cols if col not in self.data.columns]
            if missing_cols:
                print(f"Error: Missing required columns: {missing_cols}")
                return False
            
            # Convert to appropriate types
            self.data['K'] = pd.to_numeric(self.data['K'], errors='coerce')
            self.data['N_to_n_Ratio'] = pd.to_numeric(self.data['N_to_n_Ratio'], errors='coerce')
            self.data['ProveTime_ms'] = pd.to_numeric(self.data['ProveTime_ms'], errors='coerce')
            self.data['SetupTime_ms'] = pd.to_numeric(self.data['SetupTime_ms'], errors='coerce')
            self.data['VerifyTime_ms'] = pd.to_numeric(self.data['VerifyTime_ms'], errors='coerce')
            self.data['TotalTime_ms'] = pd.to_numeric(self.data['TotalTime_ms'], errors='coerce')
            self.data['ProofSize_bytes'] = pd.to_numeric(self.data['ProofSize_bytes'], errors='coerce')
            
            # Remove rows with invalid data
            initial_rows = len(self.data)
            self.data = self.data.dropna(subset=['K', 'N_to_n_Ratio', 'ProveTime_ms', 'SetupTime_ms', 'VerifyTime_ms', 'TotalTime_ms', 'ProofSize_bytes'])
            removed_rows = initial_rows - len(self.data)
            if removed_rows > 0:
                print(f"Warning: Removed {removed_rows} rows with invalid data")
            
            # Extract available systems and ratios
            self.available_systems = set(self.data['System'].unique())
            self.available_ratios = set(self.data['N_to_n_Ratio'].unique())
            
            print(f"Loaded {len(self.data)} data points")
            print(f"Available systems: {sorted(self.available_systems)}")
            print(f"Available ratios: {sorted(self.available_ratios)}")
            
            return True
            
        except FileNotFoundError:
            print(f"Error: File '{self.csv_file}' not found")
            return False
        except Exception as e:
            print(f"Error loading data: {e}")
            return False
    
    def filter_systems(self, systems: List[str]) -> List[str]:
        """Filter and validate requested systems."""
        if not systems:
            return sorted(self.available_systems)
        
        valid_systems = []
        for system in systems:
            if system in self.available_systems:
                valid_systems.append(system)
            else:
                print(f"Warning: System '{system}' not found in data. Available: {sorted(self.available_systems)}")
        
        return valid_systems
    
    def get_data_for_ratio(self, ratio: float, systems: List[str]) -> Dict[str, pd.DataFrame]:
        """Get data for a specific ratio, grouped by system."""
        ratio_data = self.data[self.data['N_to_n_Ratio'] == ratio]
        system_data = {}
        
        for system in systems:
            system_subset = ratio_data[ratio_data['System'] == system]
            if not system_subset.empty:
                # Sort by K for proper line plotting
                system_subset = system_subset.sort_values('K')
                system_data[system] = system_subset
        
        return system_data
    
    def create_plot(self, ratio: float, system_data: Dict[str, pd.DataFrame], 
                   metric: str = 'ProveTime_ms', log_scale: bool = True, output_dir: str = ".") -> str:
        """Create and save a plot for a specific ratio and metric."""
        
        # Set up the plot with professional styling
        plt.style.use('default')
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Color palette for different systems
        colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', 
                 '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
        
        lines_plotted = 0
        
        for i, (system, data) in enumerate(system_data.items()):
            if data.empty or metric not in data.columns:
                continue
                
            color = colors[i % len(colors)]
            
            # Plot the line
            ax.plot(data['K'], data[metric], 
                   marker='o', linewidth=2.5, markersize=6, 
                   label=system, color=color, alpha=0.8)
            
            lines_plotted += 1
        
        if lines_plotted == 0:
            print(f"Warning: No data to plot for ratio {ratio} and metric {metric}")
            plt.close(fig)
            return ""
        
        # Metric-specific labels and titles
        metric_labels = {
            'ProveTime_ms': ('Proving Time (ms)', 'Proving Time vs Table Size'),
            'SetupTime_ms': ('Setup Time (ms)', 'Setup Time vs Table Size'),
            'VerifyTime_ms': ('Verify Time (ms)', 'Verify Time vs Table Size'),
            'TotalTime_ms': ('Total Time (ms)', 'Total Time vs Table Size'),
            'ProofSize_bytes': ('Proof Size (bytes)', 'Proof Size vs Table Size')
        }
        
        ylabel, title_prefix = metric_labels.get(metric, (metric, metric.replace('_', ' ').title()))
        
        # Customize the plot
        ax.set_xlabel('Lookup Table Size (K)', fontsize=14, fontweight='bold')
        ax.set_ylabel(ylabel, fontsize=14, fontweight='bold')
        ax.set_title(f'{title_prefix} (N:n Ratio = {ratio})', 
                    fontsize=16, fontweight='bold', pad=20)
        
        # Set logarithmic scale if requested
        if log_scale:
            ax.set_yscale('log')
            ax.set_ylabel(f'{ylabel} - Log Scale', fontsize=14, fontweight='bold')
        
        # Add grid
        ax.grid(True, alpha=0.3, linestyle='-', linewidth=0.5)
        
        # Customize legend
        ax.legend(frameon=True, fancybox=True, shadow=True, 
                 loc='best', fontsize=12)
        
        # Set integer ticks for K axis
        k_values = []
        for data in system_data.values():
            k_values.extend(data['K'].tolist())
        if k_values:
            k_range = np.arange(int(min(k_values)), int(max(k_values)) + 1)
            ax.set_xticks(k_range)
        
        # Improve layout
        plt.tight_layout()
        
        # Save the plot
        metric_name = metric.replace('_', '').lower()
        filename = f"{metric_name}_ratio_{ratio}.png"
        filepath = os.path.join(output_dir, filename)
        plt.savefig(filepath, dpi=300, bbox_inches='tight', 
                   facecolor='white', edgecolor='none')
        plt.close(fig)
        
        return filepath
    
    def generate_plots(self, systems: List[str] = None, output_dir: str = ".", 
                      log_scale: bool = True, ratio_protocols: Dict[float, List[str]] = None,
                      metrics: List[str] = None, single_ratio: float = None) -> List[str]:
        """Generate plots for specified metrics and ratios with optional per-ratio protocol specification."""
        if self.data is None:
            print("Error: No data loaded")
            return []
        
        # Default metrics if none specified
        if metrics is None:
            metrics = ['ProveTime_ms']
        
        # Validate metrics
        available_metrics = ['ProveTime_ms', 'SetupTime_ms', 'VerifyTime_ms', 'TotalTime_ms', 'ProofSize_bytes']
        valid_metrics = [m for m in metrics if m in available_metrics]
        if not valid_metrics:
            print(f"Error: No valid metrics specified. Available: {available_metrics}")
            return []
        
        # Ensure output directory exists
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        generated_files = []
        
        # Determine which ratios to process
        ratios_to_process = [single_ratio] if single_ratio is not None else sorted(self.available_ratios)
        
        for ratio in ratios_to_process:
            if ratio not in self.available_ratios:
                print(f"Warning: Ratio {ratio} not found in data")
                continue
                
            print(f"Generating plots for ratio {ratio}...")
            
            # Determine which systems to use for this ratio
            if ratio_protocols and ratio in ratio_protocols:
                # Use ratio-specific protocols
                ratio_systems = ratio_protocols[ratio]
                print(f"  Using ratio-specific protocols: {ratio_systems}")
            elif systems:
                # Use global systems
                ratio_systems = systems
            else:
                # Use all available systems
                ratio_systems = []
            
            # Filter and validate systems for this ratio
            valid_systems = self.filter_systems(ratio_systems)
            if not valid_systems:
                print(f"Warning: No valid systems specified for ratio {ratio}")
                continue
            
            system_data = self.get_data_for_ratio(ratio, valid_systems)
            
            if not system_data:
                print(f"Warning: No data found for ratio {ratio}")
                continue
            
            # Generate plots for each metric
            for metric in valid_metrics:
                filepath = self.create_plot(ratio, system_data, metric, log_scale, output_dir)
                if filepath:
                    generated_files.append(filepath)
                    print(f"  -> Saved: {filepath}")
        
        return generated_files
    
    def plot_single_ratio_metrics(self, ratio: float, metrics: List[str], systems: List[str] = None, 
                                 output_dir: str = ".", log_scale: bool = True) -> List[str]:
        """Plot specific metrics for a single ratio."""
        return self.generate_plots(systems=systems, output_dir=output_dir, log_scale=log_scale, 
                                 metrics=metrics, single_ratio=ratio)


def parse_ratio_protocols(ratio_protocols_str: str) -> Dict[float, List[str]]:
    """Parse ratio-specific protocol specification string.
    
    Format: "ratio1:protocol1,protocol2;ratio2:protocol3,protocol4"
    Example: "2:Lasso,LogupGKR;4:CQ,Baloo;8:Plookup"
    """
    ratio_protocols = {}
    
    if not ratio_protocols_str:
        return ratio_protocols
    
    try:
        # Split by semicolon to get each ratio specification
        ratio_specs = ratio_protocols_str.split(';')
        
        for spec in ratio_specs:
            if not spec.strip():
                continue
                
            # Split by colon to get ratio and protocols
            if ':' not in spec:
                print(f"Warning: Invalid ratio specification '{spec}' - missing ':'")
                continue
                
            ratio_str, protocols_str = spec.split(':', 1)
            
            try:
                ratio = float(ratio_str.strip())
            except ValueError:
                print(f"Warning: Invalid ratio '{ratio_str}' - must be a number")
                continue
                
            # Split protocols by comma
            protocols = [p.strip() for p in protocols_str.split(',') if p.strip()]
            
            if protocols:
                ratio_protocols[ratio] = protocols
                
    except Exception as e:
        print(f"Error parsing ratio-protocols: {e}")
        
    return ratio_protocols


def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Generate performance plots from benchmark data for various metrics",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python diagram.py
  python diagram.py --protocols Lasso,LogupGKR
  python diagram.py --protocols CQ,Baloo --output-dir plots/ --linear-scale
  python diagram.py --ratio-protocols "2:Lasso,LogupGKR;4:CQ,Baloo;8:Plookup"
  python diagram.py --ratio-protocols "2:Lasso;4:CQ,Caulk;16:LogupGKR,Plookup"
  python diagram.py --metrics SetupTime_ms,TotalTime_ms --protocols Lasso,LogupGKR
  python diagram.py --metrics VerifyTime_ms,ProofSize_bytes --single-ratio 4
        """
    )
    
    parser.add_argument(
        '--input', '-i',
        default='benchmark_results.csv',
        help='Input CSV file path (default: benchmark_results.csv)'
    )
    
    parser.add_argument(
        '--protocols', '-p',
        default='',
        help='Comma-separated list of protocols to plot (default: all)'
    )
    
    parser.add_argument(
        '--ratio-protocols', '-r',
        default='',
        help='Per-ratio protocol specification: "ratio1:protocol1,protocol2;ratio2:protocol3"'
    )
    
    parser.add_argument(
        '--output-dir', '-o',
        default='.',
        help='Output directory for plots (default: current directory)'
    )
    
    parser.add_argument(
        '--linear-scale',
        action='store_true',
        help='Use linear scale for Y-axis (default: logarithmic)'
    )
    
    parser.add_argument(
        '--metrics', '-m',
        default='ProveTime_ms',
        help='Comma-separated list of metrics to plot (ProveTime_ms,SetupTime_ms,VerifyTime_ms,TotalTime_ms,ProofSize_bytes)'
    )
    
    parser.add_argument(
        '--single-ratio',
        type=float,
        help='Plot only for a specific ratio (useful for verify time and proof size)'
    )
    
    parser.add_argument(
        '--list-systems',
        action='store_true',
        help='List available systems in the data and exit'
    )
    
    return parser.parse_args()


def main():
    """Main function."""
    args = parse_arguments()
    
    # Initialize visualizer
    visualizer = BenchmarkVisualizer(args.input)
    
    # Load data
    if not visualizer.load_data():
        sys.exit(1)
    
    # List systems if requested
    if args.list_systems:
        print(f"Available systems: {sorted(visualizer.available_systems)}")
        print(f"Available ratios: {sorted(visualizer.available_ratios)}")
        return
    
    # Parse systems
    if args.protocols:
        systems = [s.strip() for s in args.protocols.split(',') if s.strip()]
    else:
        systems = []
    
    # Parse metrics
    metrics = [m.strip() for m in args.metrics.split(',') if m.strip()]
    
    # Parse ratio-specific protocols
    ratio_protocols = None
    if args.ratio_protocols:
        ratio_protocols = parse_ratio_protocols(args.ratio_protocols)
        if ratio_protocols:
            print(f"Using ratio-specific protocols: {ratio_protocols}")
        else:
            print("Warning: No valid ratio-protocols specified, using global protocols")
    
    # Check for conflicting arguments
    if args.protocols and args.ratio_protocols:
        print("Warning: Both --protocols and --ratio-protocols specified. --ratio-protocols takes precedence.")
    
    # Generate plots
    log_scale = not args.linear_scale
    generated_files = visualizer.generate_plots(systems, args.output_dir, log_scale, ratio_protocols, 
                                              metrics, args.single_ratio)
    
    # Summary
    if generated_files:
        print(f"\nSuccessfully generated {len(generated_files)} plots:")
        for filepath in generated_files:
            print(f"  - {filepath}")
    else:
        print("No plots were generated")
        sys.exit(1)


if __name__ == "__main__":
    main()