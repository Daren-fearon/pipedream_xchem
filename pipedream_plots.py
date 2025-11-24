"""
Pipedream Plotting Functions Module

Contains all plotting and visualization functions for Pipedream results.
Extracted from main collation script for better modularity.
"""

import os
import logging
import argparse
from typing import List, Dict, Any, Tuple
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# Try to import plotly for interactive plots
try:
    import plotly.graph_objects as go
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    logging.debug("plotly not available, interactive plots will not be generated")

from pipedream_thresholds import THRESHOLDS


def generate_pca_plots(results: List[Dict[str, Any]], output_dir: str, base_name: str) -> Tuple[str, str]:
    """
    Generate interactive PCA plots for crystallographic quality metrics.
    
    Args:
        results: List of result dictionaries
        output_dir: Directory to save plots
        base_name: Base name for output files
        
    Returns:
        Tuple of (ligand_metrics_plot_path, geometry_metrics_plot_path)
    """
    if not HAS_PLOTLY:
        logging.warning("plotly not available, skipping PCA plots")
        return 'NA', 'NA'
    
    df = pd.DataFrame(results)
    
    # Plot 1: Ligand RSR vs Ligand CC
    ligand_plot_path = os.path.join(output_dir, f'{base_name}_pca_ligand_metrics.html')
    
    # Extract and clean ligand metrics
    ligand_data = df[['Crystal Name', 'Ligand ID', 'Ligand RSR', 'Ligand CC']].copy()
    ligand_data[['Ligand RSR', 'Ligand CC']] = ligand_data[['Ligand RSR', 'Ligand CC']].apply(pd.to_numeric, errors='coerce')
    ligand_data = ligand_data.dropna()
    
    if len(ligand_data) >= 2:
        # Perform PCA
        scaler = StandardScaler()
        ligand_scaled = scaler.fit_transform(ligand_data[['Ligand RSR', 'Ligand CC']])
        pca_ligand = PCA(n_components=2)
        ligand_pca = pca_ligand.fit_transform(ligand_scaled)
        
        # Create hover text with Ligand ID
        hover_text_scatter = [f"<b>{name}</b><br>Ligand ID: {lig_id}<br>Ligand CC: {cc:.3f}<br>Ligand RSR: {rsr:.3f}<br>PC1: {pc1:.2f}"
                             for name, lig_id, cc, rsr, pc1 in zip(ligand_data['Crystal Name'],
                                                           ligand_data['Ligand ID'],
                                                           ligand_data['Ligand CC'],
                                                           ligand_data['Ligand RSR'],
                                                           ligand_pca[:, 0])]
        
        hover_text_pca = [f"<b>{name}</b><br>Ligand ID: {lig_id}<br>PC1: {pc1:.2f}<br>PC2: {pc2:.2f}"
                         for name, lig_id, pc1, pc2 in zip(ligand_data['Crystal Name'],
                                                   ligand_data['Ligand ID'],
                                                   ligand_pca[:, 0], ligand_pca[:, 1])]
        
        # Create subplots
        from plotly.subplots import make_subplots
        fig = make_subplots(rows=1, cols=2, 
                           horizontal_spacing=0.20,
                           subplot_titles=('Ligand Quality Metrics (colored by PC1)', 
                                         f'PCA of Ligand Metrics'))
        
        # Left plot: Scatter with thresholds and quality zones
        fig.add_trace(
            go.Scatter(x=ligand_data['Ligand RSR'], y=ligand_data['Ligand CC'],
                      mode='markers',
                      marker=dict(size=10, color=ligand_pca[:, 0], 
                                colorscale='Viridis', showscale=True,
                                colorbar=dict(title="PC1", x=0.46, len=0.75),
                                line=dict(color='white', width=1)),
                      text=hover_text_scatter,
                      hovertemplate='%{text}<extra></extra>',
                      name='Datasets'),
            row=1, col=1
        )
        
        # Add quality zones with consistent colors
        # Poor quality zone (high RSR, low CC) - light red
        fig.add_shape(type="rect", x0=THRESHOLDS.ligand_rsr_max, y0=0, 
                     x1=1.0, y1=THRESHOLDS.ligand_cc,
                     fillcolor="#FFE5E5", opacity=0.3, layer="below", line_width=0,
                     row=1, col=1)
        
        # Good quality zone (low RSR, high CC) - light green
        fig.add_shape(type="rect", x0=0, y0=THRESHOLDS.ligand_cc, 
                     x1=THRESHOLDS.ligand_rsr_max, y1=1.0,
                     fillcolor="#B8E6B8", opacity=0.3, layer="below", line_width=0,
                     row=1, col=1)
        
        # Add threshold lines
        fig.add_hline(y=THRESHOLDS.ligand_cc, line_dash="dash", line_color="darkgreen", 
                     opacity=0.6, line_width=2.5, row=1, col=1)
        fig.add_vline(x=THRESHOLDS.ligand_rsr_max, line_dash="dash", line_color="red", 
                     opacity=0.6, line_width=2.5, row=1, col=1)
        
        # Right plot: PCA biplot
        fig.add_trace(
            go.Scatter(x=ligand_pca[:, 0], y=ligand_pca[:, 1],
                      mode='markers',
                      marker=dict(size=10, color='steelblue',
                                line=dict(color='white', width=1)),
                      text=hover_text_pca,
                      hovertemplate='%{text}<extra></extra>',
                      name='PCA'),
            row=1, col=2
        )
        
        # Add loading vectors
        loadings = pca_ligand.components_.T * np.sqrt(pca_ligand.explained_variance_)
        for i, (feature, loading) in enumerate(zip(['Ligand RSR', 'Ligand CC'], loadings)):
            fig.add_annotation(x=loading[0]*3, y=loading[1]*3,
                             ax=0, ay=0, xref=f'x2', yref=f'y2',
                             axref=f'x2', ayref=f'y2',
                             showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2,
                             arrowcolor='red', opacity=0.7)
            fig.add_annotation(x=loading[0]*3.5, y=loading[1]*3.5,
                             text=feature, showarrow=False,
                             font=dict(size=11, color='red'),
                             xref=f'x2', yref=f'y2')
        
        # Update axes
        fig.update_xaxes(title_text="Ligand RSR", gridcolor='lightgray', 
                        zeroline=False, row=1, col=1)
        fig.update_yaxes(title_text="Ligand CC", gridcolor='lightgray', 
                        zeroline=False, row=1, col=1)
        fig.update_xaxes(title_text=f"PC1 ({pca_ligand.explained_variance_ratio_[0]:.1%} variance)", 
                        gridcolor='lightgray', zeroline=True, zerolinecolor='gray', row=1, col=2)
        fig.update_yaxes(title_text=f"PC2 ({pca_ligand.explained_variance_ratio_[1]:.1%} variance)", 
                        gridcolor='lightgray', zeroline=True, zerolinecolor='gray', row=1, col=2)
        
        # Update layout with consistent styling
        fig.update_layout(
            title=dict(
                text="<b>Ligand Quality Metrics - PCA Analysis</b>",
                x=0.5,
                xanchor='center',
                font=dict(size=19)
            ),
            showlegend=False,
            width=1400,
            height=600,
            hovermode='closest',
            plot_bgcolor='white',
            margin=dict(l=60, r=60, t=80, b=60)
        )
        
        # Add modebar config for consistent export
        config = {
            'displayModeBar': True,
            'displaylogo': False,
            'modeBarButtonsToRemove': ['select2d', 'lasso2d'],
            'toImageButtonOptions': {
                'format': 'png',
                'filename': 'ligand_pca_plot',
                'height': 950,
                'width': 1100,
                'scale': 2
            }
        }
        
        import plotly.io as pio
        html_string = pio.to_html(fig, 
                                  config=config,
                                  include_plotlyjs='cdn',
                                  full_html=True)
        html_string = html_string.replace('<head><meta charset="utf-8" /></head>', 
                                         f'<head><meta charset="utf-8" /><title>Ligand Quality PCA Analysis</title></head>')
        with open(ligand_plot_path, 'w', encoding='utf-8') as f:
            f.write(html_string)
        
        logging.info(f"Generated interactive ligand metrics PCA plot: {ligand_plot_path}")
    else:
        ligand_plot_path = 'NA'
        logging.warning("Insufficient data for ligand metrics PCA plot")
    
    # Plot 2: Mogul Z Bond vs Mogul Z Angle
    geometry_plot_path = os.path.join(output_dir, f'{base_name}_pca_geometry_metrics.html')
    
    # Extract and clean geometry metrics
    geometry_data = df[['Crystal Name', 'Ligand ID', 'Mogul Z Bond', 'Mogul Z Angle']].copy()
    geometry_data[['Mogul Z Bond', 'Mogul Z Angle']] = geometry_data[['Mogul Z Bond', 'Mogul Z Angle']].apply(pd.to_numeric, errors='coerce')
    geometry_data = geometry_data.dropna()
    
    if len(geometry_data) >= 2:
        # Perform PCA
        scaler = StandardScaler()
        geometry_scaled = scaler.fit_transform(geometry_data[['Mogul Z Bond', 'Mogul Z Angle']])
        pca_geometry = PCA(n_components=2)
        geometry_pca = pca_geometry.fit_transform(geometry_scaled)
        
        # Create hover text with Ligand ID
        hover_text_scatter = [f"<b>{name}</b><br>Ligand ID: {lig_id}<br>Mogul Z Bond: {bond:.2f}<br>Mogul Z Angle: {angle:.2f}<br>PC1: {pc1:.2f}"
                             for name, lig_id, bond, angle, pc1 in zip(geometry_data['Crystal Name'],
                                                               geometry_data['Ligand ID'],
                                                               geometry_data['Mogul Z Bond'],
                                                               geometry_data['Mogul Z Angle'],
                                                               geometry_pca[:, 0])]
        
        hover_text_pca = [f"<b>{name}</b><br>Ligand ID: {lig_id}<br>PC1: {pc1:.2f}<br>PC2: {pc2:.2f}"
                         for name, lig_id, pc1, pc2 in zip(geometry_data['Crystal Name'],
                                                   geometry_data['Ligand ID'],
                                                   geometry_pca[:, 0], geometry_pca[:, 1])]
        
        # Create subplots
        from plotly.subplots import make_subplots
        fig = make_subplots(rows=1, cols=2,
                           horizontal_spacing=0.20,
                           subplot_titles=('Geometry Quality Metrics (colored by PC1)', 
                                         f'PCA of Geometry Metrics'))
        
        # Left plot: Scatter with quality zones
        fig.add_trace(
            go.Scatter(x=geometry_data['Mogul Z Bond'], y=geometry_data['Mogul Z Angle'],
                      mode='markers',
                      marker=dict(size=10, color=geometry_pca[:, 0], 
                                colorscale='Plasma', showscale=True,
                                colorbar=dict(title="PC1", x=0.46, len=0.75),
                                line=dict(color='white', width=1)),
                      text=hover_text_scatter,
                      hovertemplate='%{text}<extra></extra>',
                      name='Datasets'),
            row=1, col=1
        )
        
        # Add quality zones matching the 2D quality plot
        # Poor zone (red) - both metrics > 4.0
        fig.add_shape(type="rect", x0=THRESHOLDS.mogul_z_poor, y0=THRESHOLDS.mogul_z_poor, 
                     x1=10, y1=10,
                     fillcolor="#E6B8E6", opacity=0.25, layer="below", line_width=0,
                     row=1, col=1)
        
        # OK zones (orange) - one metric between 2.5-4.0
        fig.add_shape(type="rect", x0=THRESHOLDS.mogul_z_ok, y0=0, 
                     x1=THRESHOLDS.mogul_z_poor, y1=THRESHOLDS.mogul_z_ok,
                     fillcolor="#FFE5CC", opacity=0.25, layer="below", line_width=0,
                     row=1, col=1)
        fig.add_shape(type="rect", x0=0, y0=THRESHOLDS.mogul_z_ok, 
                     x1=THRESHOLDS.mogul_z_ok, y1=THRESHOLDS.mogul_z_poor,
                     fillcolor="#FFE5CC", opacity=0.25, layer="below", line_width=0,
                     row=1, col=1)
        fig.add_shape(type="rect", x0=THRESHOLDS.mogul_z_ok, y0=THRESHOLDS.mogul_z_ok, 
                     x1=THRESHOLDS.mogul_z_poor, y1=THRESHOLDS.mogul_z_poor,
                     fillcolor="#FFE5CC", opacity=0.25, layer="below", line_width=0,
                     row=1, col=1)
        
        # Good zones (yellow) - one metric between 1.5-2.5
        fig.add_shape(type="rect", x0=THRESHOLDS.mogul_z_good, y0=0, 
                     x1=THRESHOLDS.mogul_z_ok, y1=THRESHOLDS.mogul_z_good,
                     fillcolor="#FFFFCC", opacity=0.25, layer="below", line_width=0,
                     row=1, col=1)
        fig.add_shape(type="rect", x0=0, y0=THRESHOLDS.mogul_z_good, 
                     x1=THRESHOLDS.mogul_z_good, y1=THRESHOLDS.mogul_z_ok,
                     fillcolor="#FFFFCC", opacity=0.25, layer="below", line_width=0,
                     row=1, col=1)
        fig.add_shape(type="rect", x0=THRESHOLDS.mogul_z_good, y0=THRESHOLDS.mogul_z_good, 
                     x1=THRESHOLDS.mogul_z_ok, y1=THRESHOLDS.mogul_z_ok,
                     fillcolor="#FFFFCC", opacity=0.25, layer="below", line_width=0,
                     row=1, col=1)
        
        # OK zones extending to good boundary
        fig.add_shape(type="rect", x0=THRESHOLDS.mogul_z_ok, y0=0, 
                     x1=THRESHOLDS.mogul_z_poor, y1=THRESHOLDS.mogul_z_poor,
                     fillcolor="#FFE5CC", opacity=0.25, layer="below", line_width=0,
                     row=1, col=1)
        fig.add_shape(type="rect", x0=0, y0=THRESHOLDS.mogul_z_ok, 
                     x1=THRESHOLDS.mogul_z_poor, y1=THRESHOLDS.mogul_z_poor,
                     fillcolor="#FFE5CC", opacity=0.25, layer="below", line_width=0,
                     row=1, col=1)
        
        # Excellent zone (green) - both < 1.5
        fig.add_shape(type="rect", x0=0, y0=0, 
                     x1=THRESHOLDS.mogul_z_good, y1=THRESHOLDS.mogul_z_good,
                     fillcolor="#B8E6B8", opacity=0.25, layer="below", line_width=0,
                     row=1, col=1)
        
        # Add reference lines
        fig.add_hline(y=THRESHOLDS.mogul_z_poor, line_dash="dash", line_color="purple", 
                     opacity=0.6, line_width=2.5, row=1, col=1)
        fig.add_vline(x=THRESHOLDS.mogul_z_poor, line_dash="dash", line_color="purple", 
                     opacity=0.6, line_width=2.5, row=1, col=1)
        
        fig.add_hline(y=THRESHOLDS.mogul_z_ok, line_dash="dash", line_color="orange", 
                     opacity=0.6, line_width=2.5, row=1, col=1)
        fig.add_vline(x=THRESHOLDS.mogul_z_ok, line_dash="dash", line_color="orange", 
                     opacity=0.6, line_width=2.5, row=1, col=1)
        
        fig.add_hline(y=THRESHOLDS.mogul_z_good, line_dash="dash", line_color="darkgreen", 
                     opacity=0.6, line_width=2.5, row=1, col=1)
        fig.add_vline(x=THRESHOLDS.mogul_z_good, line_dash="dash", line_color="darkgreen", 
                     opacity=0.6, line_width=2.5, row=1, col=1)
        
        # Right plot: PCA biplot
        fig.add_trace(
            go.Scatter(x=geometry_pca[:, 0], y=geometry_pca[:, 1],
                      mode='markers',
                      marker=dict(size=10, color='coral',
                                line=dict(color='white', width=1)),
                      text=hover_text_pca,
                      hovertemplate='%{text}<extra></extra>',
                      name='PCA'),
            row=1, col=2
        )
        
        # Add loading vectors
        loadings = pca_geometry.components_.T * np.sqrt(pca_geometry.explained_variance_)
        for i, (feature, loading) in enumerate(zip(['Mogul Z Bond', 'Mogul Z Angle'], loadings)):
            fig.add_annotation(x=loading[0]*3, y=loading[1]*3,
                             ax=0, ay=0, xref=f'x2', yref=f'y2',
                             axref=f'x2', ayref=f'y2',
                             showarrow=True, arrowhead=2, arrowsize=1, arrowwidth=2,
                             arrowcolor='red', opacity=0.7)
            fig.add_annotation(x=loading[0]*3.5, y=loading[1]*3.5,
                             text=feature, showarrow=False,
                             font=dict(size=11, color='red'),
                             xref=f'x2', yref=f'y2')
        
        # Update axes
        fig.update_xaxes(title_text="Mogul Z Bond", gridcolor='lightgray', 
                        zeroline=False, row=1, col=1)
        fig.update_yaxes(title_text="Mogul Z Angle", gridcolor='lightgray', 
                        zeroline=False, row=1, col=1)
        fig.update_xaxes(title_text=f"PC1 ({pca_geometry.explained_variance_ratio_[0]:.1%} variance)", 
                        gridcolor='lightgray', zeroline=True, zerolinecolor='gray', row=1, col=2)
        fig.update_yaxes(title_text=f"PC2 ({pca_geometry.explained_variance_ratio_[1]:.1%} variance)", 
                        gridcolor='lightgray', zeroline=True, zerolinecolor='gray', row=1, col=2)
        
        # Update layout with consistent styling
        fig.update_layout(
            title=dict(
                text="<b>Geometry Quality Metrics - PCA Analysis</b>",
                x=0.5,
                xanchor='center',
                font=dict(size=19)
            ),
            showlegend=False,
            width=1400,
            height=600,
            hovermode='closest',
            plot_bgcolor='white',
            margin=dict(l=60, r=60, t=80, b=60)
        )
        
        # Add modebar config for consistent export
        config = {
            'displayModeBar': True,
            'displaylogo': False,
            'modeBarButtonsToRemove': ['select2d', 'lasso2d'],
            'toImageButtonOptions': {
                'format': 'png',
                'filename': 'geometry_pca_plot',
                'height': 950,
                'width': 1100,
                'scale': 2
            }
        }
        
        # Create custom HTML with proper title
        import plotly.io as pio
        html_string = pio.to_html(fig,
                                  config=config,
                                  include_plotlyjs='cdn',
                                  full_html=True)
        # Replace the default title
        html_string = html_string.replace('<head><meta charset="utf-8" /></head>',
                                         f'<head><meta charset="utf-8" /><title>Geometry Quality PCA Analysis</title></head>')
        with open(geometry_plot_path, 'w', encoding='utf-8') as f:
            f.write(html_string)
        
        logging.info(f"Generated interactive geometry metrics PCA plot: {geometry_plot_path}")
    else:
        geometry_plot_path = 'NA'
        logging.warning("Insufficient data for geometry metrics PCA plot")
    
    return ligand_plot_path, geometry_plot_path


def generate_ligand_quality_plot(results: List[Dict[str, Any]], output_dir: str, base_name: str) -> str:
    """
    Generate an interactive 2D ligand quality plot combining fit and geometry metrics.
    
    Creates an interactive scatter plot with:
    - X-axis: Combined fit quality (Ligand CC and RSR)
    - Y-axis: Combined geometry quality (Mogul Z Bond and Angle)
    - Color: Composite quality score
    - Size: Ligand Clashes (higher = larger, indicating problems)
    - Hover: Shows Crystal Name and all quality metrics
    
    Args:
        results: List of result dictionaries
        output_dir: Directory to save plot
        base_name: Base name for output file
        
    Returns:
        Path to generated HTML plot or 'NA' if insufficient data
    """
    if not HAS_PLOTLY:
        logging.warning("plotly not available, skipping interactive quality plot")
        return 'NA'
    
    df = pd.DataFrame(results)
    
    # Extract relevant metrics
    metrics = df[['Crystal Name', 'Ligand ID', 'Ligand CC', 'Ligand RSR', 'Mogul Z Bond', 'Mogul Z Angle', 
                  'Ligand Clashes', 'R', 'Rfree', 'High Resolution (Å)', 'Composite Quality Score']].copy()
    
    # Convert 'NA' strings to numeric (pd.to_numeric handles 'NA' automatically)
    # Skip Crystal Name and Ligand ID (text columns)
    for col in metrics.columns[2:]:
        metrics[col] = pd.to_numeric(metrics[col], errors='coerce')
    
    # Remove rows with missing critical data
    metrics = metrics.dropna(subset=['Ligand CC', 'Ligand RSR', 'Mogul Z Bond', 'Mogul Z Angle'])
    
    if len(metrics) < 2:
        logging.warning("Insufficient data for ligand quality plot (need at least 2 complete datasets)")
        return 'NA'
    
    # Calculate composite scores (0-100 scale, higher is better)
    
    # Ligand CC Score
    def score_ligand_cc(cc):
        if cc >= 0.95:  # Very good fit
            return 100
        elif cc >= 0.90:  # Good fit
            return 90 + ((cc - 0.90) / 0.05) * 10
        elif cc >= 0.80:  # OK fit
            return 50 + ((cc - 0.80) / 0.10) * 40
        else:  # Poor fit (< 0.80)
            return np.maximum(((cc - 0.60) / 0.20) * 50, 0)
    
    cc_score = metrics['Ligand CC'].apply(score_ligand_cc)
    
    # Ligand RSR Score
    def score_rsr(rsr):
        if rsr < 0.2:  # Very good
            return 100 - (rsr / 0.2) * 10
        elif rsr < 0.3:  # Good
            return 90 - ((rsr - 0.2) / 0.1) * 20
        elif rsr < 0.4:  # OK
            return 70 - ((rsr - 0.3) / 0.1) * 30
        else:  # Poor (>= 0.4)
            return np.maximum(40 - ((rsr - 0.4) / 0.2) * 40, 0)
    
    rsr_score = metrics['Ligand RSR'].apply(score_rsr)
    metrics['Fit Quality'] = (cc_score + rsr_score) / 2
    
    # Geometry Quality Score
    def score_mogul_z(z_val):
        if z_val < THRESHOLDS.mogul_z_good:  # < 1.5: good
            return 100 - (z_val / THRESHOLDS.mogul_z_good) * 20
        elif z_val < THRESHOLDS.mogul_z_ok:  # 1.5-2.5: ok
            return 80 - ((z_val - THRESHOLDS.mogul_z_good) / 
                          (THRESHOLDS.mogul_z_ok - THRESHOLDS.mogul_z_good) * 30)
        elif z_val < THRESHOLDS.mogul_z_poor:  # 2.5-4.0: poor
            return 50 - ((z_val - THRESHOLDS.mogul_z_ok) / 
                         (THRESHOLDS.mogul_z_poor - THRESHOLDS.mogul_z_ok) * 30)
        else:  # > 4.0: bad
            return np.maximum(20 - ((z_val - THRESHOLDS.mogul_z_poor) / 2.0 * 20), 0)
    
    mogul_bond_score = metrics['Mogul Z Bond'].apply(score_mogul_z)
    mogul_angle_score = metrics['Mogul Z Angle'].apply(score_mogul_z)
    metrics['Geometry Quality'] = (mogul_bond_score + mogul_angle_score) / 2
    
    # Composite Quality Score (using pre-calculated values from dataframe)
    metrics['Composite Quality'] = metrics['Composite Quality Score']
    
    # Calculate point sizes based on Ligand Clashes
    if metrics['Ligand Clashes'].notna().any():
        clash_filled = metrics['Ligand Clashes'].fillna(0)
        clash_min = clash_filled.min()
        clash_max = clash_filled.max()
        
        if clash_max > clash_min:
            # Map to size range 8-25 for plotly
            metrics['marker_size'] = 8 + ((clash_filled - clash_min) / (clash_max - clash_min) * 17)
        else:
            metrics['marker_size'] = 15
    else:
        metrics['marker_size'] = 15
    
    # Create hover text with detailed information
    hover_text = []
    for idx, row in metrics.iterrows():
        ligand_id = row['Ligand ID'] if pd.notna(row['Ligand ID']) else 'NA'
        text = f"<b style='font-size:14px'>{row['Crystal Name']}</b> ({ligand_id})<br>"
        text += "─" * 30 + "<br>"
        text += f"<b>Composite Quality:</b> {row['Composite Quality']:.1f}/100<br>"
        text += "<br>"
        text += "<b>Fit Metrics:</b><br>"
        text += f"  • Ligand CC: {row['Ligand CC']:.3f}<br>"
        text += f"  • Ligand RSR: {row['Ligand RSR']:.3f}<br>"
        text += f"  → <b>Fit Quality Score: {row['Fit Quality']:.1f}/100</b><br>"
        text += "<br>"
        text += "<b>Geometry Metrics:</b><br>"
        text += f"  • Mogul Z Bond: {row['Mogul Z Bond']:.2f}<br>"
        text += f"  • Mogul Z Angle: {row['Mogul Z Angle']:.2f}<br>"
        text += f"  → <b>Geometry Quality Score: {row['Geometry Quality']:.1f}/100</b><br>"
        if pd.notna(row['Ligand Clashes']):
            text += "<br>"
            text += f"<b>Ligand Clashes:</b> {int(row['Ligand Clashes'])}"
        if pd.notna(row['High Resolution (Å)']):
            text += "<br>"
            text += f"<b>High Resolution:</b> {row['High Resolution (Å)']:.2f} Å"
        if pd.notna(row['R']) and pd.notna(row['Rfree']):
            text += "<br>"
            text += f"<b>R/Rfree:</b> {row['R']:.3f} / {row['Rfree']:.3f}"
        elif pd.notna(row['R']):
            text += "<br>"
            text += f"<b>R:</b> {row['R']:.3f}"
        elif pd.notna(row['Rfree']):
            text += "<br>"
            text += f"<b>Rfree:</b> {row['Rfree']:.3f}"
        hover_text.append(text)
    
    # Create the interactive plot
    fig = go.Figure()
    
    # Add quality zone backgrounds with harmonious pastel colors
    # Poor zone (red)
    fig.add_shape(type="rect", x0=0, y0=0, x1=50, y1=50, 
                  fillcolor="#FFB3B3", opacity=0.15, layer="below", line_width=0)
    
    # OK zones (yellow)
    fig.add_shape(type="rect", x0=50, y0=0, x1=75, y1=50, 
                  fillcolor="#FFF4C4", opacity=0.2, layer="below", line_width=0)
    fig.add_shape(type="rect", x0=0, y0=50, x1=50, y1=75, 
                  fillcolor="#FFF4C4", opacity=0.2, layer="below", line_width=0)
    fig.add_shape(type="rect", x0=50, y0=50, x1=75, y1=75, 
                  fillcolor="#FFF4C4", opacity=0.2, layer="below", line_width=0)
    
    # Good zones (light yellow/green transition)
    fig.add_shape(type="rect", x0=75, y0=0, x1=100, y1=50, 
                  fillcolor="#FFF4C4", opacity=0.15, layer="below", line_width=0)
    fig.add_shape(type="rect", x0=0, y0=75, x1=50, y1=100, 
                  fillcolor="#FFF4C4", opacity=0.15, layer="below", line_width=0)
    fig.add_shape(type="rect", x0=50, y0=75, x1=75, y1=100, 
                  fillcolor="#D4F1D4", opacity=0.2, layer="below", line_width=0)
    fig.add_shape(type="rect", x0=75, y0=50, x1=100, y1=75, 
                  fillcolor="#D4F1D4", opacity=0.2, layer="below", line_width=0)
    
    # Excellent zone (green)
    fig.add_shape(type="rect", x0=75, y0=75, x1=100, y1=100, 
                  fillcolor="#B8E6B8", opacity=0.25, layer="below", line_width=0)
    
    # Add reference lines with better styling
    # Position annotations at edges to avoid data point overlap
    fig.add_hline(y=75, line_dash="dash", line_color="darkgreen", opacity=0.6, line_width=2.5)
    fig.add_vline(x=75, line_dash="dash", line_color="darkgreen", opacity=0.6, line_width=2.5)
    fig.add_hline(y=50, line_dash="dash", line_color="orange", opacity=0.6, line_width=2.5)
    fig.add_vline(x=50, line_dash="dash", line_color="orange", opacity=0.6, line_width=2.5)
    
    # Add line labels as separate annotations positioned at edges (moved further right to avoid clashing)
    fig.add_annotation(x=103, y=87.5, text="Excellent<br>geometry<br>(>75)", 
                      showarrow=False, xanchor='left', font=dict(size=12, color="darkgreen"))
    fig.add_annotation(x=87.5, y=103, text="Excellent fit (>75)", 
                      showarrow=False, yanchor='bottom', font=dict(size=12, color="darkgreen"))
    fig.add_annotation(x=103, y=62.5, text="OK<br>geometry<br>(50-75)", 
                      showarrow=False, xanchor='left', font=dict(size=12, color="orange"))
    fig.add_annotation(x=62.5, y=103, text="OK fit (50-75)", 
                      showarrow=False, yanchor='bottom', font=dict(size=12, color="orange"))
    fig.add_annotation(x=103, y=25, text="Bad<br>geometry<br>(<50)", 
                      showarrow=False, xanchor='left', font=dict(size=12, color="red"))
    fig.add_annotation(x=25, y=103, text="Bad fit (<50)", 
                      showarrow=False, yanchor='bottom', font=dict(size=12, color="red"))

    # Main scatter plot
    fig.add_trace(go.Scatter(
        x=metrics['Fit Quality'],
        y=metrics['Geometry Quality'],
        mode='markers',
        marker=dict(
            size=metrics['marker_size'],
            color=metrics['Composite Quality'],
            colorscale='RdYlGn',
            cmin=0,
            cmax=100,
            colorbar=dict(
                title="Overall<br>Quality<br>Score",
                thickness=20,
                len=0.7,
                x=1.02,
                xanchor='left'
            ),
            line=dict(color='black', width=1),
            opacity=0.8
        ),
        text=hover_text,
        hovertemplate='%{text}<extra></extra>',
        name='Datasets'
    ))
    
    # Calculate statistics
    total = len(metrics)
    excellent = len(metrics[(metrics['Fit Quality'] >= 75) & (metrics['Geometry Quality'] >= 75)])
    good_fit = len(metrics[metrics['Fit Quality'] >= 75])
    good_geom = len(metrics[metrics['Geometry Quality'] >= 75])
    poor = len(metrics[(metrics['Fit Quality'] < 50) | (metrics['Geometry Quality'] < 50)])
    
    stats_text = f'<b style="font-size:12px">Dataset Statistics (n={total})</b><br>'
    stats_text += '─' * 25 + '<br>'
    stats_text += f'<span style="color:darkgreen">●</span> Excellent (both ≥75): <b>{excellent}</b> ({excellent/total*100:.1f}%)<br>'
    stats_text += f'<span style="color:green">●</span> Good fit (≥75): <b>{good_fit}</b> ({good_fit/total*100:.1f}%)<br>'
    stats_text += f'<span style="color:green">●</span> Good geometry (≥75): <b>{good_geom}</b> ({good_geom/total*100:.1f}%)<br>'
    stats_text += f'<span style="color:red">●</span> Poor (<50 either): <b>{poor}</b> ({poor/total*100:.1f}%)'
    
    # Add annotations - Marker Size Key at top left, Dataset Statistics at bottom left
    # Marker Size Key (top left)
    clash_legend_text = '<b style="font-size:10px">Marker Size Key</b><br>'
    clash_legend_text += '─' * 18 + '<br>'
    clash_legend_text += 'Size = <b>Ligand Clashes</b><br>'
    clash_legend_text += '<span style="font-size:9px">●</span> Low  '
    clash_legend_text += '<span style="font-size:12px">●</span> Med  '
    clash_legend_text += '<span style="font-size:15px">●</span> High'
    
    fig.add_annotation(x=2, y=98, xref="x", yref="y",
                      text=clash_legend_text, showarrow=False,
                      xanchor='left', yanchor='top',
                      bgcolor='#E6F3FF', opacity=0.90,
                      bordercolor='#4A90E2', borderwidth=2, borderpad=6,
                      font=dict(size=9))
    
    # Dataset Statistics (bottom left, very low)
    fig.add_annotation(x=2, y=2, xref="x", yref="y",
                      text=stats_text, showarrow=False,
                      xanchor='left', yanchor='bottom',
                      bgcolor='#FFF9E6', opacity=0.90,
                      bordercolor='#DAA520', borderwidth=2, borderpad=6,
                      font=dict(size=10))
    
    # Update layout
    fig.update_layout(
        title=dict(
            text='<b>2D Ligand Quality Assessment</b><br><sub>Hover over points for detailed metrics | Size indicates Ligand Clashes</sub>',
            x=0.5,
            xanchor='center',
            font=dict(size=19)
        ),
        xaxis=dict(
            title=dict(
                text='<b>Ligand Fit Quality Score</b><br><sub>(Combined Ligand CC + RSR)</sub>',
                font=dict(size=16)
            ),
            range=[-3, 110],
            gridcolor='lightgray',
            gridwidth=0.65,
            showgrid=True,
            tickmode='array',
            tickvals=[0, 25, 50, 75, 100],
            ticktext=['0', '25', '<b>50</b>', '<b>75</b>', '100'],
            tickfont=dict(size=13),
            zeroline=False
        ),
        yaxis=dict(
            title=dict(
                text='<b>Geometry Quality Score</b><br><sub>(Combined Mogul Z bond + angle)</sub>',
                font=dict(size=16)
            ),
            range=[-3, 110],
            gridcolor='lightgray',
            gridwidth=0.65,
            showgrid=True,
            tickmode='array',
            tickvals=[0, 25, 50, 75, 100],
            ticktext=['0', '25', '<b>50</b>', '<b>75</b>', '100'],
            tickfont=dict(size=13),
            zeroline=False
        ),
        width=900,
        height=600,
        hovermode='closest',
        plot_bgcolor='white',
        showlegend=False,
        margin=dict(l=80, r=80, t=120, b=80),
        hoverlabel=dict(
            bgcolor="white",
            font_size=14,
            font_family="Arial",
            bordercolor="black"
        )
    )
    
    # Add modebar buttons for better interaction
    config = {
        'displayModeBar': True,
        'displaylogo': False,
        'modeBarButtonsToAdd': ['drawline', 'drawopenpath', 'eraseshape'],
        'modeBarButtonsToRemove': ['select2d', 'lasso2d'],
        'toImageButtonOptions': {
            'format': 'png',
            'filename': 'ligand_quality_plot',
            'height': 950,
            'width': 1100,
            'scale': 2
        }
    }
    
    # Save plot as HTML with config and custom title
    output_path = os.path.join(output_dir, f'{base_name}_ligand_quality_2d.html')
    
    import plotly.io as pio
    html_string = pio.to_html(fig, config=config,
                              include_plotlyjs='cdn',
                              full_html=True)
    # Replace the default title
    html_string = html_string.replace('<head><meta charset="utf-8" /></head>',
                                     f'<head><meta charset="utf-8" /><title>2D Ligand Quality Assessment</title></head>')
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_string)
    
    logging.info(f"Generated interactive 2D ligand quality plot: {output_path}")
    print(f"Generated interactive 2D quality plot: {os.path.basename(output_path)}")
    
    return output_path


def generate_pdb_percentile_quality_plot(results: List[Dict[str, Any]], output_dir: str, base_name: str) -> str:
    """
    Generate an interactive 2D PCA Percentile Plot showing all datasets' positions relative to the PDB archive.
    
    Creates an interactive scatter plot with:
    - X-axis: Ligand Fit PCA Percentile (PDB)
    - Y-axis: Ligand Geometry PCA Percentile (PDB)
    - Hover: Shows Crystal Name and all quality metrics
    
    Args:
        results: List of result dictionaries
        output_dir: Directory to save plot
        base_name: Base name for output file
        
    Returns:
        Path to generated HTML plot or 'NA' if insufficient data
    """
    if not HAS_PLOTLY:
        logging.warning("plotly not available, skipping interactive PDB percentile plot")
        return 'NA'
    
    df = pd.DataFrame(results)
    
    # Extract relevant metrics
    metrics = df[['Crystal Name', 'Ligand ID', 'Ligand Fit PCA Percentile (PDB)', 'Ligand Geometry PCA Percentile (PDB)',
                  'Clash Score', 'Ligand CC', 'Ligand RSR', 'Mogul Z Bond', 'Mogul Z Angle', 
                  'R', 'Rfree', 'High Resolution (Å)']].copy()
    
    # Convert 'NA' strings to numeric
    for col in ['Ligand Fit PCA Percentile (PDB)', 'Ligand Geometry PCA Percentile (PDB)', 
                'Clash Score', 'Ligand CC', 'Ligand RSR', 'Mogul Z Bond', 'Mogul Z Angle', 
                'R', 'Rfree', 'High Resolution (Å)']:
        metrics[col] = pd.to_numeric(metrics[col], errors='coerce')
    
    # Remove rows with missing critical data
    metrics = metrics.dropna(subset=['Ligand Fit PCA Percentile (PDB)', 'Ligand Geometry PCA Percentile (PDB)'])
    
    if len(metrics) < 2:
        logging.warning("Insufficient data for PDB percentile plot (need at least 2 datasets with PDB percentiles)")
        return 'NA'
    
    # Create hover text with detailed information
    hover_text = []
    for idx, row in metrics.iterrows():
        ligand_id = row['Ligand ID'] if pd.notna(row['Ligand ID']) else 'NA'
        text = f"<b style='font-size:14px'>{row['Crystal Name']}</b> ({ligand_id})<br>"
        text += "─" * 30 + "<br>"
        text += "<b>PDB Percentiles:</b><br>"
        text += f"  • Fit Percentile: {row['Ligand Fit PCA Percentile (PDB)']:.1f}%<br>"
        text += f"  • Geometry Percentile: {row['Ligand Geometry PCA Percentile (PDB)']:.1f}%<br>"
        text += "<br>"
        text += "<b>Raw Metrics:</b><br>"
        if pd.notna(row['Clash Score']):
            text += f"  • Clash Score: {row['Clash Score']:.2f}<br>"
        if pd.notna(row['Ligand CC']):
            text += f"  • Ligand CC: {row['Ligand CC']:.3f}<br>"
        if pd.notna(row['Ligand RSR']):
            text += f"  • Ligand RSR: {row['Ligand RSR']:.3f}<br>"
        if pd.notna(row['Mogul Z Bond']):
            text += f"  • Mogul Z Bond: {row['Mogul Z Bond']:.2f}<br>"
        if pd.notna(row['Mogul Z Angle']):
            text += f"  • Mogul Z Angle: {row['Mogul Z Angle']:.2f}<br>"
        if pd.notna(row['High Resolution (Å)']):
            text += f"<b>High Resolution:</b> {row['High Resolution (Å)']:.2f} Å<br>"
        if pd.notna(row['R']) and pd.notna(row['Rfree']):
            text += f"<b>R/Rfree:</b> {row['R']:.3f} / {row['Rfree']:.3f}"
        elif pd.notna(row['R']):
            text += f"<b>R:</b> {row['R']:.3f}"
        elif pd.notna(row['Rfree']):
            text += f"<b>Rfree:</b> {row['Rfree']:.3f}"
        hover_text.append(text)
    
    # Create the interactive plot
    fig = go.Figure()
    
    # Add quality zone backgrounds with harmonious pastel colors
    # Poor zone (bottom-left: poor fit & poor geometry) - red
    fig.add_shape(type="rect", x0=0, y0=0, x1=50, y1=50, 
                  fillcolor="#FFB3B3", opacity=0.15, layer="below", line_width=0)
    
    # OK zones (yellow)
    fig.add_shape(type="rect", x0=50, y0=0, x1=75, y1=50, 
                  fillcolor="#FFF4C4", opacity=0.2, layer="below", line_width=0)
    fig.add_shape(type="rect", x0=0, y0=50, x1=50, y1=75, 
                  fillcolor="#FFF4C4", opacity=0.2, layer="below", line_width=0)
    fig.add_shape(type="rect", x0=50, y0=50, x1=75, y1=75, 
                  fillcolor="#FFF4C4", opacity=0.2, layer="below", line_width=0)
    
    # Good zones (light yellow/green transition)
    fig.add_shape(type="rect", x0=75, y0=0, x1=100, y1=50, 
                  fillcolor="#FFF4C4", opacity=0.15, layer="below", line_width=0)
    fig.add_shape(type="rect", x0=0, y0=75, x1=100, y1=100, 
                  fillcolor="#FFF4C4", opacity=0.15, layer="below", line_width=0)
    fig.add_shape(type="rect", x0=50, y0=75, x1=75, y1=100, 
                  fillcolor="#D4F1D4", opacity=0.2, layer="below", line_width=0)
    fig.add_shape(type="rect", x0=75, y0=50, x1=100, y1=75, 
                  fillcolor="#D4F1D4", opacity=0.2, layer="below", line_width=0)
    
    # Excellent zone (top-right: excellent fit & geometry) - green
    fig.add_shape(type="rect", x0=75, y0=75, x1=100, y1=100, 
                  fillcolor="#B8E6B8", opacity=0.25, layer="below", line_width=0)
    
    # Add reference lines
    fig.add_hline(y=75, line_dash="dash", line_color="darkgreen", opacity=0.6, line_width=2.5)
    fig.add_vline(x=75, line_dash="dash", line_color="darkgreen", opacity=0.6, line_width=2.5)
    fig.add_hline(y=50, line_dash="dash", line_color="orange", opacity=0.6, line_width=2.5)
    fig.add_vline(x=50, line_dash="dash", line_color="orange", opacity=0.6, line_width=2.5)
    
    # Add line labels as annotations positioned at edges
    fig.add_annotation(x=103, y=87.5, text="Top 25%<br>geometry", 
                      showarrow=False, xanchor='left', font=dict(size=12, color="darkgreen"))
    fig.add_annotation(x=87.5, y=103, text="Top 25% fit", 
                      showarrow=False, yanchor='bottom', font=dict(size=12, color="darkgreen"))
    fig.add_annotation(x=103, y=62.5, text="Median-<br>Top 25%<br>geometry", 
                      showarrow=False, xanchor='left', font=dict(size=12, color="orange"))
    fig.add_annotation(x=62.5, y=103, text="Median - Top 25% fit", 
                      showarrow=False, yanchor='bottom', font=dict(size=12, color="orange"))
    fig.add_annotation(x=103, y=25, text="Bottom<br>50%<br>geometry", 
                      showarrow=False, xanchor='left', font=dict(size=12, color="red"))
    fig.add_annotation(x=25, y=103, text="Bottom 50% fit", 
                      showarrow=False, yanchor='bottom', font=dict(size=12, color="red"))

    # Main scatter plot with clash-based coloring
    # Fill NaN clash scores with 0 for coloring
    clash_values = metrics['Clash Score'].fillna(0)
    
    fig.add_trace(go.Scatter(
        x=metrics['Ligand Fit PCA Percentile (PDB)'],
        y=metrics['Ligand Geometry PCA Percentile (PDB)'],
        mode='markers',
        marker=dict(
            size=12,
            color=clash_values,
            colorscale='RdYlGn_r',  # Reversed: Green (low clash score) to Red (high clash score)
            cmin=0,
            cmax=max(clash_values.max(), 5),  # At least 5 for scale
            colorbar=dict(
                title="Clash<br>Score",
                thickness=15,
                len=0.6,
                x=1.02,
                xanchor='left'
            ),
            line=dict(color='black', width=1),
            opacity=0.8
        ),
        text=hover_text,
        hovertemplate='%{text}<extra></extra>',
        name='Datasets'
    ))
    
    # Calculate statistics
    total = len(metrics)
    excellent = len(metrics[(metrics['Ligand Fit PCA Percentile (PDB)'] >= 75) & 
                           (metrics['Ligand Geometry PCA Percentile (PDB)'] >= 75)])
    good_fit = len(metrics[metrics['Ligand Fit PCA Percentile (PDB)'] >= 75])
    good_geom = len(metrics[metrics['Ligand Geometry PCA Percentile (PDB)'] >= 75])
    poor = len(metrics[(metrics['Ligand Fit PCA Percentile (PDB)'] < 50) | 
                      (metrics['Ligand Geometry PCA Percentile (PDB)'] < 50)])
    
    stats_text = f'<b style="font-size:12px">Dataset Statistics (n={total})</b><br>'
    stats_text += '─' * 25 + '<br>'
    stats_text += f'<span style="color:darkgreen">●</span> Top 25% both: <b>{excellent}</b> ({excellent/total*100:.1f}%)<br>'
    stats_text += f'<span style="color:green">●</span> Top 25% fit: <b>{good_fit}</b> ({good_fit/total*100:.1f}%)<br>'
    stats_text += f'<span style="color:green">●</span> Top 25% geometry: <b>{good_geom}</b> ({good_geom/total*100:.1f}%)<br>'
    stats_text += f'<span style="color:red">●</span> Bottom 50% either: <b>{poor}</b> ({poor/total*100:.1f}%)'
    
    # Add annotations
    # Dataset Statistics (bottom left)
    fig.add_annotation(x=2, y=2, xref="x", yref="y",
                      text=stats_text, showarrow=False,
                      xanchor='left', yanchor='bottom',
                      bgcolor='#FFF9E6', opacity=0.90,
                      bordercolor='#DAA520', borderwidth=2, borderpad=6,
                      font=dict(size=10))
    
    # Update layout
    fig.update_layout(
        title=dict(
            text='<b>PCA Percentile Plot (PDB)</b><br><sub>Position relative to PDB | Hover for details | Color indicates clashes</sub>',
            x=0.5,
            xanchor='center',
            font=dict(size=19)
        ),
        xaxis=dict(
            title=dict(
                text='<b>Ligand Fit PCA Percentile (PDB)</b><br><sub>Percentile rank vs PDB</sub>',
                font=dict(size=16)
            ),
            range=[-3, 110],
            gridcolor='lightgray',
            gridwidth=0.65,
            showgrid=True,
            tickmode='array',
            tickvals=[0, 25, 50, 75, 100],
            ticktext=['0', '25', '<b>50</b>', '<b>75</b>', '100'],
            tickfont=dict(size=13),
            zeroline=False
        ),
        yaxis=dict(
            title=dict(
                text='<b>Ligand Geometry PCA Percentile (PDB)</b><br><sub>Percentile rank vs PDB</sub>',
                font=dict(size=16)
            ),
            range=[-3, 110],
            gridcolor='lightgray',
            gridwidth=0.65,
            showgrid=True,
            tickmode='array',
            tickvals=[0, 25, 50, 75, 100],
            ticktext=['0', '25', '<b>50</b>', '<b>75</b>', '100'],
            tickfont=dict(size=13),
            zeroline=False
        ),
        width=900,
        height=600,
        hovermode='closest',
        plot_bgcolor='white',
        showlegend=False,
        margin=dict(l=80, r=80, t=120, b=80),
        hoverlabel=dict(
            bgcolor="white",
            font_size=14,
            font_family="Arial",
            bordercolor="black"
        )
    )
    
    # Add modebar buttons for better interaction
    config = {
        'displayModeBar': True,
        'displaylogo': False,
        'modeBarButtonsToAdd': ['drawline', 'drawopenpath', 'eraseshape'],
        'modeBarButtonsToRemove': ['select2d', 'lasso2d'],
        'toImageButtonOptions': {
            'format': 'png',
            'filename': 'pdb_percentile_quality_plot',
            'height': 950,
            'width': 1100,
            'scale': 2
        }
    }
    
    # Save plot as HTML with config and custom title
    output_path = os.path.join(output_dir, f'{base_name}_pdb_percentile_quality_2d.html')
    
    import plotly.io as pio
    html_string = pio.to_html(fig, config=config,
                              include_plotlyjs='cdn',
                              full_html=True)
    # Replace the default title
    html_string = html_string.replace('<head><meta charset="utf-8" /></head>',
                                     f'<head><meta charset="utf-8" /><title>PCA Percentile Plot (PDB)</title></head>')
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_string)
    
    logging.info(f"Generated interactive PCA Percentile Plot (PDB): {output_path}")
    print(f"Generated interactive PCA Percentile Plot (PDB): {os.path.basename(output_path)}")
    
    return output_path


def generate_batch_percentile_plot(results: List[Dict[str, Any]], output_dir: str, base_name: str) -> str:
    """
    Generate an interactive 2D PCA Percentile Plot showing all datasets' positions within the batch.
    
    Creates an interactive scatter plot with:
    - X-axis: Ligand Fit PCA Percentile (Batch)
    - Y-axis: Ligand Geometry PCA Percentile (Batch)
    - Hover: Shows Crystal Name and all quality metrics
    
    Args:
        results: List of result dictionaries
        output_dir: Directory to save plot
        base_name: Base name for output file
        
    Returns:
        Path to generated HTML plot or 'NA' if insufficient data
    """
    if not HAS_PLOTLY:
        logging.warning("plotly not available, skipping interactive batch percentile plot")
        return 'NA'
    
    df = pd.DataFrame(results)
    
    # Extract relevant metrics
    metrics = df[['Crystal Name', 'Ligand Fit PCA Percentile (Batch)', 'Ligand Geometry PCA Percentile (Batch)',
                  'Clash Score', 'Ligand CC', 'Ligand RSR', 'Mogul Z Bond', 'Mogul Z Angle', 
                  'R', 'Rfree', 'High Resolution (Å)']].copy()
    
    # Convert 'NA' strings to numeric
    for col in ['Ligand Fit PCA Percentile (Batch)', 'Ligand Geometry PCA Percentile (Batch)', 
                'Clash Score', 'Ligand CC', 'Ligand RSR', 'Mogul Z Bond', 'Mogul Z Angle', 
                'R', 'Rfree', 'High Resolution (Å)']:
        metrics[col] = pd.to_numeric(metrics[col], errors='coerce')
    
    # Remove rows with missing critical data
    metrics = metrics.dropna(subset=['Ligand Fit PCA Percentile (Batch)', 'Ligand Geometry PCA Percentile (Batch)'])
    
    if len(metrics) < 2:
        logging.warning("Insufficient data for batch percentile plot (need at least 2 datasets with batch percentiles)")
        return 'NA'
    
    # Create hover text with detailed information
    hover_text = []
    for idx, row in metrics.iterrows():
        ligand_id = row['Ligand ID'] if pd.notna(row['Ligand ID']) else 'NA'
        text = f"<b style='font-size:14px'>{row['Crystal Name']}</b> ({ligand_id})<br>"
        text += "─" * 30 + "<br>"
        text += "<b>Batch Percentiles:</b><br>"
        text += f"  • Fit Percentile: {row['Ligand Fit PCA Percentile (Batch)']:.1f}%<br>"
        text += f"  • Geometry Percentile: {row['Ligand Geometry PCA Percentile (Batch)']:.1f}%<br>"
        text += "<br>"
        text += "<b>Raw Metrics:</b><br>"
        if pd.notna(row['Clash Score']):
            text += f"  • Clash Score: {row['Clash Score']:.2f}<br>"
        if pd.notna(row['Ligand CC']):
            text += f"  • Ligand CC: {row['Ligand CC']:.3f}<br>"
        if pd.notna(row['Ligand RSR']):
            text += f"  • Ligand RSR: {row['Ligand RSR']:.3f}<br>"
        if pd.notna(row['Mogul Z Bond']):
            text += f"  • Mogul Z Bond: {row['Mogul Z Bond']:.2f}<br>"
        if pd.notna(row['Mogul Z Angle']):
            text += f"  • Mogul Z Angle: {row['Mogul Z Angle']:.2f}<br>"
        if pd.notna(row['High Resolution (Å)']):
            text += f"<b>High Resolution:</b> {row['High Resolution (Å)']:.2f} Å<br>"
        if pd.notna(row['R']) and pd.notna(row['Rfree']):
            text += f"<b>R/Rfree:</b> {row['R']:.3f} / {row['Rfree']:.3f}"
        elif pd.notna(row['R']):
            text += f"<b>R:</b> {row['R']:.3f}"
        elif pd.notna(row['Rfree']):
            text += f"<b>Rfree:</b> {row['Rfree']:.3f}"
        hover_text.append(text)
    
    # Create the interactive plot
    fig = go.Figure()
    
    # Add quality zone backgrounds with harmonious pastel colors
    # Poor zone (bottom-left: poor fit & poor geometry) - red
    fig.add_shape(type="rect", x0=0, y0=0, x1=50, y1=50, 
                  fillcolor="#FFB3B3", opacity=0.15, layer="below", line_width=0)
    
    # OK zones (yellow)
    fig.add_shape(type="rect", x0=50, y0=0, x1=75, y1=50, 
                  fillcolor="#FFF4C4", opacity=0.2, layer="below", line_width=0)
    fig.add_shape(type="rect", x0=0, y0=50, x1=50, y1=75, 
                  fillcolor="#FFF4C4", opacity=0.2, layer="below", line_width=0)
    fig.add_shape(type="rect", x0=50, y0=50, x1=75, y1=75, 
                  fillcolor="#FFF4C4", opacity=0.2, layer="below", line_width=0)
    
    # Good zones (light yellow/green transition)
    fig.add_shape(type="rect", x0=75, y0=0, x1=100, y1=50, 
                  fillcolor="#FFF4C4", opacity=0.15, layer="below", line_width=0)
    fig.add_shape(type="rect", x0=0, y0=75, x1=100, y1=100, 
                  fillcolor="#FFF4C4", opacity=0.15, layer="below", line_width=0)
    fig.add_shape(type="rect", x0=50, y0=75, x1=75, y1=100, 
                  fillcolor="#D4F1D4", opacity=0.2, layer="below", line_width=0)
    fig.add_shape(type="rect", x0=75, y0=50, x1=100, y1=75, 
                  fillcolor="#D4F1D4", opacity=0.2, layer="below", line_width=0)
    
    # Excellent zone (top-right: excellent fit & geometry) - green
    fig.add_shape(type="rect", x0=75, y0=75, x1=100, y1=100, 
                  fillcolor="#B8E6B8", opacity=0.25, layer="below", line_width=0)
    
    # Add reference lines
    fig.add_hline(y=75, line_dash="dash", line_color="darkgreen", opacity=0.6, line_width=2.5)
    fig.add_vline(x=75, line_dash="dash", line_color="darkgreen", opacity=0.6, line_width=2.5)
    fig.add_hline(y=50, line_dash="dash", line_color="orange", opacity=0.6, line_width=2.5)
    fig.add_vline(x=50, line_dash="dash", line_color="orange", opacity=0.6, line_width=2.5)
    
    # Add line labels as annotations positioned at edges
    fig.add_annotation(x=103, y=87.5, text="Top 25%<br>geometry", 
                      showarrow=False, xanchor='left', font=dict(size=12, color="darkgreen"))
    fig.add_annotation(x=87.5, y=103, text="Top 25% fit", 
                      showarrow=False, yanchor='bottom', font=dict(size=12, color="darkgreen"))
    fig.add_annotation(x=103, y=62.5, text="Median-<br>Top 25%<br>geometry", 
                      showarrow=False, xanchor='left', font=dict(size=12, color="orange"))
    fig.add_annotation(x=62.5, y=103, text="Median - Top 25% fit", 
                      showarrow=False, yanchor='bottom', font=dict(size=12, color="orange"))
    fig.add_annotation(x=103, y=25, text="Bottom<br>50%<br>geometry", 
                      showarrow=False, xanchor='left', font=dict(size=12, color="red"))
    fig.add_annotation(x=25, y=103, text="Bottom 50% fit", 
                      showarrow=False, yanchor='bottom', font=dict(size=12, color="red"))

    # Main scatter plot with clash-based coloring
    # Fill NaN clash scores with 0 for coloring
    clash_values = metrics['Clash Score'].fillna(0)
    
    fig.add_trace(go.Scatter(
        x=metrics['Ligand Fit PCA Percentile (Batch)'],
        y=metrics['Ligand Geometry PCA Percentile (Batch)'],
        mode='markers',
        marker=dict(
            size=12,
            color=clash_values,
            colorscale='RdYlGn_r',  # Reversed: Green (low clash score) to Red (high clash score)
            cmin=0,
            cmax=max(clash_values.max(), 5),  # At least 5 for scale
            colorbar=dict(
                title="Clash<br>Score",
                thickness=15,
                len=0.6,
                x=1.02,
                xanchor='left'
            ),
            line=dict(color='black', width=1),
            opacity=0.8
        ),
        text=hover_text,
        hovertemplate='%{text}<extra></extra>',
        name='Datasets'
    ))
    
    # Calculate statistics
    total = len(metrics)
    excellent = len(metrics[(metrics['Ligand Fit PCA Percentile (Batch)'] >= 75) & 
                           (metrics['Ligand Geometry PCA Percentile (Batch)'] >= 75)])
    good_fit = len(metrics[metrics['Ligand Fit PCA Percentile (Batch)'] >= 75])
    good_geom = len(metrics[metrics['Ligand Geometry PCA Percentile (Batch)'] >= 75])
    poor = len(metrics[(metrics['Ligand Fit PCA Percentile (Batch)'] < 50) | 
                      (metrics['Ligand Geometry PCA Percentile (Batch)'] < 50)])
    
    stats_text = f'<b style="font-size:12px">Dataset Statistics (n={total})</b><br>'
    stats_text += '─' * 25 + '<br>'
    stats_text += f'<span style="color:darkgreen">●</span> Top 25% both: <b>{excellent}</b> ({excellent/total*100:.1f}%)<br>'
    stats_text += f'<span style="color:green">●</span> Top 25% fit: <b>{good_fit}</b> ({good_fit/total*100:.1f}%)<br>'
    stats_text += f'<span style="color:green">●</span> Top 25% geometry: <b>{good_geom}</b> ({good_geom/total*100:.1f}%)<br>'
    stats_text += f'<span style="color:red">●</span> Bottom 50% either: <b>{poor}</b> ({poor/total*100:.1f}%)'
    
    # Add annotations
    # Dataset Statistics (bottom left)
    fig.add_annotation(x=2, y=2, xref="x", yref="y",
                      text=stats_text, showarrow=False,
                      xanchor='left', yanchor='bottom',
                      bgcolor='#FFF9E6', opacity=0.90,
                      bordercolor='#DAA520', borderwidth=2, borderpad=6,
                      font=dict(size=10))
    
    # Update layout
    fig.update_layout(
        title=dict(
            text='<b>PCA Percentile Plot (Batch)</b><br><sub>Position within this batch | Hover for details | Color indicates clashes</sub>',
            x=0.5,
            xanchor='center',
            font=dict(size=19)
        ),
        xaxis=dict(
            title=dict(
                text='<b>Ligand Fit PCA Percentile (Batch)</b><br><sub>Percentile rank within this batch</sub>',
                font=dict(size=16)
            ),
            range=[-3, 110],
            gridcolor='lightgray',
            gridwidth=0.65,
            showgrid=True,
            tickmode='array',
            tickvals=[0, 25, 50, 75, 100],
            ticktext=['0', '25', '<b>50</b>', '<b>75</b>', '100'],
            tickfont=dict(size=13),
            zeroline=False
        ),
        yaxis=dict(
            title=dict(
                text='<b>Ligand Geometry PCA Percentile (Batch)</b><br><sub>Percentile rank within this batch</sub>',
                font=dict(size=16)
            ),
            range=[-3, 110],
            gridcolor='lightgray',
            gridwidth=0.65,
            showgrid=True,
            tickmode='array',
            tickvals=[0, 25, 50, 75, 100],
            ticktext=['0', '25', '<b>50</b>', '<b>75</b>', '100'],
            tickfont=dict(size=13),
            zeroline=False
        ),
        width=900,
        height=600,
        hovermode='closest',
        plot_bgcolor='white',
        showlegend=False,
        margin=dict(l=80, r=80, t=120, b=80),
        hoverlabel=dict(
            bgcolor="white",
            font_size=14,
            font_family="Arial",
            bordercolor="black"
        )
    )
    
    # Add modebar buttons for better interaction
    config = {
        'displayModeBar': True,
        'displaylogo': False,
        'modeBarButtonsToAdd': ['drawline', 'drawopenpath', 'eraseshape'],
        'modeBarButtonsToRemove': ['select2d', 'lasso2d'],
        'toImageButtonOptions': {
            'format': 'png',
            'filename': 'batch_percentile_plot',
            'height': 950,
            'width': 1100,
            'scale': 2
        }
    }
    
    # Save plot as HTML with config and custom title
    output_path = os.path.join(output_dir, f'{base_name}_batch_percentile_2d.html')
    
    import plotly.io as pio
    html_string = pio.to_html(fig, config=config,
                              include_plotlyjs='cdn',
                              full_html=True)
    # Replace the default title
    html_string = html_string.replace('<head><meta charset="utf-8" /></head>',
                                     f'<head><meta charset="utf-8" /><title>PCA Percentile Plot (Batch)</title></head>')
    with open(output_path, 'w', encoding='utf-8') as f:
        f.write(html_string)
    
    logging.info(f"Generated interactive PCA Percentile Plot (Batch): {output_path}")
    print(f"Generated interactive PCA Percentile Plot (Batch): {os.path.basename(output_path)}")
    
    return output_path


def generate_individual_pdb_percentile_plot(result: Dict[str, Any], output_dir: str) -> str:
    """
    Generate an individual interactive 2D scatter plot showing PDB percentile position for a single dataset.
    
    Args:
        result: Dictionary containing dataset results
        output_dir: Directory to save the plot
        
    Returns:
        Path to the generated plot file or 'NA' if plot could not be generated
    """
    if not HAS_PLOTLY:
        logging.warning("plotly not available, skipping individual PDB percentile plot")
        return 'NA'
        
    try:
        # Extract data
        dataset_name = result.get('Crystal Name', 'Unknown')
        ligand_id = result.get('Ligand ID', 'NA')
        fit_pdb = result.get('Ligand Fit PCA Percentile (PDB)', 'NA')
        geom_pdb = result.get('Ligand Geometry PCA Percentile (PDB)', 'NA')
        
        # Check if we have valid data
        if fit_pdb == 'NA' or geom_pdb == 'NA':
            logging.debug(f"Skipping PDB percentile plot for {dataset_name}: missing PDB percentile data")
            return 'NA'
        
        try:
            fit_val = float(fit_pdb)
            geom_val = float(geom_pdb)
        except (ValueError, TypeError):
            logging.debug(f"Skipping PDB percentile plot for {dataset_name}: invalid percentile values")
            return 'NA'
        
        # Create output directory for individual plots if it doesn't exist
        plots_dir = os.path.join(output_dir, 'pdb_percentile_plots')
        os.makedirs(plots_dir, exist_ok=True)
        
        # Create safe filename from dataset name
        safe_name = "".join(c if c.isalnum() or c in ('_', '-') else '_' for c in dataset_name)
        output_path = os.path.join(plots_dir, f'{safe_name}_pdb_percentile.html')
        
        # Create the interactive plot
        fig = go.Figure()
        
        # Add quality zone backgrounds with harmonious pastel colors
        # Poor zone (bottom-left: poor fit & poor geometry) - red
        fig.add_shape(type="rect", x0=0, y0=0, x1=50, y1=50, 
                      fillcolor="#FFB3B3", opacity=0.15, layer="below", line_width=0)
        
        # OK zones (yellow)
        fig.add_shape(type="rect", x0=50, y0=0, x1=75, y1=50, 
                      fillcolor="#FFF4C4", opacity=0.2, layer="below", line_width=0)
        fig.add_shape(type="rect", x0=0, y0=50, x1=50, y1=75, 
                      fillcolor="#FFF4C4", opacity=0.2, layer="below", line_width=0)
        fig.add_shape(type="rect", x0=50, y0=50, x1=75, y1=75, 
                      fillcolor="#FFF4C4", opacity=0.2, layer="below", line_width=0)
        
        # Good zones (light yellow/green transition)
        fig.add_shape(type="rect", x0=75, y0=0, x1=100, y1=50, 
                      fillcolor="#FFF4C4", opacity=0.15, layer="below", line_width=0)
        fig.add_shape(type="rect", x0=0, y0=75, x1=100, y1=100, 
                      fillcolor="#FFF4C4", opacity=0.15, layer="below", line_width=0)
        fig.add_shape(type="rect", x0=50, y0=75, x1=75, y1=100, 
                      fillcolor="#D4F1D4", opacity=0.2, layer="below", line_width=0)
        fig.add_shape(type="rect", x0=75, y0=50, x1=100, y1=75, 
                      fillcolor="#D4F1D4", opacity=0.2, layer="below", line_width=0)
        
        # Excellent zone (top-right: excellent fit & geometry) - green
        fig.add_shape(type="rect", x0=75, y0=75, x1=100, y1=100, 
                      fillcolor="#B8E6B8", opacity=0.25, layer="below", line_width=0)
        
        # Add reference lines
        fig.add_hline(y=75, line_dash="dash", line_color="darkgreen", opacity=0.6, line_width=2.5)
        fig.add_vline(x=75, line_dash="dash", line_color="darkgreen", opacity=0.6, line_width=2.5)
        fig.add_hline(y=50, line_dash="dash", line_color="orange", opacity=0.6, line_width=2.5)
        fig.add_vline(x=50, line_dash="dash", line_color="orange", opacity=0.6, line_width=2.5)
        
        # Add line labels as annotations positioned at edges
        fig.add_annotation(x=103, y=87.5, text="Top 25%<br>geometry", 
                          showarrow=False, xanchor='left', font=dict(size=12, color="darkgreen"))
        fig.add_annotation(x=87.5, y=103, text="Top 25% fit", 
                          showarrow=False, yanchor='bottom', font=dict(size=12, color="darkgreen"))
        fig.add_annotation(x=103, y=62.5, text="Median-<br>Top 25%<br>geometry", 
                          showarrow=False, xanchor='left', font=dict(size=12, color="orange"))
        fig.add_annotation(x=62.5, y=103, text="Median - Top 25% fit", 
                          showarrow=False, yanchor='bottom', font=dict(size=12, color="orange"))
        fig.add_annotation(x=103, y=25, text="Bottom<br>50%<br>geometry", 
                          showarrow=False, xanchor='left', font=dict(size=12, color="red"))
        fig.add_annotation(x=25, y=103, text="Bottom 50% fit", 
                          showarrow=False, yanchor='bottom', font=dict(size=12, color="red"))
        
        # Create hover text
        hover_text = f"<b style='font-size:14px'>{dataset_name}</b><br>"
        if ligand_id not in ['NA', None, '']:
            hover_text += f"<b style='font-size:12px'>Ligand: {ligand_id}</b><br>"
        hover_text += "─" * 30 + "<br>"
        hover_text += "<b>PDB Archive Percentiles:</b><br>"
        hover_text += f"  • Fit Percentile: {fit_val:.1f}%<br>"
        hover_text += f"  • Geometry Percentile: {geom_val:.1f}%"
        
        # Plot the dataset point
        fig.add_trace(go.Scatter(
            x=[fit_val],
            y=[geom_val],
            mode='markers',
            marker=dict(
                size=15,
                color='#0d6efd',
                line=dict(color='black', width=2),
                opacity=0.9
            ),
            text=[hover_text],
            hovertemplate='%{text}<extra></extra>',
            name=dataset_name,
            showlegend=False
        ))
        
        # Update layout
        fig.update_layout(
            title=dict(
                text=f'<b>PCA Percentile Plot (PDB)</b><br><sub>{dataset_name}</sub>',
                x=0.5,
                xanchor='center',
                font=dict(size=18)
            ),
            xaxis=dict(
                title=dict(
                    text='<b>Ligand Fit PCA Percentile (PDB)</b><br><sub>Percentile rank vs PDB</sub>',
                    font=dict(size=14)
                ),
                range=[-3, 110],
                gridcolor='lightgray',
                gridwidth=0.65,
                showgrid=True,
                tickmode='array',
                tickvals=[0, 25, 50, 75, 100],
                ticktext=['0', '25', '<b>50</b>', '<b>75</b>', '100'],
                tickfont=dict(size=12),
                zeroline=False
            ),
            yaxis=dict(
                title=dict(
                    text='<b>Ligand Geometry PCA Percentile (PDB)</b><br><sub>Percentile rank vs PDB</sub>',
                    font=dict(size=14)
                ),
                range=[-3, 110],
                gridcolor='lightgray',
                gridwidth=0.65,
                showgrid=True,
                tickmode='array',
                tickvals=[0, 25, 50, 75, 100],
                ticktext=['0', '25', '<b>50</b>', '<b>75</b>', '100'],
                tickfont=dict(size=12),
                zeroline=False
            ),
            width=700,
            height=700,
            hovermode='closest',
            plot_bgcolor='white',
            margin=dict(l=80, r=80, t=100, b=80),
            hoverlabel=dict(
                bgcolor="white",
                font_size=13,
                font_family="Arial",
                bordercolor="black"
            )
        )
        
        # Add modebar config for consistent export options
        config = {
            'displayModeBar': True,
            'displaylogo': False,
            'modeBarButtonsToRemove': ['select2d', 'lasso2d'],
            'toImageButtonOptions': {
                'format': 'png',
                'filename': f'pdb_percentile_{safe_name}',
                'height': 800,
                'width': 800,
                'scale': 2
            }
        }
        
        # Save as HTML
        import plotly.io as pio
        html_string = pio.to_html(fig, config=config, include_plotlyjs='cdn', full_html=True)
        html_string = html_string.replace('<head><meta charset="utf-8" /></head>',
                                         f'<head><meta charset="utf-8" /><title>PCA Percentile Plot (PDB) - {dataset_name}</title></head>')
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_string)
        
        logging.debug(f"Generated PDB percentile plot for {dataset_name}: {output_path}")
        return output_path
        
    except Exception as e:
        logging.warning(f"Error generating PDB percentile plot for {result.get('Crystal Name', 'Unknown')}: {e}")
        return 'NA'


def generate_spider_plot(result: Dict[str, Any], output_dir: str) -> str:
    """
    Generate a spider/radar plot showing 6 key quality metrics for a dataset/ligand.
    
    Metrics:
    - Ligand CC (higher is better, scaled 0-1 → 0-100)
    - Ligand RSR (lower is better, inverted: (1-RSR)*100, capped at 0-1)
    - Mogul Z Angle (lower is better, inverted: max(0, 100-Z*20))
    - Mogul Z Bond (lower is better, inverted: max(0, 100-Z*20))
    - Ligand Clashes (lower is better, inverted based on max in batch)
    - B-factor Ratio (Ligand/Overall, closer to 1.0 is better)
    
    Args:
        result: Dictionary containing dataset results
        output_dir: Directory to save the plot
        
    Returns:
        Path to the generated plot file or 'NA' if plot could not be generated
    """
    if not HAS_PLOTLY:
        logging.warning("plotly not available, skipping spider plot")
        return 'NA'
    
    try:
        # Extract data
        dataset_name = result.get('Crystal Name', 'Unknown')
        ligand_id = result.get('Ligand ID', 'NA')
        
        # Get raw metric values
        ligand_cc = result.get('Ligand CC', 'NA')
        ligand_rsr = result.get('Ligand RSR', 'NA')
        mogul_z_angle = result.get('Mogul Z Angle', 'NA')
        mogul_z_bond = result.get('Mogul Z Bond', 'NA')
        ligand_clashes = result.get('Ligand Clashes', 'NA')
        ligand_mean_b = result.get('Ligand Mean B Factor', 'NA')
        overall_mean_b = result.get('Mean B Factor', 'NA')
        
        # Check if we have valid data for most metrics
        valid_count = 0
        metrics = []
        values = []
        
        # Process Ligand CC (0-1 → 0-100)
        if ligand_cc not in ['NA', None]:
            try:
                cc_val = float(ligand_cc)
                metrics.append('Ligand CC')
                values.append(cc_val * 100)  # Scale to 0-100
                valid_count += 1
            except (ValueError, TypeError):
                pass
        
        # Process Ligand RSR (invert: lower is better, typical range 0-0.5)
        if ligand_rsr not in ['NA', None]:
            try:
                rsr_val = float(ligand_rsr)
                # Invert and scale: RSR of 0 = 100, RSR of 0.5 = 0
                rsr_score = max(0, min(100, (1 - rsr_val * 2) * 100))
                metrics.append('Ligand Fit (RSR)')
                values.append(rsr_score)
                valid_count += 1
            except (ValueError, TypeError):
                pass
        
        # Process Mogul Z Angle (lower is better, typical range 0-5)
        if mogul_z_angle not in ['NA', None]:
            try:
                z_angle = float(mogul_z_angle)
                # Invert: Z=0 is 100, Z=5 is 0
                z_angle_score = max(0, min(100, 100 - z_angle * 20))
                metrics.append('Mogul Z Angle')
                values.append(z_angle_score)
                valid_count += 1
            except (ValueError, TypeError):
                pass
        
        # Process Mogul Z Bond (lower is better, typical range 0-5)
        if mogul_z_bond not in ['NA', None]:
            try:
                z_bond = float(mogul_z_bond)
                # Invert: Z=0 is 100, Z=5 is 0
                z_bond_score = max(0, min(100, 100 - z_bond * 20))
                metrics.append('Mogul Z Bond')
                values.append(z_bond_score)
                valid_count += 1
            except (ValueError, TypeError):
                pass
        
        # Process Ligand Clashes (lower is better, typical range 0-20)
        if ligand_clashes not in ['NA', None]:
            try:
                clashes = float(ligand_clashes)
                # Invert: 0 clashes = 100, 10+ clashes = 0
                clash_score = max(0, min(100, 100 - clashes * 10))
                metrics.append('Ligand Clashes')
                values.append(clash_score)
                valid_count += 1
            except (ValueError, TypeError):
                pass
        
        # Process B-factor Ratio (closer to 1.0 is better)
        if ligand_mean_b not in ['NA', None] and overall_mean_b not in ['NA', None]:
            try:
                lig_b = float(ligand_mean_b)
                overall_b = float(overall_mean_b)
                if overall_b > 0:
                    b_ratio = lig_b / overall_b
                    # Score: ratio of 1.0 = 100, ratio of 2.0 or 0.5 = 50, ratio of 3.0 or 0.33 = 0
                    if b_ratio >= 1.0:
                        b_score = max(0, 100 - (b_ratio - 1.0) * 50)
                    else:
                        b_score = max(0, 100 - (1.0 - b_ratio) * 100)
                    metrics.append('B-factor Ratio')
                    values.append(b_score)
                    valid_count += 1
            except (ValueError, TypeError, ZeroDivisionError):
                pass
        
        # Need at least 3 valid metrics to make a meaningful spider plot
        if valid_count < 3:
            logging.debug(f"Skipping spider plot for {dataset_name}: insufficient valid metrics ({valid_count}/6)")
            return 'NA'
        
        # Create output directory for spider plots if it doesn't exist
        plots_dir = os.path.join(output_dir, 'spider_plots')
        os.makedirs(plots_dir, exist_ok=True)
        
        # Create safe filename
        safe_name = "".join(c if c.isalnum() or c in ('_', '-') else '_' for c in dataset_name)
        if ligand_id not in ['NA', None, '']:
            safe_ligand = "".join(c if c.isalnum() or c in ('_', '-') else '_' for c in str(ligand_id))
            safe_name = f"{safe_name}_{safe_ligand}"
        output_path = os.path.join(plots_dir, f'{safe_name}_spider.html')
        
        # Create the spider plot
        fig = go.Figure()
        
        # Add the data trace
        fig.add_trace(go.Scatterpolar(
            r=values,
            theta=metrics,
            fill='toself',
            fillcolor='rgba(13, 110, 253, 0.3)',
            line=dict(color='#0d6efd', width=2),
            marker=dict(size=8, color='#0d6efd'),
            name=dataset_name,
            hovertemplate='<b>%{theta}</b><br>Score: %{r:.1f}<extra></extra>'
        ))
        
        # Add reference circles at 25, 50, 75
        for ref_val in [25, 50, 75]:
            fig.add_trace(go.Scatterpolar(
                r=[ref_val] * len(metrics),
                theta=metrics,
                mode='lines',
                line=dict(color='lightgray', width=1, dash='dot'),
                showlegend=False,
                hoverinfo='skip'
            ))
        
        # Update layout
        title_text = f'<b>Quality Metrics Spider Plot</b><br><sub>{dataset_name}'
        if ligand_id not in ['NA', None, '']:
            title_text += f' - {ligand_id}'
        title_text += '</sub>'
        
        fig.update_layout(
            polar=dict(
                radialaxis=dict(
                    visible=True,
                    range=[0, 100],
                    tickmode='array',
                    tickvals=[0, 25, 50, 75, 100],
                    ticktext=['0', '25', '50', '75', '100'],
                    tickfont=dict(size=12),
                    gridcolor='lightgray',
                    gridwidth=1
                ),
                angularaxis=dict(
                    tickfont=dict(size=12)
                ),
                bgcolor='rgba(240, 240, 240, 0.3)'
            ),
            title=dict(
                text=title_text,
                x=0.5,
                xanchor='center',
                font=dict(size=16)
            ),
            showlegend=False,
            width=700,
            height=700,
            margin=dict(l=80, r=80, t=100, b=80),
            hoverlabel=dict(
                bgcolor="white",
                font_size=13,
                font_family="Arial",
                bordercolor="black"
            )
        )
        
        # Add modebar config for consistent export options
        config = {
            'displayModeBar': True,
            'displaylogo': False,
            'modeBarButtonsToRemove': ['select2d', 'lasso2d'],
            'toImageButtonOptions': {
                'format': 'png',
                'filename': f'spider_{safe_name}',
                'height': 800,
                'width': 800,
                'scale': 2
            }
        }
        
        # Save as HTML
        import plotly.io as pio
        html_string = pio.to_html(fig, config=config, include_plotlyjs='cdn', full_html=True)
        html_string = html_string.replace('<head><meta charset="utf-8" /></head>',
                                         f'<head><meta charset="utf-8" /><title>Spider Plot - {dataset_name}</title></head>')
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_string)
        
        logging.debug(f"Generated spider plot for {dataset_name}: {output_path}")
        return output_path
        
    except Exception as e:
        logging.warning(f"Error generating spider plot for {result.get('Crystal Name', 'Unknown')}: {e}")
        return 'NA'


def generate_individual_batch_percentile_plot(result: Dict[str, Any], output_dir: str) -> str:
    """
    Generate an individual interactive 2D scatter plot showing batch percentile position for a single dataset.
    
    Args:
        result: Dictionary containing dataset results
        output_dir: Directory to save the plot
        
    Returns:
        Path to the generated plot file or 'NA' if plot could not be generated
    """
    if not HAS_PLOTLY:
        logging.warning("plotly not available, skipping individual batch percentile plot")
        return 'NA'
        
    try:
        # Extract data
        dataset_name = result.get('Crystal Name', 'Unknown')
        ligand_id = result.get('Ligand ID', 'NA')
        fit_batch = result.get('Ligand Fit PCA Percentile (Batch)', 'NA')
        geom_batch = result.get('Ligand Geometry PCA Percentile (Batch)', 'NA')
        
        # Check if we have valid data
        if fit_batch == 'NA' or geom_batch == 'NA':
            logging.debug(f"Skipping batch percentile plot for {dataset_name}: missing batch percentile data")
            return 'NA'
        
        try:
            fit_val = float(fit_batch)
            geom_val = float(geom_batch)
        except (ValueError, TypeError):
            logging.debug(f"Skipping batch percentile plot for {dataset_name}: invalid percentile values")
            return 'NA'
        
        # Create output directory for individual plots if it doesn't exist
        plots_dir = os.path.join(output_dir, 'batch_percentile_plots')
        os.makedirs(plots_dir, exist_ok=True)
        
        # Create safe filename from dataset name
        safe_name = "".join(c if c.isalnum() or c in ('_', '-') else '_' for c in dataset_name)
        output_path = os.path.join(plots_dir, f'{safe_name}_batch_percentile.html')
        
        # Create the interactive plot
        fig = go.Figure()
        
        # Add quality zone backgrounds with harmonious pastel colors
        # Poor zone (bottom-left: poor fit & poor geometry) - red
        fig.add_shape(type="rect", x0=0, y0=0, x1=50, y1=50, 
                      fillcolor="#FFB3B3", opacity=0.15, layer="below", line_width=0)
        
        # OK zones (yellow)
        fig.add_shape(type="rect", x0=50, y0=0, x1=75, y1=50, 
                      fillcolor="#FFF4C4", opacity=0.2, layer="below", line_width=0)
        fig.add_shape(type="rect", x0=0, y0=50, x1=50, y1=75, 
                      fillcolor="#FFF4C4", opacity=0.2, layer="below", line_width=0)
        fig.add_shape(type="rect", x0=50, y0=50, x1=75, y1=75, 
                      fillcolor="#FFF4C4", opacity=0.2, layer="below", line_width=0)
        
        # Good zones (light yellow/green transition)
        fig.add_shape(type="rect", x0=75, y0=0, x1=100, y1=50, 
                      fillcolor="#FFF4C4", opacity=0.15, layer="below", line_width=0)
        fig.add_shape(type="rect", x0=0, y0=75, x1=100, y1=100, 
                      fillcolor="#FFF4C4", opacity=0.15, layer="below", line_width=0)
        fig.add_shape(type="rect", x0=50, y0=75, x1=75, y1=100, 
                      fillcolor="#D4F1D4", opacity=0.2, layer="below", line_width=0)
        fig.add_shape(type="rect", x0=75, y0=50, x1=100, y1=75, 
                      fillcolor="#D4F1D4", opacity=0.2, layer="below", line_width=0)
        
        # Excellent zone (top-right: excellent fit & geometry) - green
        fig.add_shape(type="rect", x0=75, y0=75, x1=100, y1=100, 
                      fillcolor="#B8E6B8", opacity=0.25, layer="below", line_width=0)
        
        # Add reference lines
        fig.add_hline(y=75, line_dash="dash", line_color="darkgreen", opacity=0.6, line_width=2.5)
        fig.add_vline(x=75, line_dash="dash", line_color="darkgreen", opacity=0.6, line_width=2.5)
        fig.add_hline(y=50, line_dash="dash", line_color="orange", opacity=0.6, line_width=2.5)
        fig.add_vline(x=50, line_dash="dash", line_color="orange", opacity=0.6, line_width=2.5)
        
        # Add line labels as annotations positioned at edges
        fig.add_annotation(x=103, y=87.5, text="Top 25%<br>geometry", 
                          showarrow=False, xanchor='left', font=dict(size=12, color="darkgreen"))
        fig.add_annotation(x=87.5, y=103, text="Top 25% fit", 
                          showarrow=False, yanchor='bottom', font=dict(size=12, color="darkgreen"))
        fig.add_annotation(x=103, y=62.5, text="Median-<br>Top 25%<br>geometry", 
                          showarrow=False, xanchor='left', font=dict(size=12, color="orange"))
        fig.add_annotation(x=62.5, y=103, text="Median - Top 25% fit", 
                          showarrow=False, yanchor='bottom', font=dict(size=12, color="orange"))
        fig.add_annotation(x=103, y=25, text="Bottom<br>50%<br>geometry", 
                          showarrow=False, xanchor='left', font=dict(size=12, color="red"))
        fig.add_annotation(x=25, y=103, text="Bottom 50% fit", 
                          showarrow=False, yanchor='bottom', font=dict(size=12, color="red"))
        
        # Create hover text
        hover_text = f"<b style='font-size:14px'>{dataset_name}</b><br>"
        if ligand_id not in ['NA', None, '']:
            hover_text += f"<b style='font-size:12px'>Ligand: {ligand_id}</b><br>"
        hover_text += "─" * 30 + "<br>"
        hover_text += "<b>Batch Percentiles:</b><br>"
        hover_text += f"  • Fit Percentile: {fit_val:.1f}%<br>"
        hover_text += f"  • Geometry Percentile: {geom_val:.1f}%"
        
        # Plot the dataset point
        fig.add_trace(go.Scatter(
            x=[fit_val],
            y=[geom_val],
            mode='markers',
            marker=dict(
                size=15,
                color='#0d6efd',
                line=dict(color='black', width=2),
                opacity=0.9
            ),
            text=[hover_text],
            hovertemplate='%{text}<extra></extra>',
            name=dataset_name,
            showlegend=False
        ))
        
        # Update layout
        fig.update_layout(
            title=dict(
                text=f'<b>PCA Percentile Plot (Batch)</b><br><sub>{dataset_name}</sub>',
                x=0.5,
                xanchor='center',
                font=dict(size=18)
            ),
            xaxis=dict(
                title=dict(
                    text='<b>Ligand Fit PCA Percentile (Batch)</b><br><sub>Percentile rank within this batch</sub>',
                    font=dict(size=14)
                ),
                range=[-3, 110],
                gridcolor='lightgray',
                gridwidth=0.65,
                showgrid=True,
                tickmode='array',
                tickvals=[0, 25, 50, 75, 100],
                ticktext=['0', '25', '<b>50</b>', '<b>75</b>', '100'],
                tickfont=dict(size=12),
                zeroline=False
            ),
            yaxis=dict(
                title=dict(
                    text='<b>Ligand Geometry PCA Percentile (Batch)</b><br><sub>Percentile rank within this batch</sub>',
                    font=dict(size=14)
                ),
                range=[-3, 110],
                gridcolor='lightgray',
                gridwidth=0.65,
                showgrid=True,
                tickmode='array',
                tickvals=[0, 25, 50, 75, 100],
                ticktext=['0', '25', '<b>50</b>', '<b>75</b>', '100'],
                tickfont=dict(size=12),
                zeroline=False
            ),
            width=700,
            height=700,
            hovermode='closest',
            plot_bgcolor='white',
            margin=dict(l=80, r=80, t=100, b=80),
            hoverlabel=dict(
                bgcolor="white",
                font_size=13,
                font_family="Arial",
                bordercolor="black"
            )
        )
        
        # Add modebar config for consistent export options
        config = {
            'displayModeBar': True,
            'displaylogo': False,
            'modeBarButtonsToRemove': ['select2d', 'lasso2d'],
            'toImageButtonOptions': {
                'format': 'png',
                'filename': f'batch_percentile_{safe_name}',
                'height': 800,
                'width': 800,
                'scale': 2
            }
        }
        
        # Save as HTML
        import plotly.io as pio
        html_string = pio.to_html(fig, config=config, include_plotlyjs='cdn', full_html=True)
        html_string = html_string.replace('<head><meta charset="utf-8" /></head>',
                                         f'<head><meta charset="utf-8" /><title>PCA Percentile Plot (Batch) - {dataset_name}</title></head>')
        with open(output_path, 'w', encoding='utf-8') as f:
            f.write(html_string)
        
        logging.debug(f"Generated batch percentile plot for {dataset_name}: {output_path}")
        return output_path
        
    except Exception as e:
        logging.warning(f"Error generating batch percentile plot for {result.get('Crystal Name', 'Unknown')}: {e}")
        return 'NA'


def main() -> None:
    """Main entry point for collating Pipedream results."""
    parser = argparse.ArgumentParser(
        description='Collate Pipedream analysis results into interactive reports.',
        epilog='Example: python collate_pipedream_results.py --input results.json --output-dir /path/to/output --format both'
    )
