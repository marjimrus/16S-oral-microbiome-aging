#!/usr/bin/env python3
"""
Correlation Network Analysis for Oral Microbiome-Hormone Interactions
======================================================================

This script constructs and visualizes correlation networks between oral microbial
taxa and salivary estrogen levels (estradiol and estrone) in aging women.

Study: Salivary estrogens are associated with niche-specific oral microbiota in aging women

Methods:
    - Spearman's rank correlation coefficient for pairwise associations
    - Network construction with correlation threshold filtering
    - Force-directed graph layout algorithms (Kamada-Kawai, Fruchterman-Reingold)
    - Hub identification and network topology analysis

Input:
    - Microbial abundance data per oral niche (BM, TG, TH, GM)
    - Salivary hormone levels (estradiol, estrone)

Output:
    - Correlation matrices (Spearman, Pearson)
    - Network visualizations
    - CSV files with significant correlations

Requirements:
    - pandas >= 1.5.0
    - numpy >= 1.23.0
    - scipy >= 1.9.0
    - networkx >= 2.8.0
    - matplotlib >= 3.6.0

Author: Maria J. Rus
"""

import pandas as pd
import numpy as np
from scipy import stats
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import warnings

warnings.filterwarnings('ignore')


# =============================================================================
# CONFIGURATION
# =============================================================================

# Correlation thresholds
CORRELATION_THRESHOLD = 0.4      # Minimum correlation for network edges
STRONG_CORRELATION = 0.7         # Threshold for highlighting strong correlations
VERY_STRONG_CORRELATION = 0.9    # Threshold for very strong correlations

# Network visualization parameters
FIGURE_SIZE = (16, 12)
NODE_SIZE_MULTIPLIER = 100
SPRING_ITERATIONS = 50
SPRING_K = 0.5                   # Spring rest length

# Color scheme for oral niches
NICHE_COLORS = {
    'BM': '#E41A1C',  # Buccal Mucosa - Red
    'TG': '#377EB8',  # Tongue Dorsum - Blue
    'TH': '#4DAF4A',  # Tooth Surface (Supragingival) - Green
    'GM': '#984EA3',  # Gingival Margin (Subgingival) - Purple
    'hormone': '#FF7F00'  # Hormones - Orange
}


# =============================================================================
# DATA LOADING AND PREPROCESSING
# =============================================================================

def load_microbiome_data(filepath):
    """
    Load microbiome abundance data from CSV/Excel file.

    Expected format: rows = samples, columns = taxa
    First column should be sample IDs
    """
    if filepath.endswith('.xlsx'):
        data = pd.read_excel(filepath, index_col=0)
    else:
        data = pd.read_csv(filepath, index_col=0)
    return data


def prepare_correlation_data(abundance_data, hormone_data, niches=['BM', 'TG', 'TH', 'GM']):
    """
    Prepare data for correlation analysis by organizing taxa by niche.

    Parameters:
    -----------
    abundance_data : pd.DataFrame
        Microbial abundance data with taxa as columns
    hormone_data : pd.DataFrame
        Hormone levels (estradiol, estrone) per subject
    niches : list
        List of oral niche codes

    Returns:
    --------
    pd.DataFrame
        Combined dataset with niche-labeled taxa and hormone levels
    """
    combined_data = pd.DataFrame()

    for niche in niches:
        niche_cols = [col for col in abundance_data.columns if niche in col]
        for col in niche_cols:
            taxon_name = col.replace(f'_{niche}', '').replace(f'{niche}_', '')
            combined_data[f'{niche}_{taxon_name}'] = abundance_data[col]

    # Add hormone data
    combined_data['Estradiol'] = hormone_data['Estradiol']
    combined_data['Estrone'] = hormone_data['Estrone']

    return combined_data


# =============================================================================
# CORRELATION ANALYSIS
# =============================================================================

def calculate_correlations(data, method='spearman'):
    """
    Calculate pairwise correlations between all variables.

    Parameters:
    -----------
    data : pd.DataFrame
        Combined microbiome and hormone data
    method : str
        'spearman' or 'pearson'

    Returns:
    --------
    tuple
        (correlation_matrix, p_value_matrix)
    """
    n_vars = data.shape[1]
    corr_matrix = np.zeros((n_vars, n_vars))
    pval_matrix = np.zeros((n_vars, n_vars))

    for i in range(n_vars):
        for j in range(n_vars):
            if i == j:
                corr_matrix[i, j] = 1.0
                pval_matrix[i, j] = 0.0
            elif i < j:
                x = data.iloc[:, i].dropna()
                y = data.iloc[:, j].dropna()

                # Align indices
                common_idx = x.index.intersection(y.index)
                x = x.loc[common_idx]
                y = y.loc[common_idx]

                if len(common_idx) > 3:  # Minimum samples for correlation
                    if method == 'spearman':
                        corr, pval = stats.spearmanr(x, y)
                    else:
                        corr, pval = stats.pearsonr(x, y)

                    corr_matrix[i, j] = corr
                    corr_matrix[j, i] = corr
                    pval_matrix[i, j] = pval
                    pval_matrix[j, i] = pval

    corr_df = pd.DataFrame(corr_matrix,
                           index=data.columns,
                           columns=data.columns)
    pval_df = pd.DataFrame(pval_matrix,
                           index=data.columns,
                           columns=data.columns)

    return corr_df, pval_df


def filter_significant_correlations(corr_matrix, pval_matrix,
                                     corr_threshold=CORRELATION_THRESHOLD,
                                     pval_threshold=0.05):
    """
    Filter correlations by threshold and significance.

    Returns:
    --------
    pd.DataFrame
        Long-format dataframe with significant correlations
    """
    results = []

    for i, var1 in enumerate(corr_matrix.index):
        for j, var2 in enumerate(corr_matrix.columns):
            if i < j:  # Upper triangle only
                corr = corr_matrix.iloc[i, j]
                pval = pval_matrix.iloc[i, j]

                if abs(corr) >= corr_threshold and pval < pval_threshold:
                    results.append({
                        'Variable1': var1,
                        'Variable2': var2,
                        'Correlation': corr,
                        'P_value': pval,
                        'Abs_Correlation': abs(corr)
                    })

    return pd.DataFrame(results)


# =============================================================================
# NETWORK CONSTRUCTION
# =============================================================================

def build_correlation_network(correlations_df, corr_threshold=CORRELATION_THRESHOLD):
    """
    Build a NetworkX graph from correlation data.

    Parameters:
    -----------
    correlations_df : pd.DataFrame
        Filtered correlations with Variable1, Variable2, Correlation columns
    corr_threshold : float
        Minimum absolute correlation for edge creation

    Returns:
    --------
    nx.Graph
        Correlation network
    """
    G = nx.Graph()

    # Add all unique variables as nodes
    all_vars = set(correlations_df['Variable1']).union(set(correlations_df['Variable2']))

    for var in all_vars:
        # Determine node type (niche or hormone)
        if var in ['Estradiol', 'Estrone']:
            niche = 'hormone'
        else:
            niche = var.split('_')[0] if '_' in var else 'unknown'

        G.add_node(var, niche=niche)

    # Add edges with correlation weights
    for _, row in correlations_df.iterrows():
        if abs(row['Correlation']) >= corr_threshold:
            G.add_edge(
                row['Variable1'],
                row['Variable2'],
                weight=row['Correlation'],
                abs_weight=abs(row['Correlation']),
                direction='positive' if row['Correlation'] > 0 else 'negative'
            )

    return G


def calculate_network_metrics(G):
    """
    Calculate network topology metrics.

    Returns:
    --------
    dict
        Network statistics
    """
    metrics = {
        'n_nodes': G.number_of_nodes(),
        'n_edges': G.number_of_edges(),
        'density': nx.density(G),
        'n_components': nx.number_connected_components(G),
    }

    if G.number_of_nodes() > 0:
        metrics['avg_degree'] = sum(dict(G.degree()).values()) / G.number_of_nodes()

        # Degree centrality
        degree_cent = nx.degree_centrality(G)
        metrics['top_hubs'] = sorted(degree_cent.items(),
                                      key=lambda x: x[1],
                                      reverse=True)[:10]

        # Betweenness centrality
        if nx.is_connected(G):
            between_cent = nx.betweenness_centrality(G)
            metrics['top_bridges'] = sorted(between_cent.items(),
                                            key=lambda x: x[1],
                                            reverse=True)[:10]

    return metrics


# =============================================================================
# VISUALIZATION
# =============================================================================

def visualize_network(G, title="Microbiome-Hormone Correlation Network",
                      layout='spring', save_path=None):
    """
    Visualize the correlation network.

    Parameters:
    -----------
    G : nx.Graph
        Correlation network
    title : str
        Plot title
    layout : str
        'spring' (Fruchterman-Reingold) or 'kamada_kawai'
    save_path : str
        Path to save figure (optional)
    """
    fig, ax = plt.subplots(figsize=FIGURE_SIZE)

    # Calculate layout
    if layout == 'kamada_kawai':
        pos = nx.kamada_kawai_layout(G)
    else:
        pos = nx.spring_layout(G,
                               k=SPRING_K,
                               iterations=SPRING_ITERATIONS,
                               seed=42)

    # Node colors by niche
    node_colors = []
    for node in G.nodes():
        niche = G.nodes[node].get('niche', 'unknown')
        node_colors.append(NICHE_COLORS.get(niche, '#CCCCCC'))

    # Node sizes by degree
    degrees = dict(G.degree())
    node_sizes = [degrees[node] * NODE_SIZE_MULTIPLIER + 50 for node in G.nodes()]

    # Edge colors by correlation direction
    edge_colors = []
    edge_widths = []
    for u, v in G.edges():
        corr = G[u][v]['weight']
        if corr > STRONG_CORRELATION:
            edge_colors.append('#0000FF')  # Strong positive - Blue
            edge_widths.append(2.5)
        elif corr > 0:
            edge_colors.append('#ADD8E6')  # Weak positive - Light blue
            edge_widths.append(1.0)
        elif corr < -STRONG_CORRELATION:
            edge_colors.append('#FF0000')  # Strong negative - Red
            edge_widths.append(2.5)
        else:
            edge_colors.append('#FFB6C1')  # Weak negative - Light red
            edge_widths.append(1.0)

    # Draw network
    nx.draw_networkx_edges(G, pos,
                           edge_color=edge_colors,
                           width=edge_widths,
                           alpha=0.6,
                           ax=ax)

    nx.draw_networkx_nodes(G, pos,
                           node_color=node_colors,
                           node_size=node_sizes,
                           alpha=0.8,
                           ax=ax)

    # Labels for hormone nodes and high-degree nodes
    labels = {}
    for node in G.nodes():
        if G.nodes[node].get('niche') == 'hormone' or degrees[node] >= 5:
            labels[node] = node.split('_')[-1] if '_' in node else node

    nx.draw_networkx_labels(G, pos, labels,
                            font_size=8,
                            font_weight='bold',
                            ax=ax)

    # Legend
    legend_elements = [
        plt.scatter([], [], c=NICHE_COLORS['BM'], s=100, label='Buccal Mucosa (BM)'),
        plt.scatter([], [], c=NICHE_COLORS['TG'], s=100, label='Tongue Dorsum (TG)'),
        plt.scatter([], [], c=NICHE_COLORS['TH'], s=100, label='Supragingival (TH)'),
        plt.scatter([], [], c=NICHE_COLORS['GM'], s=100, label='Subgingival (GM)'),
        plt.scatter([], [], c=NICHE_COLORS['hormone'], s=100, label='Hormones'),
    ]
    ax.legend(handles=legend_elements, loc='upper left', fontsize=10)

    ax.set_title(title, fontsize=14, fontweight='bold')
    ax.axis('off')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved to: {save_path}")

    plt.show()

    return fig


def visualize_hormone_subnetwork(G, hormone='Estrone', save_path=None):
    """
    Visualize correlations for a specific hormone.

    Parameters:
    -----------
    G : nx.Graph
        Full correlation network
    hormone : str
        'Estradiol' or 'Estrone'
    save_path : str
        Path to save figure (optional)
    """
    # Extract subgraph with hormone and its neighbors
    if hormone not in G.nodes():
        print(f"Hormone {hormone} not found in network")
        return

    neighbors = list(G.neighbors(hormone))
    subgraph_nodes = [hormone] + neighbors
    subG = G.subgraph(subgraph_nodes).copy()

    # Organize by niche
    niche_positions = {'BM': (-1, 0), 'TG': (1, 0), 'TH': (0, 1), 'GM': (0, -1)}
    pos = {hormone: (0, 0)}

    for niche, (base_x, base_y) in niche_positions.items():
        niche_nodes = [n for n in neighbors if n.startswith(niche)]
        for i, node in enumerate(niche_nodes):
            angle = (i / max(len(niche_nodes), 1)) * np.pi / 2 - np.pi / 4
            radius = 1.5
            x = base_x * radius + np.cos(angle) * 0.5
            y = base_y * radius + np.sin(angle) * 0.5
            pos[node] = (x, y)

    fig, ax = plt.subplots(figsize=(12, 10))

    # Node colors
    node_colors = []
    for node in subG.nodes():
        niche = subG.nodes[node].get('niche', 'unknown')
        node_colors.append(NICHE_COLORS.get(niche, '#CCCCCC'))

    # Edge colors by correlation
    edge_colors = ['#0000FF' if subG[u][v]['weight'] > 0 else '#FF0000'
                   for u, v in subG.edges()]
    edge_widths = [abs(subG[u][v]['weight']) * 3 for u, v in subG.edges()]

    nx.draw_networkx_edges(subG, pos, edge_color=edge_colors,
                           width=edge_widths, alpha=0.7, ax=ax)
    nx.draw_networkx_nodes(subG, pos, node_color=node_colors,
                           node_size=500, alpha=0.9, ax=ax)

    # Labels
    labels = {node: node.split('_')[-1] if '_' in node else node
              for node in subG.nodes()}
    nx.draw_networkx_labels(subG, pos, labels, font_size=8, ax=ax)

    ax.set_title(f'{hormone} Correlation Network\n'
                 f'(Blue = positive, Red = negative)',
                 fontsize=14, fontweight='bold')
    ax.axis('off')

    plt.tight_layout()

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')

    plt.show()

    return fig


# =============================================================================
# MAIN ANALYSIS PIPELINE
# =============================================================================

def run_network_analysis(abundance_file, hormone_file, output_dir='results'):
    """
    Run complete network analysis pipeline.

    Parameters:
    -----------
    abundance_file : str
        Path to microbiome abundance data
    hormone_file : str
        Path to hormone level data
    output_dir : str
        Directory for output files
    """
    import os
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print("Microbiome-Hormone Correlation Network Analysis")
    print("=" * 60)

    # Load data
    print("\n1. Loading data...")
    abundance_data = load_microbiome_data(abundance_file)
    hormone_data = load_microbiome_data(hormone_file)
    print(f"   Loaded {abundance_data.shape[1]} microbial variables")
    print(f"   Loaded hormone data for {len(hormone_data)} subjects")

    # Prepare combined dataset
    print("\n2. Preparing correlation data...")
    # Note: This step may need adjustment based on actual data format
    combined_data = abundance_data.copy()
    combined_data['Estradiol'] = hormone_data.get('Estradiol', hormone_data.iloc[:, 0])
    combined_data['Estrone'] = hormone_data.get('Estrone', hormone_data.iloc[:, 1])
    print(f"   Combined dataset: {combined_data.shape[0]} samples, {combined_data.shape[1]} variables")

    # Calculate correlations
    print("\n3. Calculating correlations...")

    # Spearman correlations
    spearman_corr, spearman_pval = calculate_correlations(combined_data, method='spearman')
    spearman_sig = filter_significant_correlations(spearman_corr, spearman_pval,
                                                    corr_threshold=CORRELATION_THRESHOLD)
    print(f"   Spearman: {len(spearman_sig)} significant correlations (|r| >= {CORRELATION_THRESHOLD})")

    # Pearson correlations
    pearson_corr, pearson_pval = calculate_correlations(combined_data, method='pearson')
    pearson_sig = filter_significant_correlations(pearson_corr, pearson_pval,
                                                   corr_threshold=CORRELATION_THRESHOLD)
    print(f"   Pearson: {len(pearson_sig)} significant correlations (|r| >= {CORRELATION_THRESHOLD})")

    # Save correlation results
    spearman_sig.to_csv(f'{output_dir}/spearman_correlations_{CORRELATION_THRESHOLD}.csv', index=False)
    pearson_sig.to_csv(f'{output_dir}/pearson_correlations_{CORRELATION_THRESHOLD}.csv', index=False)
    print(f"   Saved correlation tables to {output_dir}/")

    # Build network
    print("\n4. Building correlation network...")
    G = build_correlation_network(spearman_sig)
    metrics = calculate_network_metrics(G)

    print(f"   Nodes: {metrics['n_nodes']}")
    print(f"   Edges: {metrics['n_edges']}")
    print(f"   Density: {metrics['density']:.4f}")
    print(f"   Average degree: {metrics.get('avg_degree', 0):.2f}")

    # Strong correlations
    strong_corr = spearman_sig[spearman_sig['Abs_Correlation'] >= STRONG_CORRELATION]
    print(f"   Strong correlations (|r| >= {STRONG_CORRELATION}): {len(strong_corr)}")

    # Save network data
    nx.write_graphml(G, f'{output_dir}/correlation_network.graphml')

    # Visualize
    print("\n5. Generating visualizations...")
    visualize_network(G,
                      title="Oral Microbiome-Hormone Correlation Network",
                      layout='spring',
                      save_path=f'{output_dir}/network_full.png')

    visualize_hormone_subnetwork(G, hormone='Estrone',
                                  save_path=f'{output_dir}/network_estrone.png')

    visualize_hormone_subnetwork(G, hormone='Estradiol',
                                  save_path=f'{output_dir}/network_estradiol.png')

    print("\n" + "=" * 60)
    print("Analysis complete!")
    print("=" * 60)

    return G, spearman_sig, metrics


# =============================================================================
# COMMAND LINE INTERFACE
# =============================================================================

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description='Correlation network analysis for microbiome-hormone data'
    )
    parser.add_argument('--abundance', '-a', required=True,
                        help='Path to microbiome abundance file')
    parser.add_argument('--hormones', '-h', required=True,
                        help='Path to hormone data file')
    parser.add_argument('--output', '-o', default='results',
                        help='Output directory')
    parser.add_argument('--threshold', '-t', type=float, default=0.4,
                        help='Correlation threshold (default: 0.4)')

    args = parser.parse_args()

    CORRELATION_THRESHOLD = args.threshold

    run_network_analysis(
        abundance_file=args.abundance,
        hormone_file=args.hormones,
        output_dir=args.output
    )
