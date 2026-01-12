
"""
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap
from pathlib import Path
import os

# Set working directory
script_dir = Path(__file__).parent
project_root = script_dir.parent
os.chdir(project_root)

print(f"Working directory: {os.getcwd()}\n")

# Define pathways with enzyme KO numbers organized by step
CARBON_PATHWAYS = {
    'Classical EMP\n(Standard Glycolysis)': {
        '1': ['K00844'], '2': ['K01810'], '3': ['K00850'], '4': ['K01623', 'K01624'],
        '6': ['K00134', 'K00150'], '7': ['K00927'], '8': ['K01834'],
        '9': ['K15633', 'K15634'], '10': ['K00873', 'K12406']
    },
    'Archaeal Modified EMP\n(ADP-dependent)': {
        '1': ['K00845'], '2': ['K01810'], '3': ['K00918'], '4': ['K01622', 'K01623', 'K01624'],
        '6': ['K00134', 'K10705'], '7': ['K00927'], '8': ['K01689'],
        '9': ['K15633', 'K15634'], '10': ['K00873', 'K12406']
    },
    'Entner-Doudoroff\n(Semi-phosphorylative)': {
        '1': ['K00844', 'K00845'], '2': ['K01810'], '3': ['K00036'], '4': ['K01057'],
        '5': ['K01625'], '6': ['K00134', 'K00150'], '7': ['K00927'], '8': ['K01834'],
        '9': ['K15633'], '10': ['K00873', 'K12406']
    },
    'Lower Glycolysis\n(GAP → Pyruvate)': {
        '5': ['K01803'], '6': ['K00134', 'K00150', 'K10705'], '7': ['K00927'],
        '8': ['K01834', 'K01689'], '9': ['K15633', 'K15634'], '10': ['K00873', 'K12406']
    },
    'Glycerol Metabolism\n(Glycerol → DHAP)': {
        '1': ['K00864'], '2': ['K00111', 'K00112', 'K00865'], 'interchange': ['K01803']
    }
}

def load_ghostkoala_file(filepath):
    """Load GhostKOALA annotation file and extract KO numbers"""
    kos = set()
    with open(filepath) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2 and parts[1]:
                kos.add(parts[1])
    return kos

def calc_completeness(strain_kos, pathway_steps):
    """Calculate pathway completeness: (steps covered / total steps) * 100"""
    steps_covered = sum(1 for ko_list in pathway_steps.values() 
                       if any(ko in strain_kos for ko in ko_list))
    return round((steps_covered / len(pathway_steps)) * 100)

def main():
    print("="*60)
    print("LOADING GHOSTKOALA DATA")
    print("="*60)
    
    # Load GhostKOALA files
    ghostkoala_dir = Path('annotation_outputs/GhostKOALA')
    strain_kos = {}
    
    for strain in ['HS', 'HM', 'HA']:
        filepath = ghostkoala_dir / f'{strain}_GhostKOALA.txt'
        if filepath.exists():
            strain_kos[strain] = load_ghostkoala_file(filepath)
            print(f"Loaded {strain}: {len(strain_kos[strain])} KOs")
        else:
            print(f"ERROR: {filepath} not found!")
            return
    
    print("\n" + "="*60)
    print("CALCULATING PATHWAY COMPLETENESS")
    print("="*60)
    
    # Calculate completeness for each pathway and strain
    results = []
    for pathway_name, pathway_steps in CARBON_PATHWAYS.items():
        pathway_short = pathway_name.split('\n')[0]
        print(f"\n{pathway_short}:")
        
        for strain in ['HS', 'HM', 'HA']:
            comp = calc_completeness(strain_kos[strain], pathway_steps)
            results.append({
                'Pathway': pathway_name,
                'Strain': strain,
                'Completeness': comp
            })
            print(f"  {strain}: {comp:3d}%")
    
    df = pd.DataFrame(results)
    
    # Pivot for visualization
    df_pivot = df.pivot(index='Pathway', columns='Strain', values='Completeness')
    df_pivot = df_pivot[['HS', 'HM', 'HA']]  # Ensure order
    
    print("\n" + "="*60)
    print("CREATING FIGURES")
    print("="*60)
    
    Path('results/figures').mkdir(parents=True, exist_ok=True)
    Path('results/tables').mkdir(parents=True, exist_ok=True)
    
    # Custom colormap
    colors_map = ['#FFFFFF', '#FFF59D', '#FFD54F', '#FFA726', '#FF6F00', '#D84315', '#B71C1C']
    cmap = LinearSegmentedColormap.from_list('custom', colors_map, N=100)
    
    # ========== BAR PLOT ==========
    fig, ax = plt.subplots(figsize=(12, 8))
    
    pathways_ordered = list(CARBON_PATHWAYS.keys())
    x = np.arange(len(pathways_ordered))
    width = 0.25
    
    colors = ['#FF6B6B', '#4ECDC4', '#95E1D3']
    labels = ['HS (H. salinarum)', 'HM (Halomicrobium sp.)', 'HA (Haloarcula sp.)']
    
    for i, (strain, color, label) in enumerate(zip(['HS', 'HM', 'HA'], colors, labels)):
        values = [df[(df['Pathway']==p) & (df['Strain']==strain)]['Completeness'].values[0] 
                 for p in pathways_ordered]
        bars = ax.bar(x + i*width, values, width, label=label, 
                     color=color, edgecolor='white', linewidth=1.5)
        
        # Add value labels
        for bar, val in zip(bars, values):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 1,
                   f'{val}%', ha='center', va='bottom', 
                   fontsize=10, fontweight='bold')
    
    # Reference lines
    ax.axhline(y=100, color='#E74C3C', linestyle='--', alpha=0.5, linewidth=1.5)
    ax.axhline(y=70, color='#F39C12', linestyle='--', alpha=0.3, linewidth=1.5)
    
    ax.set_ylabel('Pathway Completeness (%)', fontsize=12, fontweight='bold')
    ax.set_title('Carbon Metabolism Pathway Completeness', 
                 fontsize=14, fontweight='bold', pad=15)
    ax.set_xticks(x + width)
    ax.set_xticklabels(pathways_ordered, fontsize=10)
    ax.set_ylim(0, 115)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.10), ncol=3, 
              fontsize=10, frameon=True)
    
    plt.tight_layout()
    plt.subplots_adjust(bottom=0.18)
    
    barplot_file = 'results/figures/glucose_pathways_completeness_barplot.png'
    plt.savefig(barplot_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✓ {barplot_file}")
    
    # ========== HEATMAP ==========
    fig, ax = plt.subplots(figsize=(8, 7))
    
    sns.heatmap(df_pivot, annot=True, fmt='d', cmap=cmap,
                linewidths=2, linecolor='white',
                cbar_kws={'label': 'Completeness (%)', 'shrink': 0.8},
                annot_kws={'size': 14, 'weight': 'bold'},
                vmin=0, vmax=100, ax=ax)
    
    ax.set_xticklabels(['HS\n(H. salinarum)', 'HM\n(Halomicrobium sp.)', 'HA\n(Haloarcula sp.)'],
                       fontsize=11, fontweight='bold')
    ax.set_yticklabels(ax.get_yticklabels(), fontsize=10, rotation=0)
    
    ax.set_title('Carbon Metabolism Pathway Completeness (%)', 
                 fontsize=14, fontweight='bold', pad=15)
    
    plt.tight_layout()
    
    heatmap_file = 'results/figures/glucose_pathways_completeness_heatmap.png'
    plt.savefig(heatmap_file, dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✓ {heatmap_file}")
    
    # Save summary table
    table_file = 'results/tables/glucose_pathways_completeness_summary.csv'
    df.to_csv(table_file, index=False)
    print(f"✓ {table_file}")
    
    print("\n" + "="*60)
    print("COMPLETE! Numbers calculated from actual GhostKOALA data")
    print("="*60)

if __name__ == "__main__":
    main()