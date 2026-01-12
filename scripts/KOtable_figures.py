
"""

"""

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path
import os

# Set working directory to project root
script_dir = Path(__file__).parent
project_root = script_dir.parent if script_dir.name == 'scripts' else script_dir
os.chdir(project_root)

print(f"Working directory: {os.getcwd()}\n")

# Output directory
OUTPUT_DIR = Path('results/figures')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# =============================================================================
# DATA ORGANIZED BY BIOLOGICAL FUNCTION
# =============================================================================

# Format: (KO, gene, protein, HS, HM, HA, distribution_category)
functional_data = {
    "Amino Acid Biosynthesis": [
        ("K01652", "ilvB", "Acetolactate synthase large subunit", False, True, True, "HM–HA"),
        ("K01653", "ilvH", "Acetolactate synthase small subunit", False, True, True, "HM–HA"),
        ("K00053", "ilvC", "Ketol-acid reductoisomerase", False, True, True, "HM–HA"),
        ("K01687", "ilvD", "Dihydroxy-acid dehydratase (BCAA)", False, True, True, "HM–HA"),
        ("K01649", "leuA", "2-Isopropylmalate synthase", False, True, True, "HM–HA"),
        ("K01703", "leuC", "3-Isopropylmalate dehydratase large", False, True, True, "HM–HA"),
        ("K01704", "leuD", "3-Isopropylmalate dehydratase small", False, True, True, "HM–HA"),
        ("K00931", "proB", "Glutamate 5-kinase (proline)", False, True, False, "HM unique"),
        ("K00147", "proA", "Glutamate-5-semialdehyde dehydrogenase", False, True, False, "HM unique"),
        ("K00286", "proC", "Pyrroline-5-carboxylate reductase", False, True, False, "HM unique"),
        ("K02510", "hisG", "ATP phosphoribosyltransferase (His)", False, True, False, "HM unique"),
        ("K00215", "dapB", "Dihydrodipicolinate reductase (Lys)", False, True, True, "HM–HA"),
        ("K01586", "lysA", "Diaminopimelate decarboxylase (Lys)", False, True, True, "HM–HA"),
        ("K05826", "lysW", "Lysine biosynthesis carrier protein", False, True, True, "HM–HA"),
    ],
    
    "Amino Acid Catabolism": [
        ("K01745", "hutH", "Histidine ammonia-lyase", True, False, True, "HS–HA"),
        ("K01712", "hutU", "Urocanate hydratase", True, False, True, "HS–HA"),
        ("K23316", "hutI", "Imidazolonepropionase", True, False, True, "HS–HA"),
        ("K23684", "hutG", "Formiminoglutamate hydrolase", True, False, True, "HS–HA"),
        ("K01667", "tnaA", "Tryptophanase", True, False, True, "HS–HA"),
        ("K08177", "tnaB", "Tryptophan permease", True, False, False, "HS unique"),
        ("K00926", "arcC", "Carbamate kinase (Arg)", True, False, False, "HS unique"),
        ("K01620", "ltaE", "Threonine aldolase", True, True, False, "HS–HM"),
    ],
    
    "Carbon Metabolism": [
        ("K00005", "gldA", "Glycerol dehydrogenase", True, False, False, "HS unique"),
        ("K00162", "pdhB", "Pyruvate dehydrogenase E1 beta", True, False, True, "HS–HA"),
        ("K00265", "gltB", "Glutamate synthase large subunit", False, True, True, "HM–HA"),
        ("K01601", "rbcL", "RuBisCO large subunit (CO₂ fixation)", False, True, False, "HM unique"),
        ("K00697", "glgA", "Glycogen synthase", False, False, True, "HA unique"),
        ("K01178", "amyA", "Alpha-amylase", False, False, True, "HA unique"),
        ("K01191", "treA", "Trehalase", False, True, False, "HM unique"),
        ("K01192", "treF", "Trehalase", False, True, False, "HM unique"),
    ],
    
    "Nitrogen Metabolism": [
        ("K01428", "ureA", "Urease gamma subunit", False, False, True, "HA unique"),
        ("K01429", "ureB", "Urease beta subunit", False, False, True, "HA unique"),
        ("K01430", "ureC", "Urease alpha subunit", False, False, True, "HA unique"),
        ("K11959", "urtA", "Urea transporter substrate-binding", False, False, True, "HA unique"),
        ("K11960", "urtB", "Urea transporter permease", False, False, True, "HA unique"),
        ("K00376", "nosZ", "Nitrous oxide reductase", False, False, True, "HA unique"),
        ("K17316", "norB", "Nitric oxide reductase subunit B", False, False, True, "HA unique"),
    ],
    
    "Transport Systems": [
        ("K02784", "crr", "PTS glucose EIIA component", False, False, True, "HA unique"),
        ("K02770", "ptsG", "PTS glucose EIICB component", False, False, True, "HA unique"),
        ("K14445", "slc13a5", "Citrate transporter", True, False, True, "HS–HA"),
        ("K02041", "phnC", "Phosphonate transporter ATP-binding", True, False, True, "HS–HA"),
        ("K01995", "livH", "BCAA transporter permease", False, True, True, "HM–HA"),
        ("K01999", "livK", "BCAA transporter substrate-binding", False, True, True, "HM–HA"),
        ("K02052", "potA", "Polyamine transporter ATP-binding", True, True, False, "HS–HM"),
        ("K02053", "potB", "Polyamine transporter permease", True, True, False, "HS–HM"),
    ],
    
    "Energy & Phototrophy": [
        ("K03892", "bop", "Bacteriorhodopsin", True, True, False, "HS–HM"),
        ("K03469", "bat", "Bacterio-opsin activator", True, False, False, "HS unique"),
        ("K01546", "ntpA", "V-type ATPase subunit A", True, False, False, "HS unique"),
        ("K01547", "ntpB", "V-type ATPase subunit B", True, False, False, "HS unique"),
    ],
    
    "Other Functions": [
        ("K01183", "chiA", "Chitinase", True, True, False, "HS–HM"),
        ("K06901", "phrB", "DNA photolyase", True, True, False, "HS–HM"),
        ("K01495", "folE", "GTP cyclohydrolase I (folate)", False, True, False, "HM unique"),
    ],
}

# Colors matching the distribution categories
distribution_colors = {
    "HS unique": "#FADBD8",
    "HM unique": "#D5F5E3", 
    "HA unique": "#D6EAF8",
    "HS–HM": "#F5EEF8",
    "HS–HA": "#FCF3CF",
    "HM–HA": "#D1F2EB",
}

# Strain colors for dots
strain_colors = {
    "HS": "#E74C3C",
    "HM": "#27AE60",
    "HA": "#3498DB",
}

# =============================================================================
# FIGURE A: FUNCTIONAL SUMMARY TABLE
# =============================================================================

def create_functional_table():
    """Create table organized by biological function with strain presence dots."""
    
    total_gene_rows = sum(len(genes) for genes in functional_data.values())
    total_category_rows = len(functional_data)
    fig_height = 2.5 + (total_gene_rows + total_category_rows) * 0.28
    
    fig, ax = plt.subplots(figsize=(11, min(fig_height, 18)))
    ax.axis('off')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    
    headers = ['KO ID', 'Gene', 'Protein/Enzyme', 'HS', 'HM', 'HA']
    col_widths = [0.10, 0.08, 0.58, 0.06, 0.06, 0.06]
    col_positions = [0.03]
    for w in col_widths[:-1]:
        col_positions.append(col_positions[-1] + w)
    
    row_height = 0.0135
    y_start = 0.94
    
    # Title
    ax.text(0.5, 0.98, 'Functional Classification of Key KEGG Orthologs',
            ha='center', va='top', fontsize=13, fontweight='bold')
    ax.text(0.5, 0.96, 'HS: Halobacterium salinarum  |  HM: Halomicrobium kobeituziens  |  HA: Haloarcula sp.',
            ha='center', va='top', fontsize=9, style='italic', color='#555555')
    
    # Header row
    header_y = y_start
    for header, x, w in zip(headers, col_positions, col_widths):
        ax.add_patch(plt.Rectangle((x, header_y - row_height), w - 0.005, row_height,
                                   facecolor='#2C3E50', edgecolor='white', linewidth=0.5))
        ax.text(x + w/2 - 0.0025, header_y - row_height/2, header,
                ha='center', va='center', fontsize=9, fontweight='bold', color='white')
    
    current_y = header_y - row_height
    
    # Process each functional category
    for category, genes in functional_data.items():
        current_y -= row_height * 1.3
        total_width = sum(col_widths) - 0.005
        ax.add_patch(plt.Rectangle((0.03, current_y), total_width, row_height * 1.15,
                                   facecolor='#34495E', edgecolor='white', linewidth=0.5))
        ax.text(0.03 + total_width/2, current_y + row_height * 0.575, category,
                ha='center', va='center', fontsize=10, fontweight='bold', color='white')
        
        for ko, gene, protein, hs, hm, ha, dist in genes:
            current_y -= row_height
            bg_color = distribution_colors.get(dist, 'white')
            
            # Draw text cells
            row_data = [ko, gene, protein]
            for i, (val, x, w) in enumerate(zip(row_data, col_positions[:3], col_widths[:3])):
                ax.add_patch(plt.Rectangle((x, current_y), w - 0.005, row_height,
                                           facecolor=bg_color, edgecolor='#BDC3C7', linewidth=0.3))
                
                fontsize = 8
                fontstyle = 'italic' if i == 1 else 'normal'
                ha_align = 'center' if i in [0, 1] else 'left'
                x_off = w/2 - 0.0025 if i in [0, 1] else 0.006
                
                ax.text(x + x_off, current_y + row_height/2, val,
                        ha=ha_align, va='center', fontsize=fontsize, fontstyle=fontstyle)
            
            # Draw presence dots
            for j, (present, strain) in enumerate([(hs, 'HS'), (hm, 'HM'), (ha, 'HA')]):
                x = col_positions[3 + j]
                w = col_widths[3 + j]
                ax.add_patch(plt.Rectangle((x, current_y), w - 0.005, row_height,
                                           facecolor=bg_color, edgecolor='#BDC3C7', linewidth=0.3))
                if present:
                    ax.scatter(x + w/2 - 0.0025, current_y + row_height/2, 
                              s=50, color=strain_colors[strain], zorder=5)
    
    # Legend
    legend_y = current_y - 0.035
    ax.text(0.03, legend_y + 0.012, 'Distribution Categories:', fontsize=9, fontweight='bold')
    
    legend_items_row1 = [
        ("HS unique (71)", "#FADBD8"),
        ("HM unique (67)", "#D5F5E3"),
        ("HA unique (129)", "#D6EAF8"),
    ]
    legend_items_row2 = [
        ("HS–HM (22)", "#F5EEF8"),
        ("HS–HA (41)", "#FCF3CF"),
        ("HM–HA (133)", "#D1F2EB"),
    ]
    
    legend_x = 0.03
    for label, color in legend_items_row1:
        ax.add_patch(plt.Rectangle((legend_x, legend_y - 0.015), 0.018, 0.012,
                                   facecolor=color, edgecolor='black', linewidth=0.5))
        ax.text(legend_x + 0.022, legend_y - 0.009, label, ha='left', va='center', fontsize=8)
        legend_x += 0.18
    
    legend_y -= 0.022
    legend_x = 0.03
    for label, color in legend_items_row2:
        ax.add_patch(plt.Rectangle((legend_x, legend_y - 0.015), 0.018, 0.012,
                                   facecolor=color, edgecolor='black', linewidth=0.5))
        ax.text(legend_x + 0.022, legend_y - 0.009, label, ha='left', va='center', fontsize=8)
        legend_x += 0.18
    
    ax.text(0.60, legend_y - 0.009, 'Strain:', fontsize=9, fontweight='bold', va='center')
    dot_x = 0.68
    for strain, color in strain_colors.items():
        ax.scatter(dot_x, legend_y - 0.009, s=50, color=color, zorder=5)
        ax.text(dot_x + 0.015, legend_y - 0.009, strain, ha='left', va='center', fontsize=8)
        dot_x += 0.08
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'functionalKOtable.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(OUTPUT_DIR / 'functionalKOtable.pdf', bbox_inches='tight', facecolor='white')
    plt.close()
    print("✓ Figure A (Functional): functionalKOtable.png + .pdf")


# =============================================================================
# FIGURE B: DISTRIBUTION-BASED TABLE
# =============================================================================

distribution_data = [
    # HS unique
    ("HS unique", "K00926", "arcC", "Carbamate kinase", "Amino acid catabolism"),
    ("HS unique", "K08177", "tnaB", "Tryptophan permease", "Amino acid catabolism"),
    ("HS unique", "K00005", "gldA", "Glycerol dehydrogenase", "Carbon metabolism"),
    ("HS unique", "K01546", "ntpA", "V-type ATPase subunit A", "Energy"),
    ("HS unique", "K01547", "ntpB", "V-type ATPase subunit B", "Energy"),
    ("HS unique", "K03469", "bat", "Bacterio-opsin activator", "Phototrophy"),
    
    # HM unique
    ("HM unique", "K01601", "rbcL", "RuBisCO large subunit", "Carbon fixation"),
    ("HM unique", "K00931", "proB", "Glutamate 5-kinase", "AA biosynthesis (Pro)"),
    ("HM unique", "K00147", "proA", "Glutamate-5-semialdehyde dehydrogenase", "AA biosynthesis (Pro)"),
    ("HM unique", "K00286", "proC", "Pyrroline-5-carboxylate reductase", "AA biosynthesis (Pro)"),
    ("HM unique", "K02510", "hisG", "ATP phosphoribosyltransferase", "AA biosynthesis (His)"),
    ("HM unique", "K01191", "treA", "Trehalase", "Carbon metabolism"),
    ("HM unique", "K01192", "treF", "Trehalase", "Carbon metabolism"),
    ("HM unique", "K01495", "folE", "GTP cyclohydrolase I", "Cofactor biosynthesis"),
    
    # HA unique
    ("HA unique", "K01428", "ureA", "Urease gamma subunit", "Nitrogen metabolism"),
    ("HA unique", "K01429", "ureB", "Urease beta subunit", "Nitrogen metabolism"),
    ("HA unique", "K01430", "ureC", "Urease alpha subunit", "Nitrogen metabolism"),
    ("HA unique", "K11959", "urtA", "Urea transporter substrate-binding", "Transport"),
    ("HA unique", "K11960", "urtB", "Urea transporter permease", "Transport"),
    ("HA unique", "K02784", "crr", "PTS glucose EIIA component", "Transport (PTS)"),
    ("HA unique", "K02770", "ptsG", "PTS glucose EIICB component", "Transport (PTS)"),
    ("HA unique", "K00376", "nosZ", "Nitrous oxide reductase", "Denitrification"),
    ("HA unique", "K17316", "norB", "Nitric oxide reductase subunit B", "Denitrification"),
    ("HA unique", "K00697", "glgA", "Glycogen synthase", "Carbon storage"),
    ("HA unique", "K01178", "amyA", "Alpha-amylase", "Carbon metabolism"),
    
    # HS-HM shared
    ("HS–HM shared", "K03892", "bop", "Bacteriorhodopsin", "Phototrophy"),
    ("HS–HM shared", "K02052", "potA", "Polyamine transporter ATP-binding", "Transport"),
    ("HS–HM shared", "K02053", "potB", "Polyamine transporter permease", "Transport"),
    ("HS–HM shared", "K01620", "ltaE", "Threonine aldolase", "AA catabolism"),
    ("HS–HM shared", "K01183", "chiA", "Chitinase", "Polysaccharide degradation"),
    ("HS–HM shared", "K06901", "phrB", "DNA photolyase", "DNA repair"),
    
    # HS-HA shared
    ("HS–HA shared", "K01745", "hutH", "Histidine ammonia-lyase", "AA catabolism (His)"),
    ("HS–HA shared", "K01712", "hutU", "Urocanate hydratase", "AA catabolism (His)"),
    ("HS–HA shared", "K23316", "hutI", "Imidazolonepropionase", "AA catabolism (His)"),
    ("HS–HA shared", "K23684", "hutG", "Formiminoglutamate hydrolase", "AA catabolism (His)"),
    ("HS–HA shared", "K01667", "tnaA", "Tryptophanase", "AA catabolism (Trp)"),
    ("HS–HA shared", "K00162", "pdhB", "Pyruvate dehydrogenase E1 beta", "Central metabolism"),
    ("HS–HA shared", "K14445", "slc13a5", "Citrate transporter", "Transport"),
    ("HS–HA shared", "K02041", "phnC", "Phosphonate transporter ATP-binding", "Transport"),
    
    # HM-HA shared
    ("HM–HA shared", "K01652", "ilvB", "Acetolactate synthase large subunit", "AA biosynthesis (BCAA)"),
    ("HM–HA shared", "K01653", "ilvH", "Acetolactate synthase small subunit", "AA biosynthesis (BCAA)"),
    ("HM–HA shared", "K00053", "ilvC", "Ketol-acid reductoisomerase", "AA biosynthesis (BCAA)"),
    ("HM–HA shared", "K01687", "ilvD", "Dihydroxy-acid dehydratase", "AA biosynthesis (BCAA)"),
    ("HM–HA shared", "K01649", "leuA", "2-Isopropylmalate synthase", "AA biosynthesis (Leu)"),
    ("HM–HA shared", "K01703", "leuC", "3-Isopropylmalate dehydratase large", "AA biosynthesis (Leu)"),
    ("HM–HA shared", "K01704", "leuD", "3-Isopropylmalate dehydratase small", "AA biosynthesis (Leu)"),
    ("HM–HA shared", "K00215", "dapB", "Dihydrodipicolinate reductase", "AA biosynthesis (Lys)"),
    ("HM–HA shared", "K01586", "lysA", "Diaminopimelate decarboxylase", "AA biosynthesis (Lys)"),
    ("HM–HA shared", "K05826", "lysW", "Lysine biosynthesis carrier protein", "AA biosynthesis (Lys)"),
    ("HM–HA shared", "K01995", "livH", "BCAA transporter permease", "Transport"),
    ("HM–HA shared", "K01999", "livK", "BCAA transporter substrate-binding", "Transport"),
    ("HM–HA shared", "K00265", "gltB", "Glutamate synthase large subunit", "AA biosynthesis"),
]

distribution_colors_full = {
    "HS unique": "#FADBD8",
    "HM unique": "#D5F5E3",
    "HA unique": "#D6EAF8",
    "HS–HM shared": "#F5EEF8",
    "HS–HA shared": "#FCF3CF",
    "HM–HA shared": "#D1F2EB",
}


def create_distribution_table():
    """Create table organized by distribution category."""
    
    num_rows = len(distribution_data)
    fig_height = 1.5 + num_rows * 0.22
    
    fig, ax = plt.subplots(figsize=(13, min(fig_height, 16)))
    ax.axis('off')
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    
    headers = ['Distribution', 'KO ID', 'Gene', 'Protein/Enzyme', 'Functional Category']
    col_widths = [0.13, 0.08, 0.07, 0.43, 0.23]
    col_positions = [0.03]
    for w in col_widths[:-1]:
        col_positions.append(col_positions[-1] + w)
    
    row_height = 0.0155
    y_start = 0.95
    
    # Title
    ax.text(0.5, 0.99, 'KEGG Ortholog Distribution Across Halophilic Archaea Strains',
            ha='center', va='top', fontsize=13, fontweight='bold')
    ax.text(0.5, 0.97, 'HS: Halobacterium salinarum  |  HM: Halomicrobium kobeituziens  |  HA: Haloarcula sp.',
            ha='center', va='top', fontsize=9, style='italic', color='#555555')
    
    # Header row
    header_y = y_start
    for header, x, w in zip(headers, col_positions, col_widths):
        ax.add_patch(plt.Rectangle((x, header_y - row_height), w - 0.005, row_height,
                                   facecolor='#2C3E50', edgecolor='white', linewidth=0.5))
        ax.text(x + w/2 - 0.0025, header_y - row_height/2, header,
                ha='center', va='center', fontsize=9, fontweight='bold', color='white')
    
    current_y = header_y - row_height
    current_category = None
    
    for row in distribution_data:
        category, ko, gene, protein, function = row
        current_y -= row_height
        
        bg_color = distribution_colors_full.get(category, 'white')
        
        show_category = category if category != current_category else ""
        current_category = category
        
        row_values = [show_category, ko, gene, protein, function]
        
        for i, (val, x, w) in enumerate(zip(row_values, col_positions, col_widths)):
            ax.add_patch(plt.Rectangle((x, current_y), w - 0.005, row_height,
                                       facecolor=bg_color, edgecolor='#BDC3C7', linewidth=0.3))
            
            fontsize = 8
            fontweight = 'bold' if i == 0 and val else 'normal'
            fontstyle = 'italic' if i == 2 else 'normal'
            ha_align = 'left' if i in [0, 3, 4] else 'center'
            x_offset = 0.005 if i in [0, 3, 4] else w/2 - 0.0025
            
            ax.text(x + x_offset, current_y + row_height/2, val,
                    ha=ha_align, va='center', fontsize=fontsize, fontweight=fontweight,
                    fontstyle=fontstyle)
    
    # Legend
    legend_y = current_y - 0.04
    ax.text(0.03, legend_y + 0.015, 'Distribution Categories (total KOs):', fontsize=9, fontweight='bold')
    
    legend_items_row1 = [
        ("HS unique (71)", "#FADBD8"),
        ("HM unique (67)", "#D5F5E3"),
        ("HA unique (129)", "#D6EAF8"),
        ("HS–HM (22)", "#F5EEF8"),
    ]
    legend_items_row2 = [
        ("HS–HA (41)", "#FCF3CF"),
        ("HM–HA (133)", "#D1F2EB"),
        ("Core (755)", "#EAECEE"),
    ]
    
    legend_x = 0.03
    for label, color in legend_items_row1:
        ax.add_patch(plt.Rectangle((legend_x, legend_y - 0.018), 0.018, 0.014,
                                   facecolor=color, edgecolor='black', linewidth=0.5))
        ax.text(legend_x + 0.022, legend_y - 0.011, label, ha='left', va='center', fontsize=8)
        legend_x += 0.20
    
    legend_y -= 0.028
    legend_x = 0.03
    for label, color in legend_items_row2:
        ax.add_patch(plt.Rectangle((legend_x, legend_y - 0.018), 0.018, 0.014,
                                   facecolor=color, edgecolor='black', linewidth=0.5))
        ax.text(legend_x + 0.022, legend_y - 0.011, label, ha='left', va='center', fontsize=8)
        legend_x += 0.20
    
    plt.tight_layout()
    plt.savefig(OUTPUT_DIR / 'distributionKOtable.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.savefig(OUTPUT_DIR / 'distributionKOtable.pdf', bbox_inches='tight', facecolor='white')
    plt.close()
    print("✓ Figure B (Distribution): distributionKOtable.png + .pdf")


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("="*60)
    print("CREATING KO DISTRIBUTION TABLE FIGURES")
    print("="*60)
    print()
    
    create_functional_table()
    create_distribution_table()
    
    print()
    print("="*60)
    print("COMPLETE!")
    print("="*60)
    print(f"Output directory: {OUTPUT_DIR}")
    print()
    print("Files created:")
    print("  • functionalKOtable.png + .pdf")
    print("  • distributionKOtable.png + .pdf")

if __name__ == '__main__':
    main()