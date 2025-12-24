#!/usr/bin/env python3

import matplotlib.pyplot as plt

output_dir = 'C:/Users/reison/Desktop/12.12/halophilic-archaea-metabolism/results/'

ko_data = [
    ("HS unique", "K00926", "arcC", "Carbamate kinase"),
    ("HS unique", "K08177", "tnaB", "Tryptophan permease"),
    ("HS unique", "K00005", "gldA", "Glycerol dehydrogenase"),
    ("HS unique", "K01546", "ntpA", "V-type ATPase subunit A"),
    ("HS unique", "K01547", "ntpB", "V-type ATPase subunit B"),
    ("HS unique", "K03469", "bat", "Bacterio-opsin activator"),
    
    ("HM unique", "K01601", "rbcL", "RuBisCO large subunit"),
    ("HM unique", "K00931", "proB", "Glutamate 5-kinase"),
    ("HM unique", "K00147", "proA", "Glutamate-5-semialdehyde dehydrogenase"),
    ("HM unique", "K00286", "proC", "Pyrroline-5-carboxylate reductase"),
    ("HM unique", "K02510", "hisG", "ATP phosphoribosyltransferase"),
    ("HM unique", "K01191", "treA", "Trehalase"),
    ("HM unique", "K01192", "treF", "Trehalase"),
    ("HM unique", "K01495", "folE", "GTP cyclohydrolase I"),
    
    
    ("HA unique", "K01428", "ureA", "Urease gamma subunit"),
    ("HA unique", "K01429", "ureB", "Urease beta subunit"),
    ("HA unique", "K01430", "ureC", "Urease alpha subunit"),
    ("HA unique", "K11959", "urtA", "Urea ABC transporter substrate-binding protein"),
    ("HA unique", "K11960", "urtB", "Urea ABC transporter permease"),
    ("HA unique", "K02784", "crr", "PTS glucose-specific EIIA component"),
    ("HA unique", "K02770", "ptsG", "PTS glucose-specific EIICB component"),
    ("HA unique", "K00376", "nosZ", "Nitrous oxide reductase"),
    ("HA unique", "K17316", "norB", "Nitric oxide reductase subunit B"),
    ("HA unique", "K00697", "glgA", "Glycogen synthase"),
    ("HA unique", "K01178", "amyA", "Alpha-amylase"),
    
    ("HS–HM shared", "K03892", "bop", "Bacteriorhodopsin"),
    ("HS–HM shared", "K02052", "potA", "Polyamine ABC transporter ATP-binding protein"),
    ("HS–HM shared", "K02053", "potB", "Polyamine ABC transporter permease"),
    ("HS–HM shared", "K01620", "ltaE", "Threonine aldolase"),
    ("HS–HM shared", "K01183", "chiA", "Chitinase"),
    ("HS–HM shared", "K06901", "phrB", "DNA photolyase"),
    
    ("HS–HA shared", "K01745", "hutH", "Histidine ammonia-lyase"),
    ("HS–HA shared", "K01712", "hutU", "Urocanate hydratase"),
    ("HS–HA shared", "K23316", "hutI", "Imidazolonepropionase"),
    ("HS–HA shared", "K23684", "hutG", "Formiminoglutamate hydrolase"),
    ("HS–HA shared", "K01667", "tnaA", "Tryptophanase"),
    ("HS–HA shared", "K00162", "pdhB", "Pyruvate dehydrogenase E1 beta subunit"),
    ("HS–HA shared", "K14445", "slc13a5", "Citrate transporter"),
    ("HS–HA shared", "K02041", "phnC", "Phosphonate ABC transporter ATP-binding protein"),
    
    ("HM–HA shared", "K01652", "ilvB", "Acetolactate synthase large subunit"),
    ("HM–HA shared", "K01653", "ilvH", "Acetolactate synthase small subunit"),
    ("HM–HA shared", "K00053", "ilvC", "Ketol-acid reductoisomerase"),
    ("HM–HA shared", "K01687", "ilvD", "Dihydroxy-acid dehydratase"),
    ("HM–HA shared", "K01649", "leuA", "2-Isopropylmalate synthase"),
    ("HM–HA shared", "K01703", "leuC", "3-Isopropylmalate dehydratase large subunit"),
    ("HM–HA shared", "K01704", "leuD", "3-Isopropylmalate dehydratase small subunit"),
    ("HM–HA shared", "K00215", "dapB", "Dihydrodipicolinate reductase"),
    ("HM–HA shared", "K01586", "lysA", "Diaminopimelate decarboxylase"),
    ("HM–HA shared", "K05826", "lysW", "Lysine biosynthesis carrier protein"),
    ("HM–HA shared", "K01995", "livH", "BCAA ABC transporter permease"),
    ("HM–HA shared", "K01999", "livK", "BCAA ABC transporter substrate-binding protein"),
    ("HM–HA shared", "K00265", "gltB", "Glutamate synthase large subunit"),
]

category_colors = {
    "HS unique": "#FADBD8",
    "HM unique": "#D5F5E3",
    "HA unique": "#D6EAF8",
    "HS–HM shared": "#F5EEF8",
    "HS–HA shared": "#FCF3CF",
    "HM–HA shared": "#D1F2EB",
}

category_counts = {
    "HS unique": 71,
    "HM unique": 67,
    "HA unique": 129,
    "HS–HM shared": 22,
    "HS–HA shared": 41,
    "HM–HA shared": 133,
}

def create_table_figure():
    fig, ax = plt.subplots(figsize=(10, 17))
    ax.axis('off')
    
    headers = ['Distribution', 'KO ID', 'Gene', 'Protein/Enzyme']
    col_widths = [0.18, 0.12, 0.12, 0.52]
    
    y_start = 0.96
    row_height = 0.015
    x_positions = [0.03]
    for w in col_widths[:-1]:
        x_positions.append(x_positions[-1] + w)
    
    header_y = y_start
    for i, (header, x, w) in enumerate(zip(headers, x_positions, col_widths)):
        ax.add_patch(plt.Rectangle((x, header_y - row_height), w - 0.005, row_height,
                                   facecolor='#2C3E50', edgecolor='white', linewidth=0.5))
        ax.text(x + w/2 - 0.0025, header_y - row_height/2, header,
                ha='center', va='center', fontsize=9, fontweight='bold', color='white')
    
    current_y = header_y - row_height
    current_category = None
    
    for row in ko_data:
        category, ko, gene, function = row
        current_y -= row_height
        
        bg_color = category_colors.get(category, 'white')
        
        show_category = category if category != current_category else ""
        current_category = category
        row_data = [show_category, ko, gene, function]
        
        for i, (val, x, w) in enumerate(zip(row_data, x_positions, col_widths)):
            ax.add_patch(plt.Rectangle((x, current_y), w - 0.005, row_height,
                                       facecolor=bg_color, edgecolor='#BDC3C7', linewidth=0.3))
            
            fontsize = 8
            fontweight = 'bold' if i == 0 and val else 'normal'
            fontstyle = 'italic' if i == 2 else 'normal'
            ha = 'left' if i in [0, 3] else 'center'
            x_offset = 0.008 if i in [0, 3] else w/2 - 0.0025
            
            ax.text(x + x_offset, current_y + row_height/2, val,
                    ha=ha, va='center', fontsize=fontsize, fontweight=fontweight,
                    fontstyle=fontstyle, color='black')
    
    ax.text(0.5, 0.99, 'Key KEGG Ortholog Distribution Across Halophilic Archaea Strains',
            ha='center', va='top', fontsize=12, fontweight='bold', transform=ax.transAxes)
    ax.text(0.5, 0.975, 'HS: Halobacterium salinarum  |  HM: Halomicrobium kobeituziens  |  HA: Haloarcula sp.',
            ha='center', va='top', fontsize=9, style='italic', transform=ax.transAxes)
    
    legend_y = current_y - 0.022
    legend_items = [
        (f"HS unique ({category_counts['HS unique']})", "#FADBD8"),
        (f"HM unique ({category_counts['HM unique']})", "#D5F5E3"),
        (f"HA unique ({category_counts['HA unique']})", "#D6EAF8"),
        (f"HS–HM shared ({category_counts['HS–HM shared']})", "#F5EEF8"),
        (f"HS–HA shared ({category_counts['HS–HA shared']})", "#FCF3CF"),
        (f"HM–HA shared ({category_counts['HM–HA shared']})", "#D1F2EB"),
    ]
    
    legend_x = 0.03
    for label, color in legend_items:
        ax.add_patch(plt.Rectangle((legend_x, legend_y), 0.018, 0.012,
                                   facecolor=color, edgecolor='black', linewidth=0.5))
        ax.text(legend_x + 0.023, legend_y + 0.006, label, ha='left', va='center', fontsize=7.5)
        legend_x += 0.155
    
    plt.tight_layout()
    plt.savefig(output_dir + 'figure_table_KO_distribution.png', dpi=300, bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.savefig(output_dir + 'figure_table_KO_distribution.pdf', bbox_inches='tight',
                facecolor='white', edgecolor='none')
    plt.close()

if __name__ == '__main__':
    create_table_figure()