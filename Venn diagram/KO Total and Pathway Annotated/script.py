import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from pathlib import Path
import sys

SCRIPT_DIR = Path(__file__).parent.resolve()


ADJUST_TOTAL = {
    'HS only':      ( -0.002,  0.01),
    'HM only':      ( -0.02, 0.02),
    'HA only':      ( -0.1,  -0.025),
    'HS & HM':      ( 0.002, 0.02),
    'HS & HA':      ( -0.03,  -0.04),
    'HM & HA':      ( 0.08,  0.00),
    'HS & HM & HA': ( -0.03,  0.00),
}

ADJUST_PATHWAY = {
    'HS only':      ( -0.002,  0.01),
    'HM only':      ( -0.02, 0.02),
    'HA only':      ( -0.1,  -0.025),
    'HS & HM':      ( 0.002, 0.02),
    'HS & HA':      ( -0.03,  -0.04),
    'HM & HA':      ( 0.08,  0.00),
    'HS & HM & HA': ( -0.03,  0.00),
}
# ============================================================


def parse_ghostkoala(filepath):
    kos = set()
    with open(filepath, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip().replace('\r', '')
            parts = line.split('\t')
            if len(parts) >= 2 and parts[1]:
                ko = parts[1].strip()
                if ko.startswith('K') and len(ko) == 6:
                    kos.add(ko)
    return kos

def load_ko_pathway_mapping(ko_pathway_file):
    ko_to_pathways = {}
    print(f"   Reading mapping file: {ko_pathway_file.name}")
    
    with open(ko_pathway_file, 'r', encoding='utf-8') as f:
        for i, line in enumerate(f):
            line = line.strip().replace('\r', '')
            if not line: continue
            
            parts = line.split('\t')
            if len(parts) < 2:
                continue

            col_a = parts[0].strip()
            col_b = parts[1].strip()
            
            def get_ko(s):
                s = s.replace('ko:', '').replace('up:', '')
                return s if s.startswith('K') and len(s) == 6 else None

            
            ko = get_ko(col_a)
            pathway_raw = col_b
            
           
            if not ko:
                ko = get_ko(col_b)
                pathway_raw = col_a
            
            if ko:
                
                pathway = pathway_raw.replace('path:', '')
                if pathway.startswith('ko') and len(pathway) > 2 and pathway[2].isdigit():
                    pathway = 'map' + pathway[2:]
                
                
                if pathway.startswith('map'):
                    if ko not in ko_to_pathways:
                        ko_to_pathways[ko] = set()
                    ko_to_pathways[ko].add(pathway)
                
            
            if i == 0:
                print(f"   [Debug] Line 1 parsed: KO={ko}, Path={pathway}")

    return ko_to_pathways

def load_organism_pathways(pathway_file):
    pathways = set()
    with open(pathway_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip().replace('\r', '')
            parts = line.split('\t')
            if parts:
                pathways.add(parts[0].strip())
    return pathways

def filter_kos_by_organism_pathways(kos, ko_to_pathways, organism_pathways):
    filtered = set()
    for ko in kos:
        if ko in ko_to_pathways:
            if ko_to_pathways[ko] & organism_pathways:
                filtered.add(ko)
    return filtered

def create_venn(hs_kos, hm_kos, ha_kos, title, output_name, adjustments):
    only_HS = len(hs_kos - hm_kos - ha_kos)
    only_HM = len(hm_kos - hs_kos - ha_kos)
    only_HA = len(ha_kos - hs_kos - hm_kos)
    HS_HM = len((hs_kos & hm_kos) - ha_kos)
    HS_HA = len((hs_kos & ha_kos) - hm_kos)
    HM_HA = len((hm_kos & ha_kos) - hs_kos)
    all_three = len(hs_kos & hm_kos & ha_kos)
    
    total = only_HS + only_HM + only_HA + HS_HM + HS_HA + HM_HA + all_three
    
    print(f"\n{title}")
    print(f"  HS only: {only_HS}, HM only: {only_HM}, HA only: {only_HA}")
    print(f"  HS&HM: {HS_HM}, HS&HA: {HS_HA}, HM&HA: {HM_HA}, Core: {all_three}")
    
    if total == 0:
        print("  WARNING: No KOs found! Skipping this Venn diagram.")
        return
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    subsets = (only_HS, only_HM, HS_HM, only_HA, HS_HA, HM_HA, all_three)
    v = venn3(subsets=subsets, set_labels=('HS', 'HM', 'HA'), ax=ax)
    
    if v.set_labels:
        for label in v.set_labels:
            if label:
                label.set_fontsize(14)
                label.set_fontweight('bold')
    
    colors_map = {
        '100': '#ff8080',
        '010': '#80c080', 
        '001': '#8080ff',
        '110': '#c0a060',
        '101': '#c080c0',
        '011': '#80b0d0',
        '111': '#a898a8'
    }
    
    for region, color in colors_map.items():
        patch = v.get_patch_by_id(region)
        if patch:
            patch.set_color(color)
            patch.set_alpha(0.65)
            patch.set_edgecolor('none')
    
    region_to_label = {
        '100': 'HS only',
        '010': 'HM only',
        '001': 'HA only',
        '110': 'HS & HM',
        '101': 'HS & HA',
        '011': 'HM & HA',
        '111': 'HS & HM & HA'
    }
    
    counts = {
        'HS only': only_HS,
        'HM only': only_HM,
        'HA only': only_HA,
        'HS & HM': HS_HM,
        'HS & HA': HS_HA,
        'HM & HA': HM_HA,
        'HS & HM & HA': all_three
    }
    
    for region_id, label_name in region_to_label.items():
        label = v.get_label_by_id(region_id)
        if label:
            x, y = label.get_position()
            dx, dy = adjustments[label_name]
            new_x, new_y = x + dx, y + dy
            label.set_position((new_x, new_y))
            label.set_text(f'{label_name}\n{counts[label_name]}')
            label.set_fontsize(9)
    
    ax.set_title(title, fontsize=16, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / f'{output_name}.svg', bbox_inches='tight', facecolor='white')
    plt.savefig(SCRIPT_DIR / f'{output_name}.tiff', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output_name}.svg, {output_name}.tiff")


def main():
    hs_file = SCRIPT_DIR / "HS_GhostKOALA.txt"
    hm_file = SCRIPT_DIR / "HM_GhostKOALA.txt" 
    ha_file = SCRIPT_DIR / "HA_GhostKOALA.txt"
    
    # Try both filenames
    ko_pathway_file = SCRIPT_DIR / "ko_pathway.txt"
    if not ko_pathway_file.exists():
        ko_pathway_file = SCRIPT_DIR / "ko_pathways.txt"
    
    hs_pathways_file = SCRIPT_DIR / "hags.txt"
    hm_pathways_file = SCRIPT_DIR / "hmu.txt"
    ha_pathways_file = SCRIPT_DIR / "hab.txt"
    
    # Check required files
    for f in [hs_file, hm_file, ha_file]:
        if not f.exists():
            print(f"ERROR: {f.name} not found")
            sys.exit(1)
    
    hs_all = parse_ghostkoala(hs_file)
    hm_all = parse_ghostkoala(hm_file)
    ha_all = parse_ghostkoala(ha_file)
    
    print(f"Loaded GhostKOALA: HS={len(hs_all)}, HM={len(hm_all)}, HA={len(ha_all)}")
    
    create_venn(hs_all, hm_all, ha_all, 
                "Total KO Distribution", "venn_total_kos", ADJUST_TOTAL)
    
    # Pathway filtering
    pathway_files = [ko_pathway_file, hs_pathways_file, hm_pathways_file, ha_pathways_file]
    missing = [f.name for f in pathway_files if not f.exists()]
    
    if missing:
        print(f"\nMissing files for pathway Venn: {missing}")
    else:
        ko_to_pathways = load_ko_pathway_mapping(ko_pathway_file)
        print(f"\nLoaded KO-pathway mapping: {len(ko_to_pathways)} KOs")
        
        hs_pathways = load_organism_pathways(hs_pathways_file)
        hm_pathways = load_organism_pathways(hm_pathways_file)
        ha_pathways = load_organism_pathways(ha_pathways_file)
        print(f"Organism pathways: HS={len(hs_pathways)}, HM={len(hm_pathways)}, HA={len(ha_pathways)}")
        
        # Debug: show sample pathways
        sample_ko = list(ko_to_pathways.keys())[0] if ko_to_pathways else None
        sample_hs = list(hs_pathways)[:3] if hs_pathways else []
        if sample_ko:
            print(f"  Sample KO mapping: {sample_ko} -> {list(ko_to_pathways[sample_ko])[:3]}")
        print(f"  Sample HS pathways: {sample_hs}")
        
        hs_filtered = filter_kos_by_organism_pathways(hs_all, ko_to_pathways, hs_pathways)
        hm_filtered = filter_kos_by_organism_pathways(hm_all, ko_to_pathways, hm_pathways)
        ha_filtered = filter_kos_by_organism_pathways(ha_all, ko_to_pathways, ha_pathways)
        
        print(f"After filtering: HS={len(hs_filtered)}, HM={len(hm_filtered)}, HA={len(ha_filtered)}")
        
        create_venn(hs_filtered, hm_filtered, ha_filtered, 
                   "Pathway-Mapped KO Distribution", "venn_pathway_kos", ADJUST_PATHWAY)

if __name__ == '__main__':
    main()