# install matplotlib-venn (for example with pip)

import matplotlib.pyplot as plt
from matplotlib_venn import venn3, venn3_circles
from pathlib import Path
import os

script_dir = Path(__file__).parent
project_root = script_dir.parent if script_dir.name == 'scripts' else script_dir
os.chdir(project_root)

GHOSTKOALA_DIR = Path('annotation_outputs/GhostKOALA')
OUTPUT_DIR = Path('results/figures')
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

def read_ghostkoala(filepath):
    kos = set()
    with open(filepath, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2 and parts[1].startswith('K'):
                kos.add(parts[1])
    return kos

hs_kos = read_ghostkoala(GHOSTKOALA_DIR / 'HS_GhostKOALA.txt')
hm_kos = read_ghostkoala(GHOSTKOALA_DIR / 'HM_GhostKOALA.txt')
ha_kos = read_ghostkoala(GHOSTKOALA_DIR / 'HA_GhostKOALA.txt')

hs_only = hs_kos - hm_kos - ha_kos
hm_only = hm_kos - hs_kos - ha_kos
ha_only = ha_kos - hs_kos - hm_kos
hs_hm = (hs_kos & hm_kos) - ha_kos
hs_ha = (hs_kos & ha_kos) - hm_kos
hm_ha = (hm_kos & ha_kos) - hs_kos
core = hs_kos & hm_kos & ha_kos

venn_values = (len(hs_only), len(hm_only), len(hs_hm), 
               len(ha_only), len(hs_ha), len(hm_ha), len(core))

fig, ax = plt.subplots(figsize=(10, 10))

colors = ('#FF6B6B', '#4ECDC4', '#95E1D3')

v = venn3(subsets=venn_values,
          set_labels=('', '', ''),
          ax=ax,
          set_colors=colors,
          alpha=0.6)

venn3_circles(subsets=venn_values, linestyle='-', linewidth=2, color='white', ax=ax)

for idx in ['100', '010', '001', '110', '101', '011', '111']:
    if v.get_label_by_id(idx):
        v.get_label_by_id(idx).set_text('')

ax.axis('off')

ax.axis('off')
plt.tight_layout()
plt.savefig(OUTPUT_DIR / 'venn_ko_overlap.png', dpi=300, bbox_inches='tight', facecolor='white')
plt.savefig(OUTPUT_DIR / 'venn_ko_overlap.pdf', bbox_inches='tight', facecolor='white')
plt.close()

print("=" * 60)
print("VENN DIAGRAM - KO OVERLAP")
print("=" * 60)
print()
print("Colors: HS=#FF6B6B, HM=#4ECDC4, HA=#95E1D3")
print()
print("Counts for labels:")
print()
print(f"HS only:      {len(hs_only)}")
print(f"HM only:      {len(hm_only)}")
print(f"HA only:      {len(ha_only)}")
print(f"HS-HM:        {len(hs_hm)}")
print(f"HS-HA:        {len(hs_ha)}")
print(f"HM-HA:        {len(hm_ha)}")
print(f"Core:         {len(core)}")
print()
print(f"Total: HS={len(hs_kos)}, HM={len(hm_kos)}, HA={len(ha_kos)}")
print()
print(f"✓ Saved to {OUTPUT_DIR}")
print("=" * 60)