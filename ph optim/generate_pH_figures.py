"""
pH Optimization Figure Generator (3 repeats)
=============================================
Place this script in the same folder as:
  - ph_repeat_1.xlsx
  - ph_repeat_2.xlsx
  - ph_repeat_3.xlsx

Run: python generate_pH_figures.py
"""

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import os

BUFFERS = ['PIPES', 'HEPES', 'TRIS']
STRAINS = ['HS', 'HM', 'HA']
TIMES = [0, 24, 48, 72]
STRAIN_NAMES = {'HS':'H. salinarum KBTZ01','HM':'Halomicrobium sp. KBTZ05','HA':'Haloarcula sp. KBTZ06'}
STRAIN_COLORS = {'HS':'#4393c3','HM':'#d6604d','HA':'#5ab4ac'}
BUF_COLORS = {'PIPES':'#2166AC','HEPES':'#B2182B','TRIS':'#4DAF4A'}
TIME_COLORS = {0:'#fee5d9', 24:'#fcae91', 48:'#fb6a4a', 72:'#a50f15'}

PH_72H = {
    'HS': {'PIPES':[7.602,7.7,7.7,7.830],'HEPES':[7.875,7.97,7.955,7.95],'TRIS':[7.93,8.09,8.003,8.2]},
    'HM': {'PIPES':[4.790,5.229,6.370,7.001],'HEPES':[6.410,6.875,7.445,7.696],'TRIS':[7.387,7.6,7.8,8.025]},
    'HA': {'PIPES':[6.035,7.42,7.4,7.4],'HEPES':[7.65,7.8,7.9,7.9],'TRIS':[7.905,8.07,8.105,8.167]},
}
INITIAL_PH = {'PIPES':[5.0,5.5,6.0,6.5],'HEPES':[6.0,6.5,7.0,7.5],'TRIS':[6.5,7.0,7.5,8.0]}

plt.rcParams.update({'font.family':'Arial','font.size':7,'axes.linewidth':0.7,
    'xtick.major.width':0.7,'ytick.major.width':0.7,'figure.dpi':300,'savefig.dpi':300})

def parse_repeat(filepath):
    df = pd.read_excel(filepath, header=None)
    data = {}
    strain_rows = {'Control':{'PIPES':3,'HEPES':5,'TRIS':7},'HS':{'PIPES':11,'HEPES':13,'TRIS':15},
                   'HM':{'PIPES':19,'HEPES':21,'TRIS':23},'HA':{'PIPES':27,'HEPES':29,'TRIS':31}}
    time_col_start = {0:1, 24:7, 48:13, 72:19}
    for strain in ['Control'] + STRAINS:
        data[strain] = {}
        for buf in BUFFERS:
            row = strain_rows[strain][buf]
            phs = []
            for i in range(4):
                raw = df.iloc[row-1, time_col_start[0]+i]
                phs.append(float(raw.replace('pH ','')) if isinstance(raw,str) and raw.startswith('pH') else float(raw))
            data[strain][buf] = {'pH': phs}
            for t in TIMES:
                data[strain][buf][t] = [float(df.iloc[row, time_col_start[t]+i]) for i in range(4)]
    return data

script_dir = os.path.dirname(os.path.abspath(__file__))
repeats = []
print("Loading data...")
for f in ['ph_repeat_1.xlsx','ph_repeat_2.xlsx','ph_repeat_3.xlsx']:
    p = os.path.join(script_dir, f)
    repeats.append(parse_repeat(p if os.path.exists(p) else f))
    print(f"  Loaded {f}")

def corrected_mean(strain, buf, pi, t):
    return max(np.mean([r[strain][buf][t][pi] - r['Control'][buf][t][pi] for r in repeats]), 0)

def corrected_sem(strain, buf, pi, t):
    vals = [r[strain][buf][t][pi] - r['Control'][buf][t][pi] for r in repeats]
    return np.std(vals, ddof=1) / np.sqrt(len(vals))

# FIGURE 3: 3x3 grid
print("Generating Figure 3...")
fig, axes = plt.subplots(3, 3, figsize=(7.2, 7.0), sharey=True)
panel = 0
for si, strain in enumerate(STRAINS):
    for bi, buf in enumerate(BUFFERS):
        ax = axes[si][bi]
        phs = repeats[0][strain][buf]['pH']
        x = np.arange(4); w = 0.2
        for ti, t in enumerate(TIMES):
            vals = [corrected_mean(strain, buf, pi, t) for pi in range(4)]
            sems = [corrected_sem(strain, buf, pi, t) for pi in range(4)]
            ax.bar(x+(ti-1.5)*w, vals, w, yerr=sems, color=TIME_COLORS[t],
                   edgecolor='black', linewidth=0.3, alpha=0.9,
                   capsize=1.5, error_kw={'linewidth':0.4,'capthick':0.4},
                   label=f'{t}h' if (si==0 and bi==0) else None)
        ax.set_xticks(x); ax.set_xticklabels([f'{p}' for p in phs], fontsize=7)
        ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
        ax.tick_params(labelsize=6.5); ax.set_ylim(0, 1.6)
        ax.text(-0.05, 1.03, 'ABCDEFGHI'[panel], transform=ax.transAxes, fontsize=9, fontweight='bold')
        panel += 1
        if si == 0: ax.set_title(buf, fontsize=9, fontweight='bold', color=BUF_COLORS[buf])
        if si == 2: ax.set_xlabel('pH', fontsize=7.5)
        if bi == 0: ax.set_ylabel('OD$_{660}$ (corrected)', fontsize=7)
    axes[si][2].text(1.15, 0.5, STRAIN_NAMES[strain], transform=axes[si][2].transAxes,
                     fontsize=7.5, fontstyle='italic', fontweight='bold', rotation=-90, va='center', ha='center')
axes[0][0].legend(fontsize=6, frameon=False, loc='upper left', title='Time', title_fontsize=6.5)
plt.tight_layout(rect=[0, 0, 0.95, 1])
fig.savefig(os.path.join(script_dir, 'fig3_all_strains_combined.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved fig3_all_strains_combined.png")

# FIGURE 4: Buffer comparison at 48h
print("Generating Figure 4...")
fig, ax = plt.subplots(figsize=(4.0, 3.0))
x = np.arange(3); w = 0.25
for si, strain in enumerate(STRAINS):
    means, sems = [], []
    for buf in BUFFERS:
        ph_vals = [corrected_mean(strain, buf, pi, 48) for pi in range(4)]
        means.append(np.mean(ph_vals))
        sems.append(np.std(ph_vals, ddof=1) / np.sqrt(4))
    ax.bar(x+(si-1)*w, means, w, yerr=sems, color=STRAIN_COLORS[strain],
           edgecolor='black', linewidth=0.4, capsize=3, error_kw={'linewidth':0.5,'capthick':0.5},
           alpha=0.9, label=STRAIN_NAMES[strain])
ax.set_xticks(x); ax.set_xticklabels(BUFFERS, fontsize=9)
ax.set_ylabel('OD$_{660}$ at 48h\n(mean across pH, corrected)', fontsize=7.5)
ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
ax.tick_params(labelsize=7); ax.set_ylim(0, 1.4)
ax.legend(fontsize=5.5, frameon=False, loc='upper right')
ax.text(-0.08, 1.05, 'A', transform=ax.transAxes, fontsize=10, fontweight='bold')
plt.tight_layout()
fig.savefig(os.path.join(script_dir, 'fig4_buffer_comparison.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved fig4_buffer_comparison.png")

# FIGURE 5: pH shift
print("Generating Figure 5...")
BUF_MARKERS = {'PIPES':'o','HEPES':'s','TRIS':'^'}
fig, axes = plt.subplots(1, 3, figsize=(6.0, 2.5))
for si, strain in enumerate(STRAINS):
    ax = axes[si]
    for buf in BUFFERS:
        ax.scatter(INITIAL_PH[buf], PH_72H[strain][buf], s=35, color=STRAIN_COLORS[strain],
                  marker=BUF_MARKERS[buf], edgecolors='black', linewidths=0.4, zorder=3)
    ax.plot([4,9],[4,9], 'k--', linewidth=0.5, alpha=0.3)
    ax.set_xlim(4.5, 8.5); ax.set_ylim(4.5, 8.5); ax.set_aspect('equal')
    ax.spines['top'].set_visible(False); ax.spines['right'].set_visible(False)
    ax.tick_params(labelsize=6); ax.set_xlabel('Initial pH', fontsize=7)
    ax.set_title(STRAIN_NAMES[strain], fontsize=7, fontweight='bold', fontstyle='italic')
    ax.text(-0.1, 1.05, chr(65+si), transform=ax.transAxes, fontsize=10, fontweight='bold')
    if si == 0: ax.set_ylabel('Final pH (72 h)', fontsize=7)
legend_elements = [
    Line2D([0],[0],marker='o',color='w',markerfacecolor='gray',markeredgecolor='black',markersize=5,label='PIPES'),
    Line2D([0],[0],marker='s',color='w',markerfacecolor='gray',markeredgecolor='black',markersize=5,label='HEPES'),
    Line2D([0],[0],marker='^',color='w',markerfacecolor='gray',markeredgecolor='black',markersize=5,label='TRIS')]
axes[2].legend(handles=legend_elements, fontsize=5.5, frameon=False, loc='lower right')
plt.tight_layout()
fig.savefig(os.path.join(script_dir, 'fig5_pH_shift.png'), dpi=300, bbox_inches='tight')
plt.close()
print("  Saved fig5_pH_shift.png")

print(f"\nDone! All figures saved (n={len(repeats)} biological replicates).")
