"""
UNANNOTATED PROTEIN BLAST ANALYSIS
====================================
Dependencies:
  pip install matplotlib matplotlib-venn
  BLAST+ (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

This script performs reciprocal best hit (RBH) analysis on unannotated proteins
from three strains (HS, HM, HA) and generates Venn diagrams and CSV files.
"""

import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from pathlib import Path
from collections import defaultdict
import subprocess
import sys
import csv

SCRIPT_DIR = Path(__file__).parent.resolve()

# BLAST parameters
EVALUE_CUTOFF = 1e-5
MIN_IDENTITY = 60             
MIN_QUERY_COV = 70            
NUM_THREADS = 8             

# Venn diagram label position adjustments
ADJUSTMENTS = {
    'HS only':      (-0.002,  0.01),
    'HM only':      (-0.02,   0.02),
    'HA only':      ( 0.0,   -0.025),
    'HS & HM':      ( 0.002,  0.02),
    'HS & HA':      (-0.03,  -0.04),
    'HM & HA':      ( 0.08,   0.00),
    'HS & HM & HA': (-0.03,   0.00),
}
# ============================================================


def get_unannotated_ids(ghostkoala_file):
    """
    Parse GhostKOALA output to identify annotated vs unannotated proteins.
    
    Args:
        ghostkoala_file: Path to GhostKOALA annotation file
        
    Returns:
        Tuple of (unannotated_set, annotated_set)
    """
    unannotated = set()
    annotated = set()
    with open(ghostkoala_file, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip().replace('\r', '')
            if not line: continue
            parts = line.split('\t')
            protein_id = parts[0].strip()
            ko = parts[1].strip() if len(parts) >= 2 else ''
            if ko and ko.startswith('K') and len(ko) == 6:
                annotated.add(protein_id)
            else:
                unannotated.add(protein_id)
    return unannotated, annotated


def extract_sequences(fasta_file, target_ids, output_file, tag):
    """
    Extract sequences for target protein IDs from FASTA file.
    
    Args:
        fasta_file: Input FASTA file
        target_ids: Set of protein IDs to extract
        output_file: Output FASTA file
        tag: Organism tag to prepend to headers
        
    Returns:
        Number of sequences extracted
    """
    count = 0
    with open(fasta_file, 'r') as fin, open(output_file, 'w') as fout:
        writing = False
        for line in fin:
            line = line.replace('\r', '')
            if line.startswith('>'):
                header_id = line[1:].split()[0]
                writing = header_id in target_ids
                if writing:
                    fout.write(f'>{tag}|{header_id}\n')
                    count += 1
            elif writing:
                fout.write(line)
    return count


def run_blast(combined_fasta, output_file, evalue, threads):
    """
    Run BLASTp all-vs-all comparison.
    
    Args:
        combined_fasta: Combined FASTA file with all proteins
        output_file: Output file for BLAST results
        evalue: E-value cutoff
        threads: Number of threads to use
    """
    db_name = str(SCRIPT_DIR / 'unannotated_db')

    print("  Creating BLAST database...")
    try:
        subprocess.run([
            'makeblastdb',
            '-in', str(combined_fasta),
            '-dbtype', 'prot',
            '-out', db_name
        ], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: makeblastdb failed: {e.stderr}")
        sys.exit(1)

    print(f"  Running BLASTp (e-value < {evalue}, {threads} threads)...")
    print("  This may take several minutes...")
    try:
        subprocess.run([
            'blastp',
            '-query', str(combined_fasta),
            '-db', db_name,
            '-out', str(output_file),
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs',
            '-evalue', str(evalue),
            '-max_target_seqs', '10',
            '-num_threads', str(threads)
        ], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: blastp failed: {e.stderr}")
        sys.exit(1)

    hits = sum(1 for _ in open(output_file))
    print(f"  BLAST complete: {hits} hits found")


def find_rbh(blast_file, min_identity, min_qcov):
    """
    Identify reciprocal best hits from BLAST results.
    
    Args:
        blast_file: BLAST output file
        min_identity: Minimum percent identity threshold
        min_qcov: Minimum query coverage threshold
        
    Returns:
        Set of RBH pairs (tuples of protein IDs)
    """
    best_hits = {}
    with open(blast_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            qid, sid = parts[0], parts[1]
            pident = float(parts[2])
            evalue = float(parts[10])
            qcovs = float(parts[12])

            # Skip self-hits and same-organism hits
            if qid == sid: continue
            q_org = qid.split('|')[0]
            s_org = sid.split('|')[0]
            if q_org == s_org: continue
            
            # Apply filters
            if pident < min_identity: continue
            if qcovs < min_qcov: continue

            # Keep best hit per query-target organism pair
            key = (qid, s_org)
            if key not in best_hits or evalue < best_hits[key][1]:
                best_hits[key] = (sid, evalue, pident)

    # Find reciprocal best hits
    rbh_pairs = set()
    for (qid, t_org), (sid, _, _) in best_hits.items():
        q_org = qid.split('|')[0]
        reverse_key = (sid, q_org)
        if reverse_key in best_hits and best_hits[reverse_key][0] == qid:
            pair = tuple(sorted([qid, sid]))
            rbh_pairs.add(pair)

    return rbh_pairs


def create_venn(counts, title, output_name):
    """
    Create and save Venn diagram visualization.
    
    Args:
        counts: Dictionary of counts for each Venn region
        title: Plot title
        output_name: Output filename prefix
    """
    fig, ax = plt.subplots(figsize=(10, 8))

    subsets = (
        counts['hs_only'], counts['hm_only'], counts['hs_hm'],
        counts['ha_only'], counts['hs_ha'], counts['hm_ha'],
        counts['all_three']
    )
    v = venn3(subsets=subsets, set_labels=('HS', 'HM', 'HA'), ax=ax)

    if v.set_labels:
        for label in v.set_labels:
            if label:
                label.set_fontsize(14)
                label.set_fontweight('bold')

    colors_map = {
        '100': '#ff8080', '010': '#80c080', '001': '#8080ff',
        '110': '#c0a060', '101': '#c080c0', '011': '#80b0d0',
        '111': '#a898a8'
    }

    for region, color in colors_map.items():
        patch = v.get_patch_by_id(region)
        if patch:
            patch.set_color(color)
            patch.set_alpha(0.65)
            patch.set_edgecolor('none')

    region_labels = {
        '100': ('HS only',      counts['hs_only']),
        '010': ('HM only',      counts['hm_only']),
        '001': ('HA only',      counts['ha_only']),
        '110': ('HS & HM',      counts['hs_hm']),
        '101': ('HS & HA',      counts['hs_ha']),
        '011': ('HM & HA',      counts['hm_ha']),
        '111': ('HS & HM & HA', counts['all_three']),
    }

    for region_id, (label_name, count) in region_labels.items():
        label = v.get_label_by_id(region_id)
        if label:
            x, y = label.get_position()
            dx, dy = ADJUSTMENTS[label_name]
            label.set_position((x + dx, y + dy))
            label.set_text(f'{label_name}\n{count}')
            label.set_fontsize(9)

    ax.set_title(title, fontsize=16, fontweight='bold')

    plt.tight_layout()
    plt.savefig(SCRIPT_DIR / f'{output_name}.svg', bbox_inches='tight', facecolor='white')
    plt.savefig(SCRIPT_DIR / f'{output_name}.tiff', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output_name}.svg, {output_name}.tiff")


def save_protein_lists_to_csv(protein_sets, rbh_pairs, output_prefix):
    """
    Save protein lists for each category to separate CSV files.
    
    Args:
        protein_sets: Dictionary containing sets of proteins for each category
        rbh_pairs: Set of RBH pairs for finding partner proteins
        output_prefix: Prefix for output CSV files
    """
    print("\nSTEP 6: Saving protein lists to CSV files...")
    
    # Build a mapping of protein to its RBH partners (actual protein IDs)
    protein_partners = defaultdict(set)
    for p1, p2 in rbh_pairs:
        protein_partners[p1].add(p2)
        protein_partners[p2].add(p1)
    
    categories = {
        'HS_only': protein_sets['hs_only'],
        'HM_only': protein_sets['hm_only'],
        'HA_only': protein_sets['ha_only'],
        'HS_and_HM_only': protein_sets['hs_hm'],
        'HS_and_HA_only': protein_sets['hs_ha'],
        'HM_and_HA_only': protein_sets['hm_ha'],
        'Core_all_three': protein_sets['all_three']
    }
    
    for category, proteins in categories.items():
        if not proteins:
            continue
            
        csv_file = SCRIPT_DIR / f"{output_prefix}_{category}.csv"
        
        with open(csv_file, 'w', newline='', encoding='utf-8') as f:
            writer = csv.writer(f)
            
            # Write header
            writer.writerow(['Protein_ID', 'Organism', 'Original_ID', 'RBH_Partner_IDs'])
            
            # Write protein data
            for protein_id in sorted(proteins):
                org = protein_id.split('|')[0]
                original_id = '|'.join(protein_id.split('|')[1:]) if '|' in protein_id else protein_id
                
                # Get actual RBH partner protein IDs
                partners = protein_partners.get(protein_id, set())
                partner_list = '; '.join(sorted(partners)) if partners else 'None'
                
                writer.writerow([protein_id, org, original_id, partner_list])
        
        print(f"  Saved {len(proteins)} proteins to: {csv_file.name}")
    
    # Create summary CSV
    summary_file = SCRIPT_DIR / f"{output_prefix}_summary.csv"
    with open(summary_file, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(['Category', 'Count', 'Description'])
        writer.writerow(['HS_only', len(protein_sets['hs_only']), 'Unique to HS strain'])
        writer.writerow(['HM_only', len(protein_sets['hm_only']), 'Unique to HM strain'])
        writer.writerow(['HA_only', len(protein_sets['ha_only']), 'Unique to HA strain'])
        writer.writerow(['HS_and_HM_only', len(protein_sets['hs_hm']), 'Shared between HS and HM only'])
        writer.writerow(['HS_and_HA_only', len(protein_sets['hs_ha']), 'Shared between HS and HA only'])
        writer.writerow(['HM_and_HA_only', len(protein_sets['hm_ha']), 'Shared between HM and HA only'])
        writer.writerow(['Core_all_three', len(protein_sets['all_three']), 'Shared among all three strains'])
        writer.writerow(['Total', sum(len(v) for v in protein_sets.values()), 'Total unannotated proteins analyzed'])
    
    print(f"  Saved summary to: {summary_file.name}")


def main():
    print("=" * 60)
    print("  UNANNOTATED PROTEIN BLAST ANALYSIS")
    print("=" * 60)

    organisms = {
        'HS': {
            'ghostkoala': SCRIPT_DIR / 'HS_GhostKOALA.txt',
            'fasta':      SCRIPT_DIR / 'HS_protein.faa',
        },
        'HM': {
            'ghostkoala': SCRIPT_DIR / 'HM_GhostKOALA.txt',
            'fasta':      SCRIPT_DIR / 'HM_protein.faa',
        },
        'HA': {
            'ghostkoala': SCRIPT_DIR / 'HA_GhostKOALA.txt',
            'fasta':      SCRIPT_DIR / 'HA_protein.faa',
        },
    }

    # Check input files exist
    for org, files in organisms.items():
        for name, path in files.items():
            if not path.exists():
                print(f"ERROR: {path.name} not found!")
                sys.exit(1)

    # STEP 1: Identify unannotated proteins
    print("\nSTEP 1: Identifying unannotated proteins...")
    unannotated_ids = {}
    for org, files in organisms.items():
        unannot, annot = get_unannotated_ids(files['ghostkoala'])
        unannotated_ids[org] = unannot
        total = len(annot) + len(unannot)
        print(f"  {org}: {len(annot)} annotated, {len(unannot)} unannotated ({len(unannot)*100//total}%)")

    # STEP 2: Extract unannotated sequences
    print("\nSTEP 2: Extracting unannotated sequences...")
    fasta_files = {}
    for org, files in organisms.items():
        out_file = SCRIPT_DIR / f'{org}_unannotated.faa'
        count = extract_sequences(files['fasta'], unannotated_ids[org], out_file, org)
        print(f"  {org}: {count} sequences → {out_file.name}")
        fasta_files[org] = out_file

    # Combine all sequences
    combined = SCRIPT_DIR / 'all_unannotated.faa'
    with open(combined, 'w') as fout:
        for org in ['HS', 'HM', 'HA']:
            with open(fasta_files[org]) as fin:
                fout.write(fin.read())

    # Count proteins per organism
    all_proteins = defaultdict(list)
    with open(combined) as f:
        for line in f:
            if line.startswith('>'):
                prot_id = line[1:].strip().split()[0]
                org = prot_id.split('|')[0]
                all_proteins[org].append(prot_id)

    total = sum(len(v) for v in all_proteins.values())
    print(f"  Combined: {total} proteins")

    # STEP 3: Run BLAST
    print("\nSTEP 3: Running BLASTp all-vs-all...")
    blast_output = SCRIPT_DIR / 'blast_results.txt'
    run_blast(combined, blast_output, EVALUE_CUTOFF, NUM_THREADS)

    # STEP 4: Find RBH
    print(f"\nSTEP 4: Finding Reciprocal Best Hits (>={MIN_IDENTITY}% identity, >={MIN_QUERY_COV}% coverage)...")
    rbh_pairs = find_rbh(blast_output, MIN_IDENTITY, MIN_QUERY_COV)
    print(f"  RBH pairs found: {len(rbh_pairs)}")

    # STEP 5: Classify proteins
    print("\nSTEP 5: Classifying proteins and creating Venn diagram...")

    # Build connections map (organism level, not protein IDs)
    connections = defaultdict(set)
    for p1, p2 in rbh_pairs:
        org1 = p1.split('|')[0]
        org2 = p2.split('|')[0]
        connections[p1].add(org2)
        connections[p2].add(org1)

    # Get all protein sets
    hs_set = set(all_proteins['HS'])
    hm_set = set(all_proteins['HM'])
    ha_set = set(all_proteins['HA'])

    # Identify proteins with connections to each organism
    hs_has_hm = {p for p in hs_set if 'HM' in connections.get(p, set())}
    hs_has_ha = {p for p in hs_set if 'HA' in connections.get(p, set())}
    hm_has_hs = {p for p in hm_set if 'HS' in connections.get(p, set())}
    hm_has_ha = {p for p in hm_set if 'HA' in connections.get(p, set())}
    ha_has_hs = {p for p in ha_set if 'HS' in connections.get(p, set())}
    ha_has_hm = {p for p in ha_set if 'HM' in connections.get(p, set())}

    # Create sets for unique proteins
    hs_only_set = hs_set - hs_has_hm - hs_has_ha
    hm_only_set = hm_set - hm_has_hs - hm_has_ha
    ha_only_set = ha_set - ha_has_hs - ha_has_hm

    # Create sets for shared proteins (include all proteins from both/all organisms)
    # For CSV export, we want ALL proteins in each category
    hs_hm_only_set = (hs_has_hm - hs_has_ha) | (hm_has_hs - hm_has_ha)
    hs_ha_only_set = (hs_has_ha - hs_has_hm) | (ha_has_hs - ha_has_hm)
    hm_ha_only_set = (hm_has_ha - hm_has_hs) | (ha_has_hm - ha_has_hs)
    all_three_set = (hs_has_hm & hs_has_ha) | (hm_has_hs & hm_has_ha) | (ha_has_hs & ha_has_hm)

    # For Venn diagram, count relationships not individual proteins
    # Use just one organism's perspective to avoid double counting
    venn_counts = {
        'hs_only':    len(hs_only_set),
        'hm_only':    len(hm_only_set),
        'ha_only':    len(ha_only_set),
        'hs_hm':      len(hs_has_hm - hs_has_ha),      # Count from HS perspective
        'hs_ha':      len(hs_has_ha - hs_has_hm),      # Count from HS perspective
        'hm_ha':      len(hm_has_ha - hm_has_hs),      # Count from HM perspective
        'all_three':  len(hs_has_hm & hs_has_ha),      # Count from HS perspective
    }

    # For CSV export, keep track of all proteins
    protein_sets = {
        'hs_only': hs_only_set,
        'hm_only': hm_only_set,
        'ha_only': ha_only_set,
        'hs_hm': hs_hm_only_set,
        'hs_ha': hs_ha_only_set,
        'hm_ha': hm_ha_only_set,
        'all_three': all_three_set,
    }

    # Sanity checks for all organisms
    hs_check = len(hs_only_set) + len(hs_has_hm - hs_has_ha) + len(hs_has_ha - hs_has_hm) + len(hs_has_hm & hs_has_ha)
    hm_check = len(hm_only_set) + len(hm_has_hs - hm_has_ha) + len(hm_has_ha - hm_has_hs) + len(hm_has_hs & hm_has_ha)
    ha_check = len(ha_only_set) + len(ha_has_hs - ha_has_hm) + len(ha_has_hm - ha_has_hs) + len(ha_has_hs & ha_has_hm)

    print(f"\n  Results:")
    print(f"    HS only:      {venn_counts['hs_only']}")
    print(f"    HM only:      {venn_counts['hm_only']}")
    print(f"    HA only:      {venn_counts['ha_only']}")
    print(f"    HS & HM:      {venn_counts['hs_hm']}")
    print(f"    HS & HA:      {venn_counts['hs_ha']}")
    print(f"    HM & HA:      {venn_counts['hm_ha']}")
    print(f"    Core (all 3): {venn_counts['all_three']}")
    
    print(f"\n  CSV Export Counts (includes proteins from all organisms in category):")
    print(f"    HS & HM proteins: {len(hs_hm_only_set)} ({len(hs_has_hm - hs_has_ha)} HS + {len(hm_has_hs - hm_has_ha)} HM)")
    print(f"    HS & HA proteins: {len(hs_ha_only_set)} ({len(hs_has_ha - hs_has_hm)} HS + {len(ha_has_hs - ha_has_hm)} HA)")
    print(f"    HM & HA proteins: {len(hm_ha_only_set)} ({len(hm_has_ha - hm_has_hs)} HM + {len(ha_has_hm - ha_has_hs)} HA)")
    print(f"    Core proteins: {len(all_three_set)}")
    
    print(f"\n  Sanity checks:")
    print(f"    HS: {hs_check} (expected {len(hs_set)}) {'✓' if hs_check == len(hs_set) else '✗'}")
    print(f"    HM: {hm_check} (expected {len(hm_set)}) {'✓' if hm_check == len(hm_set) else '✗'}")
    print(f"    HA: {ha_check} (expected {len(ha_set)}) {'✓' if ha_check == len(ha_set) else '✗'}")

    # Create Venn diagram with correct counts
    title = f"Unannotated Proteins (BLASTp RBH, ≥{MIN_IDENTITY}% identity)"
    create_venn(venn_counts, title, "venn_unannotated_blast")

    # Save protein lists to CSV files (includes all proteins from shared categories)
    save_protein_lists_to_csv(protein_sets, rbh_pairs, "unannotated_proteins")

    print("\nDone!")
    print("=" * 60)


if __name__ == '__main__':
    main()