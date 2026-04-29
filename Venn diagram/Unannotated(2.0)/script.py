"""
UNANNOTATED PROTEIN BLAST ANALYSIS
====================================
Dependencies:
  pip install matplotlib matplotlib-venn
  BLAST+ (https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/)

This script performs reciprocal best hit (RBH) analysis on unannotated proteins
from three strains (HS, HM, HA) and generates Venn diagrams and CSV files.
Venn diagram counts represent ortholog groups (connected components of RBH pairs).
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
    
    Uses a temporary directory for the BLAST database to avoid issues
    with spaces in file paths (BLAST+ doesn't handle them well).
    
    Args:
        combined_fasta: Combined FASTA file with all proteins
        output_file: Output file for BLAST results
        evalue: E-value cutoff
        threads: Number of threads to use
    """
    import tempfile, shutil

    # Use a temp directory for the BLAST DB to avoid space-in-path issues
    tmp_dir = Path(tempfile.mkdtemp(prefix='blast_'))
    db_name = str(tmp_dir / 'unannotated_db')
    # Also copy the input FASTA to the temp dir for safety
    tmp_fasta = tmp_dir / 'query.faa'
    shutil.copy2(str(combined_fasta), str(tmp_fasta))

    print(f"  Using temp directory: {tmp_dir}")
    print("  Creating BLAST database...")
    try:
        subprocess.run([
            'makeblastdb',
            '-in', str(tmp_fasta),
            '-dbtype', 'prot',
            '-out', db_name
        ], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: makeblastdb failed: {e.stderr}")
        shutil.rmtree(tmp_dir, ignore_errors=True)
        sys.exit(1)

    tmp_output = tmp_dir / 'blast_results.txt'
    print(f"  Running BLASTp (e-value < {evalue}, {threads} threads)...")
    print("  This may take several minutes...")
    try:
        subprocess.run([
            'blastp',
            '-query', str(tmp_fasta),
            '-db', db_name,
            '-out', str(tmp_output),
            '-outfmt', '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs',
            '-evalue', str(evalue),
            '-max_target_seqs', '10',
            '-num_threads', str(threads)
        ], check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as e:
        print(f"ERROR: blastp failed: {e.stderr}")
        shutil.rmtree(tmp_dir, ignore_errors=True)
        sys.exit(1)

    # Copy results back to the original output location
    shutil.copy2(str(tmp_output), str(output_file))
    
    # Clean up temp directory
    shutil.rmtree(tmp_dir, ignore_errors=True)

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


def build_ortholog_groups(rbh_pairs, all_proteins):
    """
    Build ortholog groups (clusters) from RBH pairs using connected components.
    
    Each group contains proteins from one or more organisms that are connected
    by RBH relationships. This ensures consistent counting across all organisms.
    
    Args:
        rbh_pairs: Set of RBH pairs (tuples of protein IDs)
        all_proteins: Dict mapping organism tag to list of protein IDs
        
    Returns:
        Tuple of (groups, protein_to_group) where:
          - groups: dict mapping group_id to set of protein IDs
          - protein_to_group: dict mapping protein_id to group_id
    """
    # Union-Find data structure for clustering
    parent = {}
    
    def find(x):
        if x not in parent:
            parent[x] = x
        while parent[x] != x:
            parent[x] = parent[parent[x]]  # path compression
            x = parent[x]
        return x
    
    def union(a, b):
        ra, rb = find(a), find(b)
        if ra != rb:
            parent[ra] = rb
    
    # Initialize all proteins
    for org, proteins in all_proteins.items():
        for p in proteins:
            find(p)
    
    # Merge RBH pairs
    for p1, p2 in rbh_pairs:
        union(p1, p2)
    
    # Build groups
    groups = defaultdict(set)
    for p in parent:
        groups[find(p)].add(p)
    
    # Map protein to group
    protein_to_group = {}
    for gid, members in groups.items():
        for p in members:
            protein_to_group[p] = gid
    
    return dict(groups), protein_to_group


def classify_by_ortholog_groups(groups, protein_to_group):
    """
    Classify ortholog groups by which organisms they span and assign
    proteins to Venn categories accordingly.
    """
    category_names = ['hs_only', 'hm_only', 'ha_only', 'hs_hm', 'hs_ha', 'hm_ha', 'all_three']
    
    org_coverage_to_category = {
        frozenset(['HS']):             'hs_only',
        frozenset(['HM']):             'hm_only',
        frozenset(['HA']):             'ha_only',
        frozenset(['HS', 'HM']):       'hs_hm',
        frozenset(['HS', 'HA']):       'hs_ha',
        frozenset(['HM', 'HA']):       'hm_ha',
        frozenset(['HS', 'HM', 'HA']): 'all_three',
    }
    
    # Classify each GROUP by its organism coverage
    group_category = {}     # group_id → category
    group_sets = {c: set() for c in category_names}  # category → set of group_ids
    
    for gid, members in groups.items():
        orgs = frozenset(p.split('|')[0] for p in members)
        category = org_coverage_to_category.get(orgs)
        if category:
            group_category[gid] = category
            group_sets[category].add(gid)
    
    # Assign each PROTEIN to its group's category
    protein_sets = {c: set() for c in category_names}
    for gid, members in groups.items():
        cat = group_category.get(gid)
        if cat:
            for p in members:
                protein_sets[cat].add(p)
    
    # Venn counts = number of ortholog groups per region
    venn_counts = {c: len(group_sets[c]) for c in category_names}
    
    # Track paralogs (groups with >1 protein from same organism)
    paralog_info = {org: [] for org in ['HS', 'HM', 'HA']}
    for gid, members in groups.items():
        org_members = defaultdict(list)
        for p in members:
            org_members[p.split('|')[0]].append(p)
        for org, plist in org_members.items():
            if len(plist) > 1:
                cat = group_category.get(gid, 'unclassified')
                paralog_info[org].append({
                    'group_id': gid,
                    'category': cat,
                    'proteins': plist,
                    'count': len(plist)
                })
    
    return venn_counts, protein_sets, paralog_info


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
    plt.savefig(SCRIPT_DIR / f'{output_name}.png', dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"  Saved: {output_name}.svg, {output_name}.tiff, {output_name}.png")


def save_protein_lists_to_csv(protein_sets, rbh_pairs, output_prefix):
    """
    Save protein lists for each category to separate CSV files.
    
    Args:
        protein_sets: Dictionary containing sets of proteins for each category
        rbh_pairs: Set of RBH pairs for finding partner proteins
        output_prefix: Prefix for output CSV files
    """
    print("\nSTEP 7: Saving protein lists to CSV files...")
    
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
        writer.writerow(['Category', 'Ortholog_Groups', 'Total_Proteins', 'HS_Proteins', 'HM_Proteins', 'HA_Proteins', 'Description'])
        
        for key, label, desc in [
            ('hs_only',    'HS_only',        'Unique to HS strain'),
            ('hm_only',    'HM_only',        'Unique to HM strain'),
            ('ha_only',    'HA_only',        'Unique to HA strain'),
            ('hs_hm',     'HS_and_HM_only', 'Shared between HS and HM only'),
            ('hs_ha',     'HS_and_HA_only', 'Shared between HS and HA only'),
            ('hm_ha',     'HM_and_HA_only', 'Shared between HM and HA only'),
            ('all_three', 'Core_all_three', 'Shared among all three strains'),
        ]:
            proteins = protein_sets[key]
            hs_count = sum(1 for p in proteins if p.split('|')[0] == 'HS')
            hm_count = sum(1 for p in proteins if p.split('|')[0] == 'HM')
            ha_count = sum(1 for p in proteins if p.split('|')[0] == 'HA')
            # For single-org categories, groups = proteins; for shared, need to count groups
            # We approximate groups as max per-org count (since RBH is 1-to-1 within org pairs)
            org_counts = [c for c in [hs_count, hm_count, ha_count] if c > 0]
            n_groups = max(org_counts) if org_counts else 0
            writer.writerow([label, n_groups, len(proteins), hs_count, hm_count, ha_count, desc])
        
        total_proteins = sum(len(v) for v in protein_sets.values())
        writer.writerow(['Total', '', total_proteins, '', '', '', 'Total unannotated proteins analyzed'])
    
    print(f"  Saved summary to: {summary_file.name}")


def main():
    print("=" * 60)
    print("  UNANNOTATED PROTEIN BLAST ANALYSIS (v2 - fixed)")
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
    for org in ['HS', 'HM', 'HA']:
        expected = len(unannotated_ids[org])
        actual = len(all_proteins[org])
        status = '✓' if actual == expected else f'✗ (expected {expected})'
        print(f"    {org}: {actual} extracted {status}")

    # STEP 3: Run BLAST
    print("\nSTEP 3: Running BLASTp all-vs-all...")
    blast_output = SCRIPT_DIR / 'blast_results.txt'
    run_blast(combined, blast_output, EVALUE_CUTOFF, NUM_THREADS)

    # STEP 4: Find RBH
    print(f"\nSTEP 4: Finding Reciprocal Best Hits (>={MIN_IDENTITY}% identity, >={MIN_QUERY_COV}% coverage)...")
    rbh_pairs = find_rbh(blast_output, MIN_IDENTITY, MIN_QUERY_COV)
    print(f"  RBH pairs found: {len(rbh_pairs)}")

    # STEP 5: Build ortholog groups
    print("\nSTEP 5: Building ortholog groups from RBH pairs...")
    groups, protein_to_group = build_ortholog_groups(rbh_pairs, all_proteins)
    
    # Count group types
    group_type_counts = defaultdict(int)
    for gid, members in groups.items():
        orgs = frozenset(p.split('|')[0] for p in members)
        group_type_counts[orgs] += 1
    
    print(f"  Total ortholog groups: {len(groups)}")
    for orgs, count in sorted(group_type_counts.items(), key=lambda x: (len(x[0]), x[0])):
        print(f"    {' & '.join(sorted(orgs))}: {count} groups")

    # STEP 6: Classify proteins and create Venn diagram
    print("\nSTEP 6: Classifying proteins and creating Venn diagram...")
    venn_counts, protein_sets, paralog_info = classify_by_ortholog_groups(groups, protein_to_group)

    print(f"\n  Ortholog Group Counts (Venn diagram values):")
    print(f"    HS only:      {venn_counts['hs_only']}")
    print(f"    HM only:      {venn_counts['hm_only']}")
    print(f"    HA only:      {venn_counts['ha_only']}")
    print(f"    HS & HM:      {venn_counts['hs_hm']}")
    print(f"    HS & HA:      {venn_counts['hs_ha']}")
    print(f"    HM & HA:      {venn_counts['hm_ha']}")
    print(f"    Core (all 3): {venn_counts['all_three']}")

    # Sanity checks
    print(f"\n  Per-organism protein totals:")
    for org in ['HS', 'HM', 'HA']:
        org_total = sum(1 for key in protein_sets for p in protein_sets[key] if p.split('|')[0] == org)
        expected = len(all_proteins[org])
        status = '✓' if org_total == expected else '✗'
        print(f"    {org}: {org_total} (expected {expected}) {status}")

    # Create Venn diagram
    title = f"Unannotated Proteins (BLASTp RBH, ≥{MIN_IDENTITY}% identity)"
    create_venn(venn_counts, title, "venn_unannotated_blast")

    # Save protein lists to CSV files
    save_protein_lists_to_csv(protein_sets, rbh_pairs, "unannotated_proteins")

    print("\nDone!")
    print("=" * 60)


if __name__ == '__main__':
    main()