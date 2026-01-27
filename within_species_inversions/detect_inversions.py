#!/usr/bin/env python3
"""
Script to detect inversions between birds of the same species using Mummer and Lastz alignments.

This script:
1. Reads summary_features.csv to determine which haplotypes belong to the same species
2. Processes Mummer alignments to detect inversions
3. Processes Lastz alignments to detect inversions
4. Compares results from both methods
"""

import os
import sys
import pandas as pd
import argparse
from pathlib import Path
from collections import defaultdict
import json
import itertools
from concurrent.futures import ProcessPoolExecutor, as_completed



def parse_args():
    parser = argparse.ArgumentParser(
        description='Detect inversions between birds of the same species using Mummer and Lastz alignments'
    )
    parser.add_argument(
        '--summary-csv', '-s',
        required=True,
        help='Path to summary_features.csv file'
    )
    parser.add_argument(
        '--mummer-dir', '-m',
        required=True,
        help='Directory containing Mummer alignment files'
    )
    parser.add_argument(
        '--lastz-dir', '-l',
        required=True,
        help='Directory containing Lastz alignment files'
    )
    parser.add_argument(
        '--output-dir', '-o',
        default='./inversion_results',
        help='Output directory for results (default: ./inversion_results)'
    )
    parser.add_argument(
        '--min-length', '-L',
        type=int,
        default=1000,
        help='Minimum inversion length in bp (default: 1000)'
    )
    parser.add_argument(
        '--min-identity', '-i',
        type=float,
        default=80.0,
        help='Minimum alignment identity percentage (default: 80.0)'
    )
    parser.add_argument("-c", "--cores", type=int, default=1,
                        help="Number of parallel processes (default: 1)")
    
    return parser.parse_args()

def process_pair(args_tuple,min_length=1000, min_identity=80.0):
    (
        args,
        order,
        species,
        h1,
        h2
    ) = args_tuple

    results = {
        'mummer': defaultdict(list),
        'lastz': defaultdict(list),
        'comparison': defaultdict(dict)
    }

    # mummer alignments
    mummer_prefix = os.path.join(
        args.mummer_dir,
        f"{species}_{h1['Haplotype']}_vs_{h2['Haplotype']}"
    )

    delta = mummer_prefix + ".delta"
    filtered_delta = mummer_prefix + ".filtered.delta"
    coords = mummer_prefix + ".coords"

    if os.path.exists(delta):
        # Process Mummer alignments
        mummer_inversions = []
        alignments = parse_mummer_delta(delta)
        inversions = detect_inversions(alignments, min_length=min_length, min_identity=min_identity)
        mummer_inversions.extend(inversions)
    
    # lastz alignments
    lastz_prefix = os.path.join(
        args.lastz_dir,
        f"{h1['Haplotype']}_{h2['Haplotype']}"
    )

    txt = lastz_prefix + ".txt"
    lastz_inversions = []
    alignments = parse_lastz_txt(txt)
    inversions = detect_inversions(alignments, min_length=min_length, min_identity=min_identity)
    lastz_inversions.extend(inversions)
    pair_key = f"{h1['Haplotype']}_vs_{h2['Haplotype']}"
    results['mummer'][species].extend(mummer_inversions)
    results['lastz'][species].extend(lastz_inversions)

    results['comparison'][species][pair_key] = {
                    'mummer_count': len(mummer_inversions),
                    'lastz_count': len(lastz_inversions),
                    'mummer_inversions': mummer_inversions,
                    'lastz_inversions': lastz_inversions
                }

    return results
    
def parse_mummer_coords(coords_file):
    """
    Parse Mummer .coords file to extract alignment information.
    
    Format: [S1] [E1] [S2] [E2] [LEN1] [LEN2] [%IDY] [LENR] [LENQ] [COVR] [COVQ] [TAGS]
    
    Note: In coords files, if s2 > e2, the query is on the reverse strand (inversion).
    
    Returns:
        list: List of alignment dictionaries with is_reverse flag
    """
    alignments = []
    
    try:
        with open(coords_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('=') or line.startswith('['):
                    continue
                
                parts = line.split()
                if len(parts) < 12:
                    continue
                
                try:
                    s1 = int(parts[0])
                    e1 = int(parts[1])
                    s2 = int(parts[2])
                    e2 = int(parts[3])
                    len1 = int(parts[4])
                    len2 = int(parts[5])
                    pct_idy = float(parts[6])
                    lenr = int(parts[7])
                    lenq = int(parts[8])
                    
                    # Check if query is on reverse strand (s2 > e2 indicates reverse)
                    is_reverse = s2 > e2
                    
                    # Normalize coordinates (always use smaller value as start)
                    if is_reverse:
                        query_start = min(s2, e2)
                        query_end = max(s2, e2)
                    else:
                        query_start = s2
                        query_end = e2
                    
                    # Extract query and reference names from tags
                    if len(parts) >= 12:
                        ref_tag = parts[10]
                        query_tag = parts[11]
                    else:
                        # Try to extract from filename or header
                        ref_tag = "unknown"
                        query_tag = "unknown"
                    
                    alignments.append({
                        'ref_start': s1,
                        'ref_end': e1,
                        'query_start': query_start,
                        'query_end': query_end,
                        'ref_len': lenr,
                        'query_len': lenq,
                        'identity': pct_idy,
                        'ref_name': ref_tag,
                        'query_name': query_tag,
                        'is_reverse': is_reverse
                    })
                except (ValueError, IndexError) as e:
                    continue
        
        return alignments
    
    except Exception as e:
        print(f"Error parsing Mummer coords file {coords_file}: {e}", file=sys.stderr)
        return []


def parse_mummer_delta(delta_file):
    """
    Parse Mummer .delta file to extract alignment information including strand orientation.
    
    Returns:
        list: List of alignment dictionaries with strand information
    """
    alignments = []
    
    try:
        with open(delta_file, 'r') as f:
            current_ref = None
            current_query = None
            ref_len = None
            query_len = None
            
            for line in f:
                line = line.strip()
                
                if line.startswith('>'):
                    # Header line: >ref_name query_name ref_len query_len
                    parts = line[1:].split()
                    if len(parts) >= 4:
                        current_ref = parts[0]
                        current_query = parts[1]
                        ref_len = int(parts[2])
                        query_len = int(parts[3])
                
                elif line and not line.startswith('NUCMER') and not line.startswith('PROMER'):
                    # Alignment data line
                    parts = line.split()
                    if len(parts) >= 7:
                        try:
                            ref_start = int(parts[0])
                            ref_end = int(parts[1])
                            query_start = int(parts[2])
                            query_end = int(parts[3])
                            
                            # In delta format, if query_start > query_end, it's on reverse strand
                            is_reverse = query_start > query_end
                            
                            if is_reverse:
                                # Swap to get correct coordinates
                                query_start, query_end = query_end, query_start
                            
                            alignments.append({
                                'ref_start': ref_start,
                                'ref_end': ref_end,
                                'query_start': query_start,
                                'query_end': query_end,
                                'ref_name': current_ref,
                                'query_name': current_query,
                                'is_reverse': is_reverse,
                                'ref_len': ref_len,
                                'query_len': query_len
                            })
                        except (ValueError, IndexError):
                            continue
        
        return alignments
    
    except Exception as e:
        print(f"Error parsing Mummer delta file {delta_file}: {e}", file=sys.stderr)
        return []


def parse_lastz_txt(txt_file):
    """
    Parse Lastz .txt file to extract alignment information.
    
    Format: #name1	strand1	start1	end1	length1	name2	strand2	start2+	end2+	length2	id%
    
    Returns:
        list: List of alignment dictionaries
    """
    alignments = []
    
    try:
        with open(txt_file, 'r') as f:
            for line in f:
                line = line.strip()
                
                # Skip header line
                if line.startswith('#'):
                    continue
                
                if not line:
                    continue
                
                # Parse tab-separated values
                parts = line.split('\t')
                if len(parts) < 11:
                    continue
                
                try:
                    name1 = parts[0]
                    strand1 = parts[1]
                    start1 = int(parts[2])
                    end1 = int(parts[3])
                    length1 = int(parts[4])
                    name2 = parts[5]
                    strand2 = parts[6]
                    start2 = int(parts[7])
                    end2 = int(parts[8])
                    length2 = int(parts[9])
                    identity = float(parts[10].rstrip('%'))
                    
                    # Determine if this is a reverse alignment (inversion)
                    # If strands are different, it's a reverse alignment
                    is_reverse = (strand1 == '+' and strand2 == '-') or (strand1 == '-' and strand2 == '+')
                    
                    alignments.append({
                        'ref_start': start1,
                        'ref_end': end1,
                        'query_start': start2,
                        'query_end': end2,
                        'ref_name': name1,
                        'query_name': name2,
                        'is_reverse': is_reverse,
                        'ref_strand': strand1,
                        'query_strand': strand2,
                        'ref_len': length1,
                        'query_len': length2,
                        'identity': identity
                    })
                except (ValueError, IndexError) as e:
                    continue
        
        return alignments
    
    except Exception as e:
        print(f"Error parsing Lastz txt file {txt_file}: {e}", file=sys.stderr)
        return []

def detect_inversions(alignments, min_length=1000, min_identity=80.0):
    """
    Detect inversions from alignment data.
    
    An inversion is detected when:
    1. The query alignment is on the reverse strand (is_reverse=True)
    2. The alignment is sufficiently long
    3. The alignment has sufficient identity
    
    Returns:
        list: List of inversion dictionaries
    """
    inversions = []
    
    for aln in alignments:
        # Check if this is a reverse alignment (potential inversion)
        is_reverse = aln.get('is_reverse', False)
        
        if not is_reverse:
            continue
        
        # Calculate alignment length
        ref_len = abs(aln['ref_end'] - aln['ref_start'])
        query_len = abs(aln['query_end'] - aln['query_start'])
        aln_length = min(ref_len, query_len)
        
        # Check minimum length
        if aln_length < min_length:
            continue
        
        # Check identity if available
        identity = aln.get('identity', 100.0)
        if identity < min_identity:
            continue
        
        inversions.append({
            'ref_name': aln['ref_name'],
            'query_name': aln['query_name'],
            'ref_start': aln['ref_start'],
            'ref_end': aln['ref_end'],
            'query_start': aln['query_start'],
            'query_end': aln['query_end'],
            'length': aln_length,
            'identity': identity,
            'ref_len': aln.get('ref_len', 0),
            'query_len': aln.get('query_len', 0)
        })
    
    return inversions



def compare_results(species, comparisons, output_dir):
    """
    Compare Mummer and Lastz inversion results for ONE species.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    summary_data = []

    for pair_key, data in comparisons.items():
        mummer_inv = data['mummer_inversions']
        lastz_inv = data['lastz_inversions']

        overlap_count = 0
        matched_lastz = set()

        for m_inv in mummer_inv:
            for idx, l_inv in enumerate(lastz_inv):
                if idx in matched_lastz:
                    continue

                # Reference overlap
                m_ref_range = (
                    min(m_inv['ref_start'], m_inv['ref_end']),
                    max(m_inv['ref_start'], m_inv['ref_end'])
                )
                l_ref_range = (
                    min(l_inv['ref_start'], l_inv['ref_end']),
                    max(l_inv['ref_start'], l_inv['ref_end'])
                )

                ref_overlap = (
                    m_ref_range[0] <= l_ref_range[1] and
                    m_ref_range[1] >= l_ref_range[0]
                )

                # Query overlap (optional)
                m_query_range = (
                    min(m_inv.get('query_start', 0), m_inv.get('query_end', 0)),
                    max(m_inv.get('query_start', 0), m_inv.get('query_end', 0))
                )
                l_query_range = (
                    min(l_inv.get('query_start', 0), l_inv.get('query_end', 0)),
                    max(l_inv.get('query_start', 0), l_inv.get('query_end', 0))
                )

                query_overlap = (
                    m_query_range[0] <= l_query_range[1] and
                    m_query_range[1] >= l_query_range[0]
                )

                if ref_overlap and (not m_query_range[0] or query_overlap):
                    overlap_count += 1
                    matched_lastz.add(idx)
                    break

        summary_data.append({
            'Species': species,
            'Pair': pair_key,
            'Mummer_Inversions': data['mummer_count'],
            'Lastz_Inversions': data['lastz_count'],
            'Overlap_Count': overlap_count,
            'Mummer_Only': data['mummer_count'] - overlap_count,
            'Lastz_Only': data['lastz_count'] - overlap_count
        })

    return pd.DataFrame(summary_data)

def main():
    args = parse_args()
    
    print("="*80)
    print("INVERSION DETECTION SCRIPT")
    print("="*80)
    print(f"\nSummary CSV: {args.summary_csv}")
    print(f"Mummer directory: {args.mummer_dir}")
    print(f"Lastz directory: {args.lastz_dir}")
    print(f"Output directory: {args.output_dir}")
    os.makedirs(args.output_dir, exist_ok=True)
    print(f"Minimum inversion length: {args.min_length} bp")
    print(f"Minimum identity: {args.min_identity}%")
    
    # Load species mapping
    df = pd.read_csv(args.summary_csv, sep="\t")

    jobs = []
    for order, order_df in df.groupby("Order"):
        for species, sp_df in order_df.groupby("Species"):
            haplotypes = sp_df.to_dict(orient="records")
            for h1, h2 in itertools.combinations(haplotypes, 2):
                jobs.append((args, order, species, h1, h2))

    print(f"Total pairwise jobs: {len(jobs)}")
    print(f"Using {args.cores} cores")

    with ProcessPoolExecutor(max_workers=args.cores) as executor:
        futures = [executor.submit(process_pair, job) for job in jobs]
        total_results = []
        # save results to tsv
        for f in as_completed(futures):
            results = f.result()
            total_results.append(results)
    

    species_results = defaultdict(dict)

    for res in total_results:
        for species, comparisons in res["comparison"].items():
            species_results[species].update(comparisons)

    for species, comparisons in species_results.items():
        # Apply your function
        result_df = compare_results(species,comparisons, args.output_dir)

        # Write species-specific TSV
        out_tsv = os.path.join(args.output_dir, f"{species}_comparison.tsv")

        result_df.to_csv(out_tsv, sep="\t", index=False)
    

    print("All jobs completed.")


if __name__ == '__main__':
    main()
