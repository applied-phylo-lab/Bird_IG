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


def load_species_mapping(csv_path):
    """
    Load the summary_features.csv and create a mapping of haplotypes to species.
    
    Returns:
        dict: {species: [list of haplotypes]}
    """
    try:
        df = pd.read_csv(csv_path)
        
        # Try to find species column (common names: Species, species, Species_name, etc.)
        species_col = None
        for col in df.columns:
            if 'species' in col.lower():
                species_col = col
                break
        
        if species_col is None:
            raise ValueError("Could not find species column in CSV. Available columns: " + ", ".join(df.columns))
        
        # Try to find haplotype/sample column
        haplotype_col = None
        for col in df.columns:
            if any(term in col.lower() for term in ['haplotype']):
                haplotype_col = col
                break
        
        if haplotype_col is None:
            # Use first column as fallback
            haplotype_col = df.columns[0]
        
        # Group haplotypes by species
        species_map = defaultdict(list)
        for _, row in df.iterrows():
            species = row[species_col]
            haplotype = str(row[haplotype_col])
            species_map[species].append(haplotype)
        
        print(f"Loaded {len(species_map)} species from {csv_path}")
        for species, haplotypes in species_map.items():
            print(f"  {species}: {len(haplotypes)} haplotypes")
        
        return dict(species_map)
    
    except Exception as e:
        print(f"Error loading species mapping: {e}", file=sys.stderr)
        raise


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


def parse_lastz_maf(maf_file):
    """
    Parse Lastz .maf file to extract alignment information.
    
    MAF format:
    a score=...
    s ref_name start length strand total_length sequence
    s query_name start length strand total_length sequence
    
    Returns:
        list: List of alignment dictionaries
    """
    alignments = []
    
    try:
        with open(maf_file, 'r') as f:
            current_block = {}
            
            for line in f:
                line = line.strip()
                
                if line.startswith('a'):
                    # Alignment header
                    current_block = {}
                
                elif line.startswith('s'):
                    # Sequence line
                    parts = line.split()
                    if len(parts) >= 7:
                        name = parts[1]
                        start = int(parts[2])
                        length = int(parts[3])
                        strand = parts[4]
                        total_len = int(parts[5])
                        
                        if 'ref' not in current_block:
                            current_block['ref'] = {
                                'name': name,
                                'start': start,
                                'end': start + length,
                                'strand': strand,
                                'total_len': total_len
                            }
                        else:
                            current_block['query'] = {
                                'name': name,
                                'start': start,
                                'end': start + length,
                                'strand': strand,
                                'total_len': total_len
                            }
                            
                            # Complete alignment block
                            if 'ref' in current_block and 'query' in current_block:
                                ref = current_block['ref']
                                query = current_block['query']
                                
                                is_reverse = (ref['strand'] == '-' and query['strand'] == '+') or \
                                            (ref['strand'] == '+' and query['strand'] == '-')
                                
                                alignments.append({
                                    'ref_start': ref['start'],
                                    'ref_end': ref['end'],
                                    'query_start': query['start'],
                                    'query_end': query['end'],
                                    'ref_name': ref['name'],
                                    'query_name': query['name'],
                                    'is_reverse': is_reverse,
                                    'ref_strand': ref['strand'],
                                    'query_strand': query['strand'],
                                    'ref_len': ref['total_len'],
                                    'query_len': query['total_len']
                                })
                                current_block = {}
        
        return alignments
    
    except Exception as e:
        print(f"Error parsing Lastz MAF file {maf_file}: {e}", file=sys.stderr)
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


def find_alignment_files(directory, extensions):
    """
    Find all alignment files with given extensions in directory.
    
    Returns:
        dict: {file_path: (ref_name, query_name)}
    """
    files = {}
    directory = Path(directory)
    
    if not directory.exists():
        print(f"Warning: Directory {directory} does not exist", file=sys.stderr)
        return files
    
    for ext in extensions:
        for file_path in directory.rglob(f"*{ext}"):
            # Try to extract ref and query names from filename
            # Common patterns: ref_vs_query.ext, ref-query.ext, etc.
            stem = file_path.stem
            for sep in ['_vs_', '-vs-', '_', '-']:
                if sep in stem:
                    parts = stem.split(sep, 1)
                    if len(parts) == 2:
                        files[str(file_path)] = (parts[0], parts[1])
                        break
            else:
                # If no separator found, use filename as both
                files[str(file_path)] = (stem, stem)
    
    return files


def process_species_inversions(species_map, mummer_dir, lastz_dir, output_dir, min_length=1000, min_identity=80.0):
    """
    Process alignments for each species and detect inversions.
    """
    results = {
        'mummer': defaultdict(list),
        'lastz': defaultdict(list),
        'comparison': defaultdict(dict)
    }
    
    # Find Mummer alignment files
    print("\nFinding Mummer alignment files...")
    mummer_files = find_alignment_files(mummer_dir, ['.coords', '.delta'])
    print(f"Found {len(mummer_files)} Mummer alignment files")
    
    # Find Lastz alignment files
    print("\nFinding Lastz alignment files...")
    lastz_files = find_alignment_files(lastz_dir, ['.txt', '.maf', '.lav'])
    print(f"Found {len(lastz_files)} Lastz alignment files")
    
    # Process each species
    for species, haplotypes in species_map.items():
        print(f"\nProcessing species: {species}")
        print(f"  Haplotypes: {', '.join(haplotypes)}")
        
        # Find pairwise comparisons within species
        for i, haplo1 in enumerate(haplotypes):
            for haplo2 in haplotypes[i+1:]:
                print(f"  Comparing {haplo1} vs {haplo2}")
                
                # Process Mummer alignments
                mummer_inversions = []
                for file_path, (ref, query) in mummer_files.items():
                    # Check if this file matches our haplotypes (flexible matching)
                    # Try exact match first, then substring match
                    ref_match = (haplo1 == ref or haplo2 == ref or haplo1 in ref or haplo2 in ref)
                    query_match = (haplo1 == query or haplo2 == query or haplo1 in query or haplo2 in query)
                    
                    if ref_match and query_match and ref != query:
                        print(f"    Processing Mummer file: {Path(file_path).name}")
                        
                        if file_path.endswith('.delta'):
                            alignments = parse_mummer_delta(file_path)
                        else:
                            alignments = parse_mummer_coords(file_path)
                        
                        inversions = detect_inversions(alignments, min_length=min_length, min_identity=min_identity)
                        mummer_inversions.extend(inversions)
                
                # Process Lastz alignments
                lastz_inversions = []
                for file_path, (ref, query) in lastz_files.items():
                    # Check if this file matches our haplotypes (flexible matching)
                    ref_match = (haplo1 == ref or haplo2 == ref or haplo1 in ref or haplo2 in ref)
                    query_match = (haplo1 == query or haplo2 == query or haplo1 in query or haplo2 in query)
                    
                    if ref_match and query_match and ref != query:
                        print(f"    Processing Lastz file: {Path(file_path).name}")
                        
                        if file_path.endswith('.txt'):
                            alignments = parse_lastz_txt(file_path)
                        elif file_path.endswith('.maf'):
                            alignments = parse_lastz_maf(file_path)
                        else:
                            # For .lav files, would need different parser
                            print(f"    Warning: .lav format not yet supported, skipping")
                            continue
                        
                        inversions = detect_inversions(alignments, min_length=min_length, min_identity=min_identity)
                        lastz_inversions.extend(inversions)
                
                # Store results
                pair_key = f"{haplo1}_vs_{haplo2}"
                results['mummer'][species].extend(mummer_inversions)
                results['lastz'][species].extend(lastz_inversions)
                
                # Compare results
                results['comparison'][species][pair_key] = {
                    'mummer_count': len(mummer_inversions),
                    'lastz_count': len(lastz_inversions),
                    'mummer_inversions': mummer_inversions,
                    'lastz_inversions': lastz_inversions
                }
                
                print(f"    Found {len(mummer_inversions)} Mummer inversions, {len(lastz_inversions)} Lastz inversions")
    
    return results


def compare_results(results, output_dir):
    """
    Compare Mummer and Lastz inversion results and generate summary.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Create summary DataFrame
    summary_data = []
    
    for species, comparisons in results['comparison'].items():
        for pair_key, data in comparisons.items():
            mummer_inv = data['mummer_inversions']
            lastz_inv = data['lastz_inversions']
            
            # Calculate overlap (inversions detected by both methods)
            # Check if inversions overlap in both reference and query coordinates
            overlap_count = 0
            matched_lastz = set()  # Track which lastz inversions have been matched
            
            for m_inv in mummer_inv:
                for idx, l_inv in enumerate(lastz_inv):
                    if idx in matched_lastz:
                        continue
                    
                    # Check if inversions overlap in reference coordinates
                    m_ref_range = (min(m_inv['ref_start'], m_inv['ref_end']), 
                                  max(m_inv['ref_start'], m_inv['ref_end']))
                    l_ref_range = (min(l_inv['ref_start'], l_inv['ref_end']), 
                                  max(l_inv['ref_start'], l_inv['ref_end']))
                    
                    ref_overlap = (m_ref_range[0] <= l_ref_range[1] and 
                                  m_ref_range[1] >= l_ref_range[0])
                    
                    # Also check query coordinates if available
                    m_query_range = (min(m_inv.get('query_start', 0), m_inv.get('query_end', 0)), 
                                    max(m_inv.get('query_start', 0), m_inv.get('query_end', 0)))
                    l_query_range = (min(l_inv.get('query_start', 0), l_inv.get('query_end', 0)), 
                                    max(l_inv.get('query_start', 0), l_inv.get('query_end', 0)))
                    
                    query_overlap = (m_query_range[0] <= l_query_range[1] and 
                                    m_query_range[1] >= l_query_range[0])
                    
                    # Consider it an overlap if either ref or query coordinates overlap
                    # (or both if both are available)
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
    
    summary_df = pd.DataFrame(summary_data)
    summary_path = output_dir / 'inversion_summary.csv'
    summary_df.to_csv(summary_path, index=False)
    print(f"\nSummary saved to {summary_path}")
    
    # Save detailed results
    detailed_path = output_dir / 'inversion_details.json'
    # Convert to JSON-serializable format
    json_results = {}
    for species, comparisons in results['comparison'].items():
        json_results[species] = {}
        for pair_key, data in comparisons.items():
            json_results[species][pair_key] = {
                'mummer_count': data['mummer_count'],
                'lastz_count': data['lastz_count'],
                'mummer_inversions': data['mummer_inversions'],
                'lastz_inversions': data['lastz_inversions']
            }
    
    with open(detailed_path, 'w') as f:
        json.dump(json_results, f, indent=2)
    print(f"Detailed results saved to {detailed_path}")
    
    # Print summary statistics
    print("\n" + "="*80)
    print("SUMMARY STATISTICS")
    print("="*80)
    print(f"\nTotal species analyzed: {len(results['comparison'])}")
    print(f"Total pairwise comparisons: {len(summary_data)}")
    print(f"\nTotal Mummer inversions: {summary_df['Mummer_Inversions'].sum()}")
    print(f"Total Lastz inversions: {summary_df['Lastz_Inversions'].sum()}")
    print(f"Total overlapping inversions: {summary_df['Overlap_Count'].sum()}")
    print(f"\nMummer-only inversions: {summary_df['Mummer_Only'].sum()}")
    print(f"Lastz-only inversions: {summary_df['Lastz_Only'].sum()}")
    
    return summary_df


def main():
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
    
    args = parser.parse_args()
    
    print("="*80)
    print("INVERSION DETECTION SCRIPT")
    print("="*80)
    print(f"\nSummary CSV: {args.summary_csv}")
    print(f"Mummer directory: {args.mummer_dir}")
    print(f"Lastz directory: {args.lastz_dir}")
    print(f"Output directory: {args.output_dir}")
    print(f"Minimum inversion length: {args.min_length} bp")
    print(f"Minimum identity: {args.min_identity}%")
    
    # Load species mapping
    species_map = load_species_mapping(args.summary_csv)
    
    # Process inversions
    results = process_species_inversions(
        species_map,
        args.mummer_dir,
        args.lastz_dir,
        args.output_dir,
        args.min_length,
        args.min_identity
    )
    
    # Compare and save results
    summary_df = compare_results(results, args.output_dir)
    
    print("\n" + "="*80)
    print("Analysis complete!")
    print("="*80)


if __name__ == '__main__':
    main()
