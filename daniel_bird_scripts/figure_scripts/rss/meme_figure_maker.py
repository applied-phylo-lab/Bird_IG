import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.path import Path
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties
import numpy as np
import sys
import os

# Get input file from command line argument
if len(sys.argv) > 1:
    input_file = sys.argv[1]
else:
    print("Usage: python script.py <input_csv_file>")
    print("Example: python script.py heptamer_list.csv")
    print("         python script.py IGH_heptamer_analysis.csv")
    sys.exit(1)

if not os.path.exists(input_file):
    print(f"Error: File '{input_file}' not found")
    sys.exit(1)

# Read the data file
seq_df = pd.read_csv(input_file)

# Detect sequence type and length from the column name
if 'Heptamer' in seq_df.columns:
    seq_type = 'Heptamer'
    seq_type_lower = 'heptamer'
    seq_length = 7
elif 'Nonamer' in seq_df.columns:
    seq_type = 'Nonamer'
    seq_type_lower = 'nonamer'
    seq_length = 9
else:
    print("Error: CSV must contain either 'Heptamer' or 'Nonamer' column")
    sys.exit(1)

# Detect the weight column (look for various occurrence column names)
weight_column = None
possible_weight_columns = ['Haplotypes', 'Number of Haplotypes']
for col in possible_weight_columns:
    if col in seq_df.columns:
        weight_column = col
        break

if weight_column is None:
    print(f"Error: Could not find occurrence/count column. Available columns: {list(seq_df.columns)}")
    sys.exit(1)

print(f"Detected {seq_type} data (sequence length: {seq_length})")
print(f"Using weight column: '{weight_column}'")

# ============================================================================
# CREATE MEME SEQUENCE LOGO
# ============================================================================

# Create figure
fig = plt.figure(figsize=(14, 4))
ax = fig.add_subplot(111)

bases = ['A', 'C', 'G', 'T']
pwm = {base: [0] * seq_length for base in bases}

# Build position weight matrix weighted by occurrences
for _, row in seq_df.iterrows():
    sequence = row[seq_type].upper()
    weight = row[weight_column]
    for pos, nucleotide in enumerate(sequence):
        if nucleotide in bases:
            pwm[nucleotide][pos] += weight

# Convert counts to frequencies
total_counts = [sum(pwm[base][pos] for base in bases) for pos in range(seq_length)]
ppm = {base: [pwm[base][pos] / total_counts[pos] if total_counts[pos] > 0 else 0 
              for pos in range(seq_length)] for base in bases}

ax.set_xlim(0, seq_length)
ax.set_ylim(0, 1)
ax.axis('off')

# Colors for bases
colors = {'A': '#F8CD9C', 'C': '#172869', 'G': '#088BBE', 'T': '#EA7580'}
fp = FontProperties(family='monospace', weight='bold')

# Draw sequence logo
for pos in range(seq_length):
    freqs = [(base, ppm[base][pos]) for base in bases]
    freqs.sort(key=lambda x: x[1])
    
    y_offset = 0
    for base, freq in freqs:
        if freq > 0.01:
            letter_height = freq
            text_path = TextPath((0, 0), base, size=1, prop=fp)
            bbox = text_path.get_extents()
            
            width_scale = 0.85
            height_scale = letter_height / (bbox.height)
            x_offset = pos + 0.5 - (bbox.width * width_scale / 2)
            
            vertices = text_path.vertices
            vertices_scaled = vertices.copy()
            vertices_scaled[:, 0] = vertices[:, 0] * width_scale + x_offset
            vertices_scaled[:, 1] = vertices[:, 1] * height_scale + y_offset
            
            path = Path(vertices_scaled, text_path.codes)
            patch = patches.PathPatch(path, facecolor=colors[base], 
                                     edgecolor='black', linewidth=1.5)
            ax.add_patch(patch)
            
            y_offset += letter_height

plt.tight_layout(pad=0)

# Create output filename based on input filename
input_basename = os.path.splitext(os.path.basename(input_file))[0]
output_filename = input_basename+"_meme_figure.svg"
plt.savefig(os.path.dirname(input_file)+"/"+output_filename, dpi=300, bbox_inches='tight', facecolor='white')

print(f"\nMEME sequence logo created successfully: {output_filename}")
print(f"Total {weight_column.lower()}: {seq_df[weight_column].sum()}")
print(f"\nMost frequent base at each position:")
for pos in range(seq_length):
    max_base = max(bases, key=lambda b: ppm[b][pos])
    max_freq = ppm[max_base][pos]
    print(f"Position {pos+1}: {max_base} ({max_freq:.1%})")