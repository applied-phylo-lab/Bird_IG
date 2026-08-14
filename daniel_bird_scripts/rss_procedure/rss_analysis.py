import os
import csv
import sys
import subprocess
from collections import Counter

locus=sys.argv[1]
bird_data=["/local/storage/kav67/clean_birds"]
input_file = sys.argv[2]
rss_data=[]
with open(input_file,"r") as read:
    reader=csv.reader(read)
    header=next(reader)
    for row in reader:
        rss_data.append(row)
    read.close()
heptamers=[]
nonamers=[]
for r in rss_data:
    if [r[9],0,0,0,0,0,0,0] not in heptamers:
        heptamers.append([r[9],0,0,0,0,0,0,0]) #[seq,less25,less0,more25,more0,occur,species,hap]
    if [r[10],0,0,0,0,0,0,0] not in nonamers:
        nonamers.append([r[10],0,0,0,0,0,0,0])

type=sys.argv[3] #-h for heptamers, -n for nonamers
def analyze_rss_positions(gene_list):
    # Filter genes that have an RSS (not False)
    rss_genes = [(i, gene) for i, gene in enumerate(gene_list) if gene[7] is not False and gene[7] != "False" and gene[7] is not None]
    
    # If no genes with RSS, return empty strings
    if not rss_genes:
        return ["", "", "", "", "", "", "", ""]
    
    total_genes = len(gene_list)
    first_25_threshold = total_genes * 0.25
    last_25_threshold = total_genes * 0.75
    
    # First and last RSS gene with their original indices
    first_rss_index, first_rss_gene = rss_genes[0]
    last_rss_index, last_rss_gene = rss_genes[-1]
    
    # Check conditions for first RSS gene
    first_rss_seq = first_rss_gene[7]
    first_rss_position = first_rss_gene[2]
    first_is_first_gene = first_rss_index == 0
    first_in_first_25 = first_is_first_gene or first_rss_index < first_25_threshold

    # Check conditions for last RSS gene
    last_rss_seq = last_rss_gene[7]
    last_rss_position = last_rss_gene[2]
    last_is_last_gene = last_rss_index == total_genes - 1
    last_in_last_25 = last_is_last_gene or last_rss_index >= last_25_threshold
    
    return [
        first_rss_seq,
        first_rss_position,
        first_is_first_gene,
        first_in_first_25,
        last_rss_seq,
        last_rss_position,
        last_is_last_gene,
        last_in_last_25
    ]

def is_low_complexity(seq, threshold=0.7):
    """
    Returns True if the sequence is considered low complexity.
    A sequence is low complexity if the most frequent base or the sum
    of the two most frequent bases exceed the threshold proportion.
    """
    seq = str(seq).upper()
    if len(seq) == 0:
        return True

    counts = Counter(seq)
    freqs = sorted(counts.values(), reverse=True)
    
    if freqs[0] / len(seq) >= threshold:
        return True
    if len(freqs) > 1 and sum(freqs[:2]) / len(seq) >= threshold:
        return True
    return False

order_names=["Cormorants","Cranes","Doves","Eagles","Falcons","finches","house_finches","Hummingbirds","Ibises","Landfowl","MiscBirds","Owls","Parrots","Plovers","Songbirds","Suboscines","Waterfowl","Woodpeckers"] #names of all bird order folders

for q in bird_data:
    for f in os.listdir(q):
        if f in order_names and os.path.isdir(q+"/"+f):
            for f1 in os.listdir(q+"/"+f):
                if os.path.isdir(q+"/"+f+"/"+f1):
                    species_rss=[]
                    for f2 in os.listdir(q+"/"+f+"/"+f1):
                        if os.path.isdir(q+"/"+f+"/"+f1+"/"+f2):
                            for f3 in os.listdir(q+"/"+f+"/"+f1+"/"+f2):
                                if f3=="combined_genes_"+locus+".txt":
                                    genes=[]
                                    with open(q+"/"+f+"/"+f1+"/"+f2+"/"+f3,"r") as genes_read:
                                        reader=csv.reader(genes_read,delimiter="\t")
                                        header=next(reader)
                                        for row in reader:
                                            if len(row[4])>250 and is_low_complexity(row[4])==False: #and str(row[5])=="True":
                                                genes.append(row)
                                        genes_read.close()
                                    for g in genes:
                                        found=False
                                        for r in rss_data:
                                            if g[0:5]==r[1:6]:
                                                found=True
                                                if type=="-h":
                                                    g.append(r[9])
                                                elif type=="-n":
                                                    g.append(r[10])
                                        if found==False:
                                            g.append(False)
                                    out=analyze_rss_positions(genes)
                                    if type=="-n":
                                        for n in nonamers:
                                            if n[0]==out[0]:
                                                if out[2]==True:
                                                    n[2]+=1
                                                if out[3]==True:
                                                    n[1]+=1
                                            if n[0]==out[4]:
                                                if out[6]==True:
                                                    n[4]+=1
                                                if out[7]==True:
                                                    n[3]+=1
                                        hap_found=[]
                                        for g in genes:
                                            for h in nonamers:
                                                if g[7]==h[0]:
                                                    if h[0] not in hap_found:
                                                        h[7]+=1
                                                        hap_found.append(h[0])
                                                    if h[0] not in species_rss:
                                                        h[6]+=1
                                                        species_rss.append(h[0])
                                                    h[5]+=1
                                    if type=="-h":
                                        for n in heptamers:
                                            if n[0]==out[0]:
                                                if out[2]==True:
                                                    n[2]+=1
                                                if out[3]==True:
                                                    n[1]+=1
                                            if n[0]==out[4]:
                                                if out[6]==True:
                                                    n[4]+=1
                                                if out[7]==True:
                                                    n[3]+=1
                                        hap_found=[]
                                        for g in genes:
                                            for h in heptamers:
                                                if g[7]==h[0]:
                                                    if h[0] not in hap_found:
                                                        h[7]+=1
                                                        hap_found.append(h[0])
                                                    if h[0] not in species_rss:
                                                        h[6]+=1
                                                        species_rss.append(h[0])
                                                    h[5]+=1

filter_threshold=10 #threshold for filtering out heptamers with low number of haplotypes
if type=="-h":
    with open(os.path.dirname(input_file)+"/"+locus+"_heptamer_analysis.csv","w",newline="") as write:
        writer=csv.writer(write)
        writer.writerow(["Heptamer","Number in the First 25% of Genes","Number of First Genes","Number in the Last 25% of Genes","Number of Last Genes","Number of Occurences","Number of Species","Number of Haplotypes"])
        for hep in heptamers:
            if int(hep[-1])>=filter_threshold:
                writer.writerow(hep)
        write.close()

    subprocess.run(["python","meme_figure_maker.py",os.path.dirname(input_file)+"/"+locus+"_heptamer_analysis.csv"])
    subprocess.run(["python","horizontal_bar_chart.py",os.path.dirname(input_file)+"/"+locus+"_heptamer_analysis.csv"])

elif type=="-n":
    with open(os.path.dirname(input_file)+"/"+locus+"_nonamer_analysis.csv","w",newline="") as write:
        writer=csv.writer(write)
        writer.writerow(["Nonamer","Number in the First 25% of Genes","Number of First Genes","Number in the Last 25% of Genes","Number of Last Genes","Number of Occurences","Number of Species","Number of Haplotypes"])
        for hep in nonamers:
            if int(hep[-1])>=filter_threshold:
                writer.writerow(hep)
        write.close()

    subprocess.run(["python","meme_figure_maker.py",os.path.dirname(input_file)+"/"+locus+"_nonamer_analysis.csv"])
    subprocess.run(["python","horizontal_bar_chart.py",os.path.dirname(input_file)+"/"+locus+"_nonamer_analysis.csv"])
    