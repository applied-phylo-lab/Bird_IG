import os
import sys
import subprocess
import shutil

location = sys.argv[1]
locus = sys.argv[2]
out_folder = sys.argv[3]
spacing = sys.argv[4] #-12 for 12bp spacer, -23 for 23 bp spacer
limit = sys.argv[5] #-o for orders, -s for species, -h for haplotype
if limit!="-n":
    limit_value = sys.argv[6]
else:
    limit_value=""

target_genes = "target_genes_extrator.py"
extract_rss = "rss_extrator.py"
combine_rss = "rss_result_table_creator.py"
analyze_rss = "rss_analysis.py"

if os.path.isfile("target_genes.csv"):
    os.remove("target_genes.csv")
if os.path.isfile("extracted_rss_zones_directdown.csv"):
    os.remove("extracted_rss_zones_directdown.csv")
if os.path.isfile("extracted_rss_zones_onedown.csv"):
    os.remove("extracted_rss_zones_onedown.csv")
if os.path.isfile("extracted_rss_zones_twodown.csv"):
    os.remove("extracted_rss_zones_twodown.csv")
if os.path.isfile("combined_rss_zones.csv"):
    os.remove("combined_rss_zones.csv")
if os.path.isfile("heptamer_list.csv"):
    os.remove("heptamer_list.csv")
if os.path.isfile("nonamer_list.csv"):
    os.remove("nonamer_list.csv")

subprocess.run(["python",target_genes,location,locus,limit,limit_value])
subprocess.run(["python",extract_rss,"-dd",locus,spacing])
subprocess.run(["python",extract_rss,"-od",locus,spacing])
subprocess.run(["python",extract_rss,"-td",locus,spacing])
subprocess.run(["python",combine_rss,"extracted_rss_zones_directdown.csv","extracted_rss_zones_onedown.csv","extracted_rss_zones_twodown.csv"])
if os.path.isdir(out_folder)==False:
    os.mkdir(out_folder)
else:
    shutil.rmtree(out_folder)
    os.mkdir(out_folder)
shutil.move("target_genes.csv",out_folder)
shutil.move("extracted_rss_zones_directdown.csv",out_folder)
shutil.move("extracted_rss_zones_onedown.csv",out_folder)
shutil.move("extracted_rss_zones_twodown.csv",out_folder)
shutil.move("combined_rss_zones.csv",out_folder)
shutil.move("heptamer_list.csv",out_folder)
shutil.move("nonamer_list.csv",out_folder)

if os.path.isfile(locus+"_gene_fasta.fasta"):
    os.remove(locus+"_gene_fasta.fasta")
if os.path.isfile("samtools_out.fasta"):
    os.remove("samtools_out.fasta")

subprocess.run(["python","rss_analysis.py",locus,out_folder+"/combined_rss_zones.csv","-h"])
subprocess.run(["python","rss_analysis.py",locus,out_folder+"/combined_rss_zones.csv","-n"])