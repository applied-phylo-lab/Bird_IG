import os
import copy
import csv
import time
import subprocess
import sys

def reverse_complement(sequence):
    complement = {
        'A': 'T', 'T': 'A',
        'G': 'C', 'C': 'G',
        'a': 't', 't': 'a',
        'g': 'c', 'c': 'g',
        'N': 'N', 'n': 'n' 
    }
    rev_comp = ''.join(complement.get(base, base) for base in reversed(sequence))
    
    return rev_comp

def heptamer_point_mutation(str1, str2, num_p):
    #checks to see if target sequence is num_p bp or more different from reference(s), returns True if this is False for any reference
    #change 'num_p' to whatever number of point mutations you want to limit
    check=[]
    for s in str1:
        if len(str2)!=len(s):
            return False
        differences = 0
        n=0
        for c in s:
            if c != str2[n]:
                differences+=1
            n+=1
        if differences <= num_p:
            check.append(True)
        else:
            check.append(False)
    if True in check:
        return True
    else:
        return False

def nonamer_point_mutation(str1, str2, num_p):
    #checks to see if target sequence is num_p bp or more different from reference(s), returns True if this is False for any reference
    #change 'num_p' to whatever number of point mutations you want to limit
    check=[]
    for s in str1:
        if len(str2)!=len(s):
            return False
        differences = 0
        n=0
        for c in s:
            if c != str2[n]:
                differences+=1
            n+=1
        if differences <= num_p:
            check.append(True)
        else:
            check.append(False)
    if True in check:
        return True
    else:
        return False
    
if sys.argv[2]=="IGH":
    #[functional imgt chicken IGH, highly expressed tufted duck IGH]
    ref_heptamers=["CACGGTG", "GATAGTG"]
    hep_num=2
    #[functional imgt chicken IGH, highly expressed tufted duck IGH]
    ref_nonamers=["CACAAAACC", "GGTGGGGTT"]
    non_num=6
    
elif sys.argv[2]=="IGL":
    #[main human IGL]
    ref_heptamers=["CACAGTG"]
    hep_num=2
    #[main human IGL]
    ref_nonamers=["ACAGAAACC"]
    non_num=6

elif sys.argv[2]=="TRA":
    #top from tcr paper (https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1359169/full)
    
    ref_heptamers=["CACGGGA","CACCCTG","CACGGGG"]
    hep_num=2

    ref_nonamers=["GCAACAACC","GCACGAACC"]
    non_num=6

elif sys.argv[2]=="TRB":
    #top from tcr paper (https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1359169/full)

    ref_heptamers=["CACTGTG"]
    hep_num=2

    ref_nonamers=["CGCAAACCT"]
    non_num=6

elif sys.argv[2]=="TRG":
    #top from tcr paper (https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1359169/full)

    ref_heptamers=["CACTATG","CACAGCG","CACAGCA","CACCATG","CACCGAG"]
    hep_num=2

    ref_nonamers=["GCAAATACT","ACAAAAAGG","ACAAAGACA","CCACAAAAC"]
    non_num=6

elif sys.argv[2]=="TRD":
    #top from tcr paper (https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2024.1359169/full)

    ref_heptamers=["CACGGGA","CACAGGA"]
    hep_num=2

    ref_nonamers=["ACAAAAATA","ACAAAAATC","ACAAAAAAC"]
    non_num=6

if sys.argv[1]=="-dd":
    title="directdown"
elif sys.argv[1]=="-od":
    title="onedown"
elif sys.argv[1]=="-td":
    title="twodown"

def rss_calc(downstream,spacer,section,delete_size,strand):
    if downstream=="-dd":
        ds_spacer=0
    elif downstream=="-od":
        ds_spacer=1
    elif downstream=="-td":
        ds_spacer=2
    
    if strand=="-":
        ds_spacer+=1


    if spacer=="-12":
        v_spacer=12
    elif spacer=="-23":
        v_spacer=23


    heptamer=section[delete_size+ds_spacer:delete_size+(ds_spacer+7)]
    spacing=section[delete_size+(ds_spacer+7):delete_size+(ds_spacer+7+v_spacer)]
    nonamer=section[delete_size+(ds_spacer+7+v_spacer):delete_size+(ds_spacer+7+v_spacer+9)]

    return [heptamer, spacing, nonamer]

spacing_term=str(sys.argv[3]) #-12 for 12bp spacer, -23 for 23 bp spacer

with open("extracted_rss_zones_"+title+".csv","w",newline="") as rss_write:
    writer=csv.writer(rss_write)
    writer.writerow(["Source","GeneType","Contig","Pos","Strand","Sequence","Productive","Locus","Extracted Section (gene +-50bp)","Heptamer","Nonamer"])
    rss_write.close()

genomes=[]
with open("bird_genome_paths.csv","r") as gp:
    reader=csv.reader(gp)
    header=next(reader)
    for row in reader:
        genomes.append(row)
    gp.close()

with open("target_genes.csv","r") as tg:
    reader=csv.reader(tg)
    header=next(reader)
    for row in reader:
        start=int(row[3])-50
        end=int(row[3])+len(row[5])+50
        for g in genomes:
            if row[0].split("/")[-1]==g[0].split("/")[-1]:
                input_fasta=g[1]
        region=str(row[2]+":"+str(start)+"-"+str(end))
        process_args=["samtools","faidx",input_fasta,region,"-o","samtools_out.fasta"]
        print(process_args)
        process=subprocess.run(process_args)
        if os.path.getsize("samtools_out.fasta")>50:
            with open("samtools_out.fasta","r") as samout:
                reader=csv.reader(samout)
                header=next(reader)
                section=""
                for row1 in reader:
                    if row1[0].startswith(">")==False:
                        section=section+str(row1[0])
                if row[4]=="-":
                    section=reverse_complement(section)
                delete_size=50+len(str(row[5]))
                rss=rss_calc(str(sys.argv[1]),spacing_term,section,delete_size,str(row[4]))
                heptamer=rss[0]
                nonamer=rss[2]
                samout.close()
            if heptamer_point_mutation(ref_heptamers,heptamer,hep_num)==True and nonamer_point_mutation(ref_nonamers,nonamer,non_num)==True:
                with open("extracted_rss_zones_"+title+".csv","a",newline="") as rss_write:
                    writer=csv.writer(rss_write)
                    write_line=copy.deepcopy(row)
                    write_line.append(section)
                    write_line.append(heptamer)
                    write_line.append(nonamer)
                    writer.writerow(write_line)
                    rss_write.close()
                os.remove("samtools_out.fasta")