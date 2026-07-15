import os
import sys
import time
import csv
csv.field_size_limit(100000000000)
import subprocess
import itertools
import shutil

reference=sys.argv[1]
#-file  :  default (2 hep, 6 non) / exact (0 hep, 0 non) /  heptamer threshold (num) - nonamer threshold (num)  :  path to csv in the format [Sequence, Locus, RSS Type (Heptamer/Nonamer), Relative Location (Upstream/Downstream)]  :  locus
#-base  :  default (2 hep, 6 non) / exact (0 hep, 0 non) /  heptamer threshold (num) - nonamer threshold (num)  :  locus
if reference.startswith("-file:")==False and reference.startswith("-base:")==False:
    print("enter valid reference argument")
    sys.exit(1)
if reference.startswith("-file:"):
    locus=reference.split(":")[3]
elif reference.startswith("-base:"):
    locus=str(reference.split(":")[2])

plot=sys.argv[2]
#-plot (makes locus plot and saves genes to clean birds) / nplot (doesn't make locus plot or save genes to clean birds)  :  ysum (run summary stats script on results) / nsum (doesn't run summary stats script on results)
if plot.endswith("nsum")==False and plot.startswith("-nplot")==False and plot.endswith("ysum")==False and plot.startswith("-plot")==False:
    print("enter valid plotting argument")
    sys.exit(1)

#-both_same  :  heptamer threshold (num) - nonamer threshold (num)  |  -none (no up and down heptamers and nonamers have to be equal)
both_same=sys.argv[3]
if both_same.startswith("-both_same:")==False and both_same!="-none":
    print("enter valid both_same argument")
    sys.exit(1)

# -motif : heptamer threshold (num) - nonamer threshold (num) : prop/top  |  -none
motif_filter=sys.argv[4]
if motif_filter.startswith("-motif:")==False and motif_filter!="-none":
    print("enter valid motif_filter argument")
    sys.exit(1)

# -size : bp threshold (num)  |  -none
size_filter=sys.argv[5]
if size_filter.startswith("-size:")==False and size_filter!="-none":
    print('enter valid size filtering argument')
    sys.exit(1)


# -outdel / outkeep : folder name |  -none
save_folder=sys.argv[6]
if save_folder.startswith("-outdel:")==False and save_folder.startswith("-outkeep:")==False and save_folder!="-none":
    print('enter valid save folder argument')
    sys.exit(1)
if save_folder.startswith("-outdel:"):
    out_folder=save_folder.replace("-outdel:","")
    save_type="w"
elif save_folder.startswith("-outkeep:"):
    out_folder=save_folder.replace("-outkeep:","")
    save_type="a"
else:
    out_folder=""
if os.path.isdir(out_folder)==False and out_folder!="":
    os.mkdir(out_folder)

# -f : Text in the relative path of the target (ex: "Landfowl" to get all landfowl, "IMGT_Chicken" to get IMGT_Chicken, ect)  |  -n
filter=sys.argv[7]
if filter!="-n" and filter.startswith("-f:")==False:
    print("enter valid data filtering argument")
    sys.exit(1)

params=[reference,plot,both_same,motif_filter,size_filter,filter]
with open(out_folder+"/params.txt",save_type,newline="") as write:
    writer=csv.writer(write)
    if save_type=="a":
        writer.writerow([ ])
    for p in params:
        writer.writerow([p])
    write.close()

contigs=[]
with open("/local/storage/kav67/clean_birds/ig_contig_list.csv","r") as contig_read:
    reader=csv.reader(contig_read)
    header=next(reader)
    for row in reader:
        contigs.append(row[:3])
    contig_read.close()

genomes=[]
with open("/local/storage/kav67/clean_birds/bird_genome_paths.csv","r") as genome_read:
    reader=csv.reader(genome_read)
    header=next(reader)
    for row in reader:
        genomes.append(row)
    genome_read.close()

genes=[]
with open("/local/storage/kav67/clean_birds/gene_list.csv","r") as gene_read:
    reader=csv.reader(gene_read)
    header=next(reader)
    for row in reader:
        genes.append(row)
    gene_read.close()

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

def downstream_extract(contig,genes,genome,out_folder):
    gene_pos=[]
    for gene in genes:
        gene_pos.append(int(gene[3]))
    first_gene=min(gene_pos)
    last_gene=max(gene_pos)
    for g in genes:
        if int(last_gene)==int(g[3]):
            last_gene=int(last_gene)+len(g[5])
    subprocess.run(["samtools","faidx",genome])
    with open(genome+".fai","r") as samtools_read:
        reader=csv.reader(samtools_read,delimiter="\t")
        for row in reader:
            if row[0]==contig:
                contig_size=row[1]
        samtools_read.close()
    upstream_section=contig+":1-"+str(first_gene)
    v_cluster=contig+":"+str(first_gene)+"-"+str(last_gene)
    downstream_section=contig+":"+str(last_gene)+"-"+str(contig_size)
    subprocess.run(["samtools","faidx",genome,upstream_section,"-n","0","-o",out_folder+"/upstream.fasta"])
    subprocess.run(["samtools","faidx",genome,v_cluster,"-n","0","-o",out_folder+"/v_cluster.fasta"])
    subprocess.run(["samtools","faidx",genome,downstream_section,"-n","0","-o",out_folder+"/downstream.fasta"])
    
    sections=[]
    for f in [out_folder+"/upstream.fasta",out_folder+"/v_cluster.fasta",out_folder+"/downstream.fasta"]:
        with open(f,"r") as read:
            reader=csv.reader(read)
            header=next(reader)
            for row in reader:
                section=str(row[0])
            read.close()
        reverse_section=reverse_complement(section)
        if "upstream.fasta" in f:
            start_num=0
        elif "v_cluster.fasta" in f:
            start_num=(first_gene-1)
        elif "downstream.fasta" in f:
            start_num=(last_gene-1)
        sections.append([f.replace(out_folder,"").replace(".fasta",""),section,reverse_section,start_num])
    return sections

def d_gene_extract(section,downstream_heptamers,upstream_heptamers,downstream_nonamers,upstream_nonamers,thresholds,section_start_num):
    global both_same
    global locus
    n=0
    if locus=="IGH":
        min_d_size=20
        max_d_size=71
        up_spacer=12
        down_spacer=12
    elif locus=="TRB":
        min_d_size=10
        max_d_size=40
        up_spacer=23
        down_spacer=12
    elif locus=="TRD":
        min_d_size=10
        max_d_size=40
        up_spacer=23
        down_spacer=12

    d_genes=[]
    for s in section:
        upstream_nonamer_zone=section[n:(n+9)]
        gene=[]
        if thresholds==[]:
            if upstream_nonamer_zone in upstream_nonamers:
                upstream_heptamer_zone=section[(n+9+up_spacer):(n+9+up_spacer+7)]
                if upstream_heptamer_zone in upstream_heptamers:
                    #if there are multiple potential downstream rss matches from an upstream rss using different d sizes, this will take the match from the largest potential d_gene
                    for d_size in range(min_d_size, max_d_size):
                        d_gene=section[(n+9+up_spacer+7):((n+9+up_spacer+7)+d_size)]
                        downstream_heptamer_zone=section[((n+9+up_spacer+7)+d_size):(((n+9+up_spacer+7)+d_size)+7)]
                        if downstream_heptamer_zone in downstream_heptamers:
                            downstream_nonamer_zone=section[(((n+9+up_spacer+7)+d_size)+(7+down_spacer)):(((n+9+up_spacer+7)+d_size)+(7+down_spacer)+9)]
                            if downstream_nonamer_zone in downstream_nonamers:
                                if both_same!="-none":
                                    both_same_hep=both_same.split(":")[1].split("-")[0]
                                    both_same_non=both_same.split(":")[1].split("-")[1]
                                    if heptamer_point_mutation([upstream_heptamer_zone],downstream_heptamer_zone,both_same_hep)==False:
                                        continue
                                    if nonamer_point_mutation([upstream_nonamer_zone],downstream_nonamer_zone,both_same_non)==False:
                                        continue
                                whole_gene=section[n:(((n+9+up_spacer+7)+d_size)+(7+down_spacer)+9)]
                                if section_start_num<0:
                                    gene_pos=abs(section_start_num+(n+9+up_spacer+7)-1)
                                else:
                                    gene_pos=abs(section_start_num+(n+9+up_spacer+7)+1)
                                gene=[gene_pos,upstream_heptamer_zone,upstream_nonamer_zone,d_gene,downstream_heptamer_zone,downstream_nonamer_zone,whole_gene]
        else:
            if nonamer_point_mutation(upstream_nonamers,upstream_nonamer_zone,thresholds[1])==True:
                upstream_heptamer_zone=section[(n+9+up_spacer):(n+9+up_spacer+7)]
                if heptamer_point_mutation(upstream_heptamers,upstream_heptamer_zone,thresholds[0])==True:
                    #if there are multiple potential downstream rss matches from an upstream rss using different d sizes, this will take the match from the largest potential d_gene
                    for d_size in range(min_d_size, max_d_size):
                        d_gene=section[(n+9+up_spacer+7):((n+9+up_spacer+7)+d_size)]
                        downstream_heptamer_zone=section[((n+9+up_spacer+7)+d_size):(((n+9+up_spacer+7)+d_size)+7)]
                        if heptamer_point_mutation(downstream_heptamers,downstream_heptamer_zone,thresholds[0])==True:
                            downstream_nonamer_zone=section[(((n+9+up_spacer+7)+d_size)+(7+down_spacer)):(((n+9+up_spacer+7)+d_size)+(7+down_spacer)+9)]
                            if nonamer_point_mutation(downstream_nonamers,downstream_nonamer_zone,thresholds[1])==True:
                                if both_same!="-none":
                                    both_same_hep=int(both_same.split(":")[1].split("-")[0])
                                    both_same_non=int(both_same.split(":")[1].split("-")[1])
                                    if heptamer_point_mutation([upstream_heptamer_zone],downstream_heptamer_zone,both_same_hep)==False:
                                        continue
                                    if nonamer_point_mutation([upstream_nonamer_zone],downstream_nonamer_zone,both_same_non)==False:
                                        continue
                                whole_gene=section[n:(((n+9+up_spacer+7)+d_size)+(7+down_spacer)+9)]
                                if section_start_num<0:
                                    gene_pos=abs(section_start_num+(n+9+up_spacer+7)-1)
                                else:
                                    gene_pos=abs(section_start_num+(n+9+up_spacer+7)+1)
                                gene=[gene_pos,upstream_heptamer_zone,upstream_nonamer_zone,d_gene,downstream_heptamer_zone,downstream_nonamer_zone,whole_gene,n]
        if gene!=[]:
            d_genes.append(gene)
            print("\n\n")
            for g in gene:
                print(g)
            print("\n\n")
        n+=1
    return d_genes

if out_folder!="":
    if save_type=="w":
        if os.path.isfile(out_folder+"/bird_d_genes.csv"):
            os.remove(out_folder+"/bird_d_genes.csv")
        with open(out_folder+"/bird_d_genes.csv","w",newline="") as write:
            writer=csv.writer(write)
            header=["Source","GeneType","Contig","Pos","Strand","Sequence","Productive","Locus","Upstream Heptamer","Upstream Nonamer","Downstream Heptamer","Downstream Nonamer","Whole Section","Location Relative to V-Cluster"]
            writer.writerow(header)
            write.close()

spec_list=[]
for contig in contigs:
    if filter.startswith("-f:"):
        found=False
        for fil in filter.replace("-f:","").split("-"):
            if fil in contig[1]:
                found=True
        if found==False:
            continue
    
    if str(contig[2])!=locus:
        continue
    
    if contig[1] not in spec_list and plot.startswith("-plot"):
        spec_list.append(contig[1])
        with open("/local/storage/kav67/clean_birds/"+contig[1]+"/"+locus+"D.csv","w",newline="") as write:
            writer=csv.writer(write)
            header=["Source","GeneType","Contig","Pos","Strand","Sequence","Productive","Locus","Upstream Heptamer","Upstream Nonamer","Downstream Heptamer","Downstream Nonamer","Whole Section","Location Relative to V-Cluster"]
            writer.writerow(header)
            write.close()

    target_genes=[]
    for gene in genes:
        if gene[0]==contig[1] and gene[2]==contig[0]:
            target_genes.append(gene)
    for g in genomes:
        if g[0]==contig[1]:
            genome=g[1]
    sections=downstream_extract(contig[0],target_genes,genome,out_folder)    
    os.remove(out_folder+"/downstream.fasta")
    os.remove(out_folder+"/upstream.fasta")
    os.remove(out_folder+"/v_cluster.fasta")

    #   sections == [[label(upstream),pos,neg],[label(downstream),pos,neg],[label(v_cluster),pos,neg]]
    print("\n",contig[1],contig[0])
    
    within_contig_heptamers=[]
    within_contig_nonamers=[]
    for t in target_genes:
        if t[8]!="" and t[8] not in within_contig_heptamers:
            within_contig_heptamers.append(t[8])
        if t[9]!="" and t[9] not in within_contig_nonamers:
            within_contig_nonamers.append(t[9])

    if reference.startswith("-file:"):
        if reference.split(":")[1]=="default":
            thresholds=[2,6]
        elif reference.split(":")[1]=="exact":
            thresholds=[]
        elif len(reference.split(":")[1].split("-"))==2:
            thresholds=[int(reference.split(":")[1].split("-")[0]),int(reference.split(":")[1].split("-")[1])]
        
        heptamers=[]
        upstream_heptamers=[]
        nonamers=[]
        upstream_nonamers=[]

        with open(reference.split(":")[2],"r") as rss_read:
            reader=csv.reader(rss_read)
            header=next(reader)
            if row[1]==locus:
                if row[2]=="Heptamer":
                    if row[3]=="Upstream":
                        upstream_heptamers.append(row[0])
                    if row[3]=="Downstream":
                        heptamers.append(row[0])
                if row[2]=="Nonamer":
                    if row[3]=="Upstream":
                        upstream_nonamers.append(row[0])
                    if row[3]=="Downstream":
                        nonamers.append(row[0])
    
    elif reference.startswith("-base:"):
        if reference.split(":")[1]=="default":
            thresholds=[2,6]
        elif reference.split(":")[1]=="exact":
            thresholds=[]
        elif len(reference.split(":")[1].split("-"))==2:
            thresholds=[int(reference.split(":")[1].split("-")[0]),int(reference.split(":")[1].split("-")[1])]
        
        if locus=="IGH":
            heptamers=['CACGGTG']
            upstream_heptamers=["CACTGTG","CACCGTG"]
            nonamers=["ACAAAAACC"]
            upstream_nonamers=["GGATTTTGG"]
        
        elif locus=="TRB":
            heptamers=["TGTTTTT"]
            upstream_heptamers=["CACAATG"]
            nonamers=["ACCGTTGTG"]
            upstream_nonamers=["CAAAAACCT"]
        
        elif locus=="TRD":
            heptamers=["ATTTTTT","TGTTTTT"]
            upstream_heptamers=["CACAATG","CACTGCG"]
            nonamers=["ATCACTGTG","ACTACTGTG"]
            upstream_nonamers=["ACAAAAACC","ACACAAACG"]

    print(upstream_heptamers)
    print(heptamers)
    print("\n")
    print(upstream_nonamers)
    print(nonamers)
    
    write_list=[]
    for s in sections:
        print(s[0]," positive: ",s[1][:50]," negative: ",s[2][:50])
        positive_d_genes=d_gene_extract(s[1],heptamers,upstream_heptamers,nonamers,upstream_nonamers,thresholds,int(s[3])) #positive
        negative_d_genes=d_gene_extract(s[2],heptamers,upstream_heptamers,nonamers,upstream_nonamers,thresholds,-(int(s[3])+int(len(s[2])))) #negative

        for p in positive_d_genes:
            wr=[contig[1],"D",contig[0],int(p[0]),"+",p[3],"",contig[2],p[1],p[2],p[4],p[5],p[6],s[0]]
            write_list.append(wr)
        for n in negative_d_genes:
            wr=[contig[1],"D",contig[0],int(n[0]),"-",n[3],"",contig[2],n[1],n[2],n[4],n[5],n[6],s[0]]
            write_list.append(wr)
    
    write_list.sort(key=lambda x: x[3]) 
    if motif_filter.startswith("-motif:") and len(write_list)>5:
        up_heptamers=[]
        down_heptamers=[]
        up_nonamers=[]
        down_nonamers=[]
        for w in write_list:
            up_heptamers.append(w[8])
            up_nonamers.append(w[9])
            down_heptamers.append(w[10])
            down_nonamers.append(w[11])
        
        def motif_create(seqs,motif_filter):
            filter_type=motif_filter.split(":")[2]
            num=len(seqs[0])
            motif=[]
            for n in range(num):
                motif.append([["A",0],["C",0],["T",0],["G",0]])
            for s in seqs:
                v=0
                for n in s:
                    if n.capitalize()=="A":
                        b=0
                    elif n.capitalize()=="C":
                        b=1
                    elif n.capitalize()=="T":
                        b=2
                    elif n.capitalize()=="G":
                        b=3
                    motif[v][b][1]+=1
                    v+=1
            if filter_type=="top":
                top=[]
                for m in motif:
                    tot=int(m[0][1])+int(m[1][1])+int(m[2][1])+int(m[3][1])
                    high=["","",0]
                    for m1 in m:
                        m1.append(round((int(m1[1])/tot)*100,2))
                        if m1[2]>high[2]:
                            high=m1
                    top.append([high[0]])
            elif filter_type=="prop":
                top=[]
                for m in motif:
                    tot=int(m[0][1])+int(m[1][1])+int(m[2][1])+int(m[3][1])
                    high=["","",0]
                    for m1 in m:
                        m1.append(round((int(m1[1])/tot)*100,2))
                    m.sort(key=lambda x: x[2], reverse=True) 

                    if int(m[1][2])==0 or int(m[0][2])/int(m[1][2]) >2:
                        top.append([m[0][0]])
                    else:
                        if int(m[2][2])==0 or int(m[1][2])/int(m[2][2]) > 2:
                            top.append([m[0][0],m[1][0]])
                        else:
                            if int(m[3][2])==0 or int(m[2][2])/int(m[3][2]) >2:
                                top.append([m[0][0],m[1][0],m[2][0]])
                            else:
                                top.append([m[0][0],m[1][0],m[2][0],m[3][0]])
            print(motif)
            print(top)
            def get_string_combinations(nested_list):
                combinations = itertools.product(*nested_list)
                return ["".join(comb) for comb in combinations]
            seq_list=get_string_combinations(top)
            def motif_score(top):
                score=0
                for t in top:
                    score = score + len(t)*0.25
                if len(top)==7:
                    if score<2.75:
                        thresh=3
                    elif 2.75<=score<3.5:
                        thresh=2
                    elif score>=3.5:
                        thresh=1
                elif len(top)==9:
                    if score<2.75:
                        thresh=5
                    elif 2.75<=score<3.25:
                        thresh=4
                    elif 3.25<=score<4:
                        thresh=3
                    elif score>=4:
                        thresh=2
                return thresh
            
            if motif_filter.split(":")[1]=="scoring":
                thresh=motif_score(top)
            else:
                if len(top)==7:
                    thresh=int(motif_filter.split(":")[1].split("-")[0])
                elif len(top)==9:
                    thresh=int(motif_filter.split(":")[1].split("-")[1])   
            print(seq_list)
            return [seq_list,motif,thresh]

        up_hep_motif=motif_create(up_heptamers,motif_filter)
        up_non_motif=motif_create(up_nonamers,motif_filter)
        down_hep_motif=motif_create(down_heptamers,motif_filter)
        down_non_motif=motif_create(down_nonamers,motif_filter)  
        
        up_heptamer_thresh=int(up_hep_motif[2])
        up_nonamer_thresh=int(up_non_motif[2])
        down_heptamer_thresh=int(down_hep_motif[2])
        down_nonamer_thresh=int(down_non_motif[2])

        write_list_real=[]
        for w in write_list:
            if heptamer_point_mutation(up_hep_motif[0],w[8],up_heptamer_thresh)==True and nonamer_point_mutation(up_non_motif[0],w[9],up_nonamer_thresh)==True and heptamer_point_mutation(down_hep_motif[0],w[10],down_heptamer_thresh)==True and nonamer_point_mutation(down_non_motif[0],w[11],down_nonamer_thresh)==True:
                write_list_real.append(w)
        
        write_list=write_list_real
    
    if size_filter.startswith("-size:") and len(write_list)>10:
        size_thresh=0.55
        
        gene_size_genes=[]
        for r in range(51):
            gene_size_genes.append([])
        for w in write_list:
            length=len(w[5])
            gene_size_genes[length-20].append(w[5])
        
        gene_sizes=[]
        for g in gene_size_genes:
            gene_sizes.append(len(g))
        print(gene_sizes)
        
        size_strictness=int(size_filter.split(":")[1])
        n=20
        gene_sizes_index=[]
        for g in gene_sizes:
            gls=[n,g]
            gene_sizes_index.append(gls)
            n+=1
        gene_sizes_index.sort(key=lambda x: x[1],reverse=True)
        n1=0
        keep_sizes=[]
        for g in gene_sizes_index:
            print(g)
            if n1==0:
                keep_sizes.append(gene_sizes_index[n1])
            try:
                if int(gene_sizes_index[n1+1][1])!=0 and gene_sizes_index[n1+1][1]/gene_sizes_index[n1][1] >= size_thresh:
                    keep_sizes.append(gene_sizes_index[n1+1])
                else:
                    break
            except IndexError:
                keep_sizes
            n1+=1
        print(keep_sizes)
        write_list_real=[]
        for w in write_list:
            for k in keep_sizes:
                if abs(len(w[5])-k[0])<=size_strictness:
                    if w not in write_list_real:
                        write_list_real.append(w)
        write_list=write_list_real
        print(len(write_list))

    '''for w in write_list:
        if w[4]=="-":
            w[3]=int(w[3])-len(w[5])'''

    with open(out_folder+"/bird_d_genes.csv","a",newline="") as write:
        writer=csv.writer(write)
        for w in write_list:
            writer.writerow(w) 
        write.close()
    
    if plot.startswith("-plot"):
        subprocess.run(["python","ig_contig_plotter.py","-d","-n","-contig:"+str(contig[0]),"-h",contig[1].split("/")[-1]],cwd="../")
        with open("/local/storage/kav67/clean_birds/"+contig[1]+"/"+locus+"D.csv","a",newline="") as write:
            writer=csv.writer(write)
            for w in write_list:
                writer.writerow(w)
            write.close()

if plot.endswith("ysum"):
    shutil.copy("summary_stats.py",out_folder)
    shutil.copy("double_gene_checker.py",out_folder)
    time.sleep(1)
    if filter=="-n":
        filter="-none"
    subprocess.run(["python", "summary_stats.py","-thresh:1-1:occur",filter],cwd=out_folder)
    os.remove(out_folder+"/summary_stats.py")
    os.remove(out_folder+"/double_gene_checker.py")