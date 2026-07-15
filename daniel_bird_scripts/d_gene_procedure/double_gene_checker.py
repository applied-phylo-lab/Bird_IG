import os
import sys
import subprocess
import csv

genes=[]
with open("bird_d_genes.csv","r") as read:
    reader=csv.reader(read)
    header=next(reader)
    for row in reader:
        genes.append(row) 
    read.close()

haplotypes=[]
for g in genes:
    if haplotypes!=[]:
        found=False
        for h in haplotypes:
            if g[0]==h[0]:
                found=True
                h[1].append(g)
        if found==False:
            haplotypes.append([g[0],[g]])
    else:
        haplotypes.append([g[0],[g]])

species=[]
for g in genes:
    if species!=[]:
        found=False
        for h in species:
            if g[0].split("/")[0]+"/"+g[0].split("/")[1] == h[0]:
                found=True
                h[1].append(g)
        if found==False:
            species.append([g[0].split("/")[0]+"/"+g[0].split("/")[1],[g]])
    else:
        species.append([g[0].split("/")[0]+"/"+g[0].split("/")[1],[g]])

orders=[]
for g in genes:
    if orders!=[]:
        found=False
        for h in orders:
            if g[0].split("/")[0] == h[0]:
                found=True
                h[1].append(g)
        if found==False:
            orders.append([g[0].split("/")[0],[g]])
    else:
        orders.append([g[0].split("/")[0],[g]])

for h in haplotypes:
    pairs=0
    n=0
    for h1 in h[1]:
        if n+1!=len(h[1]) and h[1][n][4]!=h[1][n+1][4]:
            if h[1][n][4]=="+":
                current=int(h[1][n][3])
                next_gene=int(h[1][n+1][3])-len(h[1][n+1][5])
            elif h[1][n][4]=="-":
                current=int(h[1][n][3])-len(h[1][n][5])
                next_gene=int(h[1][n+1][3])
            
            if current==next_gene:
                pairs+=1
        n+=1
    h.append(pairs)
    #h.append(round((h[2]*2)/len(h[1])*100,2))
    print("\nHaplotype: ",h[0], "\nNumber of Genes: ", len(h[1]), "\nNumber of Genes with an Opposite Twin: ", h[2]*2, "Percentage of Total: ", round((h[2]*2)/len(h[1])*100,2))

for s in species:
    pairs=0
    for h in haplotypes:
        if s[0] in h[0]:
            pairs=pairs+h[2]
    s.append(pairs)
    #s.append(round((s[2]*2)/len(s[1])*100,2))
    print("\nSpecies: ",s[0], "\nNumber of Genes: ", len(s[1]), "\nNumber of Genes with an Opposite Twin: ", s[2]*2, "Percentage of Total: ", round((s[2]*2)/len(s[1])*100,2))

for o in orders:
    pairs=0
    for h in haplotypes:
        if o[0] in h[0]:
            pairs=pairs+h[2]
    o.append(pairs)
    #o.append(round((o[2]*2)/len(o[1])*100,2))
    print("\nOrder: ", o[0], "\nNumber of Genes: ", len(o[1]), "\nNumber of Genes with an Opposite Twin: ", o[2]*2, "Percentage of Total: ", round((o[2]*2)/len(o[1])*100,2))

try:
    hap_data=[]
    with open("haplotypes_d_genes.csv","r") as read:
        reader=csv.reader(read)
        header=next(reader)
        for row in reader:
            hap_data.append(row)
        read.close()

    for hap in hap_data:
        found=False
        for h in haplotypes:
            if hap[2] in h[0]:
                found=True
                hap.append(int(h[2])*2)
        if found==False:
            hap.append(0)

    with open("haplotypes_d_genes.csv","w",newline="") as write:
        writer=csv.writer(write)
        writer.writerow(["Order","Species","Haplotype","Number of V Genes","Number of D Genes","Number +","Number -","Number with an Opposite Twin"])
        for row in hap_data:
            writer.writerow(row)
        write.close()
except:
    print("Error saving results to haplotypes_d_genes.csv")


try:
    hap_data=[]
    with open("species_d_genes.csv","r") as read:
        reader=csv.reader(read)
        header=next(reader)
        for row in reader:
            hap_data.append(row)
        read.close()

    for hap in hap_data:
        found=False
        for h in species:
            if hap[1] in h[0]:
                found=True
                hap.append(int(h[2])*2)
        if found==False:
            hap.append(0)

    with open("species_d_genes.csv","w",newline="") as write:
        writer=csv.writer(write)
        writer.writerow(["Order","Species","Number of V Genes","Number of D Genes","Number +","Number -","Number with an Opposite Twin"])
        for row in hap_data:
            writer.writerow(row)
        write.close()
except:
    print("Error saving results to species_d_genes.csv")


try:
    hap_data=[]
    with open("order_d_genes.csv","r") as read:
        reader=csv.reader(read)
        header=next(reader)
        for row in reader:
            hap_data.append(row)
        read.close()

    for hap in hap_data:
        found=False
        for h in orders:
            if hap[0] in h[0]:
                found=True
                hap.append(int(h[2])*2)
        if found==False:
            hap.append(0)

    with open("order_d_genes.csv","w",newline="") as write:
        writer=csv.writer(write)
        writer.writerow(["Order","Number of V Genes","Number of D Genes","Number +","Number -","Number with an Opposite Twin"])
        for row in hap_data:
            writer.writerow(row)
        write.close()
except:
    print("Error saving results to order_d_genes.csv")
