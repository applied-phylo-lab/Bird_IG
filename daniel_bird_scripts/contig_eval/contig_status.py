import os
import sys
import csv
import copy

order_list=["Cormorants","Cranes","Doves","Eagles","Falcons","Hummingbirds","Ibises","Landfowl","MiscBirds","Owls","Parrots","Plovers","Songbirds","Suboscines","Waterfowl","Woodpeckers"] #names of all bird order folders

folder=sys.argv[1]

write_list=[]
good4=0
none4=0
more4=0
for f in os.listdir(folder):
    if f in order_list and os.path.isdir(folder+"/"+f):
        good3=0
        more3=0
        none3=0
        order=[]
        for f1 in os.listdir(folder+"/"+f):
            if os.path.isdir(folder+"/"+f+"/"+f1):
                good2=0
                more2=0
                none2=0
                species=[]
                for f2 in os.listdir(folder+"/"+f+"/"+f1):
                    if os.path.isdir(folder+"/"+f+"/"+f1+"/"+f2):
                        good1=0
                        more1=0
                        none1=0
                        for f3 in os.listdir(folder+"/"+f+"/"+f1+"/"+f2):
                            good=0
                            more=0
                            none=0
                            if f3.replace("combined_genes_","").replace(".txt","") in ["TRA","TRB","TRG","TRD"]:
                                with open(folder+"/"+f+"/"+f1+"/"+f2+"/"+f3,"r") as read:
                                    reader=csv.reader(read,delimiter="\t")
                                    header=next(reader)
                                    contigs=[]
                                    list=[]
                                    contig_list=[]
                                    for row in reader:
                                        list.append(row[1])
                                        if row[1] not in contigs:
                                            contigs.append(row[1])
                                    greatest=0
                                    for c in contigs:
                                        contig_num=list.count(c)
                                        if contig_num>greatest:
                                            greatest=contig_num
                                        contig_list.append([c,contig_num])
                                    contig_list_real=[]
                                    for c in contig_list:
                                        if ((c[1]-greatest)/greatest)*100>=-80:
                                            contig_list_real.append(c)
                                    locus=f3.replace("combined_genes_","").replace(".txt","")
                                    wri=[f+"/"+f1+"/"+f2,locus,len(contig_list_real),str(contig_list_real).replace("[[","[").replace("]]","]")]
                                    #print(wri)
                                    write_list.append(wri)
                                    species.append(wri)
                                    order.append(wri)
                                    if len(contig_list_real)==1:
                                        good+=1
                                    elif len(contig_list_real)==0:
                                        none+=1
                                    elif len(contig_list_real)>1:
                                        more+=1
                                good1=good1+good
                                none1=none1+none
                                more1=more1+more
                        #print(f2+":\n   -Good: "+str(good1)+"\n   -Multiple: "+str(more1)+"\n   -None: "+str(none1))    
                        good2=good2+good1
                        none2=none2+none1
                        more2=more2+more1   
                #print(f1+":\n   -Good: "+str(good2)+"\n   -Multiple: "+str(more2)+"\n   -None: "+str(none2))  
                '''with open(folder+"/"+f+"/"+f1+"/contig_data.csv","w",newline="") as cd:
                    writer=csv.writer(cd)
                    writer.writerow(["Source","Locus","Number of Potential Contigs","[Contig Name, Number of Genes]"])
                    for sp in species:
                        writer.writerow(sp)
                    cd.close()'''
                good3=good3+good2
                none3=none3+none2
                more3=more3+more2
        print(f+":\n   -Good: "+str(good3)+"\n   -Multiple: "+str(more3)+"\n   -None: "+str(none3))
        '''with open(folder+"/"+f+"/contig_data.csv","w",newline="") as cd:
            writer=csv.writer(cd)
            writer.writerow(["Source","Locus","Number of Potential Contigs","[Contig Name, Number of Genes]"])
            for od in order:
                writer.writerow(od)
            cd.close()'''
        good4=good4+good3
        none4=none4+none3
        more4=more4+more3
print("Birds:\n   -Good: "+str(good4)+"\n   -Multiple: "+str(more4)+"\n   -None: "+str(none4))
with open("contig_summary.csv","w",newline="") as write:
    writer=csv.writer(write)
    writer.writerow(["Source","Locus","Number of Potential Contigs","[Contig Name, Number of Genes]"])
    for r in write_list:
        writer.writerow(r)
    write.close()

if os.path.isdir("contig_data")==True:
    os.rmdir("contig_data")
os.mkdir("contig_data")    

with open("contig_data/missing_contigs.csv","w",newline="") as mch:
    writer=csv.writer(mch)
    writer.writerow(["Source","Locus","Number of Potential Contigs","[Contig Name, Number of Genes]"])
    for r in write_list:
        if str(r[2])=="0":
            writer.writerow(r)
    mch.close()
with open("contig_data/multiple_contigs.csv","w",newline="") as uch:
    writer=csv.writer(uch)
    writer.writerow(["Source","Locus","Number of Potential Contigs","[Contig Name, Number of Genes]"])
    for r in write_list:
        if int(r[2])>1:
            writer.writerow(r)
    uch.close()
with open("contig_data/good_contigs.csv","w",newline="") as gch:
    writer=csv.writer(gch)
    writer.writerow(["Source","Locus","Number of Potential Contigs","[Contig Name, Number of Genes]"])
    for r in write_list:
        if str(r[2])=="1":
            writer.writerow(r)
    gch.close()