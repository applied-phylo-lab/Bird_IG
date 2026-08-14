import os
import sys
import csv

contig_list=[]
write_list=[]

with open("contig_data/good_contigs.csv","r") as good:
    reader=csv.reader(good)
    header=next(reader)
    for row in reader:  
        contig_list.append([row[3].split(",")[0].replace("'","").replace("[","").replace(" ",""),row[0],row[1],row[3].split(",")[1].replace("'","").replace("]","").replace(" ","")])
        write_list.append(row)
    good.close()

with open("contig_data/multiple_contigs.csv","r") as multiple:
    reader=csv.reader(multiple)
    header=next(reader)
    for row in reader:
        cont=row[3]
        contigs=[]
        for c in cont.split("],"):
            contigs.append([c.split(",")[0].replace('[','').replace(']','').replace("'","").replace(" ",""),c.split(",")[1].replace('[','').replace(']','').replace("'",""),row[0]])
        temp_contig_list=[]
        for contig in contigs:
            found=False
            for q in temp_contig_list:
                if contig[0]==q[0]:
                    found=True
            if found==False:
                temp_contig_list.append([contig[0],int(contig[1])])
            else:
                for cont1 in temp_contig_list:
                    if contig[0]==cont1[0]:
                        cont1[1]=int(cont1[1])+int(contig[1])
        for contig in temp_contig_list:
            contig_list.append([contig[0],row[0],row[1],contig[1]])
        write_list.append([row[0],row[1],len(temp_contig_list),str(temp_contig_list).replace('[[','[').replace(']]',']')])

with open("contig_data/missing_contigs.csv","r") as missing:
    reader=csv.reader(missing)
    header=next(reader)
    for row in reader:
        write_list.append(row)
    missing.close()

with open("contig_data/finalized_contig_data.csv","w") as write_file:
    writer=csv.writer(write_file)
    writer.writerow(["Source","Locus","Number of Contigs","[Contig Name, Number of Genes]"])
    for w in write_list:
        writer.writerow(w)
    write_file.close()
with open("contig_data/ig_contig_list.csv","w") as write_file:
    writer=csv.writer(write_file)
    writer.writerow(["Contig","Source","Locus","Number of Genes (before filtering)"])
    for w in contig_list:
        writer.writerow(w)
    write_file.close()