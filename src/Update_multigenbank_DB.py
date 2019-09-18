#!/usr/bin/env python3
#Plasmid Script 2018/05/30 (LUDO: Ludovic V. Mallet)
#./complement_multigenbank.py -i 20180613_Plasmid_DB/20180611_plasmidDB_sequence.gb -d 20180613_Plasmid_DB/20180613_sequence22425.seq


import Bio
import time
import copy
import csv
import ete3
import argparse
from Bio import SeqIO
from Bio import Entrez
import random

Entrez.email = "pierre-emmanuel.douarre@anses.fr"

#integration of the command line argument
# -i & -t correspond to the input files and -s is an option for the the maximum size
parser = argparse.ArgumentParser()
parser.add_argument("-i", action="store", dest="filin", required=True, help="The input multigenbank file")
parser.add_argument("-d", action="store", dest="ids", required=True, help="the list of ids file")
params = parser.parse_args()
#required=False means the argument is optional


idss={}
Parsed_entry=0
rec_genbank = SeqIO.parse(params.filin, "gb")
for rec_gb in rec_genbank:
        Parsed_entry=Parsed_entry+1
        idss[str(rec_gb.id)] = str(rec_gb.id)


filids=open(params.ids,"r")    #open the tsv in a read mode
expected_ids={}

for i in filids.readlines():
   iid=i.strip()
   expected_ids[str(iid)]= str(iid)
filids.close()

rest = list(i if i not in idss.keys() else 'NA' for i in expected_ids.keys())
templist= [item for item in rest if item != 'NA']

print(" ".join(templist))


chunck_size=100
ids=templist

for i in range(0, len(ids), chunck_size):
	chunk=ids[i:i + chunck_size]
	print(chunk)
	print("\n\n\n\n current #"+str(i)+" out of "+str(len(ids))+"\n")
	chunk_gb_file="chunk"+str(i)+".gb"
	
        #if os.path.isfile(last_gb_file_name_in_chunck) and os.stat(last_gb_file_name_in_chunck).st_size != 0:
                #rec_chunck_genbank= (SeqIO.read(str(y)+".gb", "genbank") for y in chunk)
                ##rec_chunck_xml = (Entrez.parse(str(z)+".xml") for z in chunk)
                ##rec_gb = SeqIO.read(gb_file_name, "genbank")
                ##xml_file_handler=open(xml_file_name,"r")
                ##rec_xml = Entrez.read(xml_file_handler)
        #else:
	filin_gb = Entrez.efetch(db="nucleotide", id=chunk, rettype="gbwithparts", retmode="text")
	rec_chunck_genbank = SeqIO.parse(filin_gb, "genbank")
	filout=open(chunk_gb_file,"w")  #take the rec.id for each genbank and add .gbk
	for j in rec_chunck_genbank:
		SeqIO.write(j, filout,"genbank")
	#filout.write(filin_gb)
	time.sleep(random.uniform(1,6))
	filout.close()


