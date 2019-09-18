#!/usr/bin/env python3
#Plasmid Script 2018/05/30 (LUDO: Ludovic V. Mallet)
#~/parse_multigenbank_V2.py -i Esearch_genbank_file.gb -t Accession_topology_file.tsv > ids.log


import Bio
import time
import copy
import csv
import ete3 
import argparse
from Bio import SeqIO


#integration of the command line argument
# -i & -t correspond to the input files and -s is an option for the the maximum size
parser = argparse.ArgumentParser()
parser.add_argument("-i", action="store", dest="filin", required=True, help="The input multigenbank file")
parser.add_argument("-t", action="store", dest="filetopology", required=True, help="The input topology file",default="accession_topology.tsv")
parser.add_argument("-s", action="store", dest="entry_size", required=False, help="The maximum size of sequences to keep",default=3000000)
params = parser.parse_args()
#required=False means the argument is optional 


ncbi = ete3.NCBITaxa()
desired_ranks = ['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']

def get_desired_ranks(taxid, desired_ranks):
  try :
    lineage = ncbi.get_lineage(taxid)
    #print(lineage)
    names = ncbi.get_taxid_translator(lineage)
    #print(names)
    lineage2ranks = ncbi.get_rank(names)
  	
    m=dict(zip(ncbi.get_rank(names).values(),names.values()))
    #list(names[i] if i in m else '' for i in lineage2ranks.values())
    templist=list(m[i] if i in m else 'NA' for i in desired_ranks)
   
  except:
     templist = ['NA','NA','NA','NA','NA','NA','NA']
     print(taxid)
  #templist[:] = [item for item in templist if item != '']
  return(templist)
	
# ~/parsePlasmidDB.pl Complete_Plasmid_DB.gb > accession_topology.tsv
filtsv=open(params.filetopology,"r")	#open the tsv in a read mode 
topology={}	

for i in filtsv.readlines():
   fields=i.split("\t");
   topology[fields[0]]=fields[2]
   topology[fields[1]]=fields[2]
filtsv.close()

#Create the csv and multifasta files "all" & filtered
filcsv=open("plasmid_db_headers.csv","w")	#open the csv file in write mode
filcsv.write("")	#empty the empty the csv
filcsv.close()	#close the csv		
filcsv=open("plasmid_db_headers.csv","a")	#add the headers in the csv file

fildball=open("plasmid_db.fa","w")	
fildball.write("")	
fildball.close()
fildball=open("plasmid_db.fa","a")	

fildbfiltered=open("plasmid_db_filtered.fa","w")	
fildbfiltered.write("")	
fildbfiltered.close()
fildbfiltered=open("plasmid_db_filtered.fa","a")	

filcsvfiltered=open("plasmid_db_headers_filtered.csv","w")	
filcsvfiltered.write("")	
filcsvfiltered.close()
filcsvfiltered=open("plasmid_db_headers_filtered.csv","a")	

#create dictionnaries
genbank={}
filtered_genbank={}
idss={}
topos={}
sizes={}
taxos={}
descs={}

#Parsed_entry displays the number of hits being parsed
Parsed_entry=0
rec_genbank = SeqIO.parse(params.filin, "gb")
for rec_gb in rec_genbank:
	Parsed_entry=Parsed_entry+1
	print(str(Parsed_entry))
	for i in rec_gb.features:
		if i.type=='source':
			dbxref=i.qualifiers['db_xref']
			bd_xref_tax =  filter(lambda x:'taxon' in x, dbxref)
			taxid=str(list(bd_xref_tax)[0]).strip('taxon:')
			break
	taxo= get_desired_ranks(taxid, desired_ranks)
	genbank[str(rec_gb.id)]=rec_gb
	topo=topology[str(rec_gb.id).split(".")[0]].rstrip()
	
	rec_gb2=copy.deepcopy(rec_gb)
	rec_gb2.description= "\t".join([str(rec_gb.id), str(topo), str(len(rec_gb)), str(taxo), str(rec_gb.description)])
	
	#store data
	topos[str(rec_gb.id)] = str(topo)
	idss[str(rec_gb.id)] = str(rec_gb.id)
	sizes[str(rec_gb.id)]= len(rec_gb)
	taxos[str(rec_gb.id)]= str(taxo)
	descs[str(rec_gb.id)]= str(rec_gb.description)
	
	SeqIO.write(rec_gb2, fildball,"fasta")
	filcsv.write(rec_gb2.description + "\n")
		
filcsv.close()
fildball.close()

#--------------------------------------------------------------------------------------

for i in genbank.keys():
	duplicate=0
	for k in genbank.keys():
		if idss[i]!= idss[k] and topos[i] == topos[k] and sizes[i] == sizes[k] and taxos[i]==taxos[k] and descs[i]==descs[k] :
			duplicate=1
			print("ids "+str(i)+" and "+str(k)+" have identical description, length, taxo and topology")
			newdate1 = time.strptime( genbank[i].annotations['date'], "%d-%b-%Y")
			newdate2 = time.strptime( genbank[k].annotations['date'], "%d-%b-%Y")
			if newdate1 > newdate2:
				filtered_genbank[i]=genbank[i]
			else:
				filtered_genbank[k]=genbank[k]
	if duplicate==0:
		filtered_genbank[i]=genbank[i]
		
for i in filtered_genbank.keys():		
	if sizes[i]< params.entry_size:
		rec_gb2=copy.deepcopy(genbank[i])
		rec_gb2.description= "\t".join([str(idss[i]), str(topos[i]), str(sizes[i]), str(taxos[i]), str(descs[i])])
		filcsvfiltered.write(rec_gb2.description + "\n")
		SeqIO.write(filtered_genbank[i], fildbfiltered,"fasta")
	else:
		print("entry "+str(i)+" was removed because its size "+str(sizes[i])+" was over the limit "+ str(params.entry_size)+"\n")

fildbfiltered.close()	
	
	
	
	
