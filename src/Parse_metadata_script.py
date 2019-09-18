#!/usr/bin/python
##ARNAUD
import os
import sys

file_list = "PlasmidDB12084_Acclist.txt"
f_file=open(file_list, 'r')
lines = f_file.readlines()
f_file.close()

plasmidId = []
for line in lines:
	line = line.rstrip()
	plasmidId.append(line)

metadata_file=open('PlasmidDB_metadata.txt','r')
#metadata_file=open('test.txt','r')
lines=metadata_file.readlines()
metadata_file.close()

dico_result = {}
field_list = []

for line in lines :

	
	if 'type: source' in line:
		#pId = line.split('\t')[0]
		pId = line.split('\t')[0].split(' ')[0]
		if pId in plasmidId :
			flag = True
			dico_result[pId]={}
			#print(pId)
		else :
			flag = False

	elif flag == True :

		if len(line)>1 and ("location:" not in line) and ("qualifiers:" not in line):

			line = line.rstrip()
			#line = line.replace(' ','')
			#line = line.replace('taxon:','')
			line = line.replace('\t','')
			line = line.replace('Key: ','@@')
			line = line.replace(', Value: ','@@')
			line_list=line.split('@@')
			value = line_list[2].replace('[','').replace(']','').replace("'",'')		
			key = line_list[1]		
			dico_result[pId][key]=value
			if key not in field_list :
				field_list.append(key)	

#print(field_list)
#os.exit()

outfile = open('metadata_tab.tsv','w')
outfile.write('id' + '\t' + '\t'.join(field_list) + "\n")

for pId in dico_result :
	outfile.write(pId)
	for element in field_list:
		if element in dico_result[pId]:
			outfile.write('\t' + dico_result[pId][element])
		else:
			outfile.write('\tNA')	
	outfile.write("\n")	
		
	






