#! /usr/bin/python

import MySQLdb
import re
import pdb

databaseID=raw_input("Input database ID:")
user=raw_input("Input database user:")
passwd=raw_input("Input database password:")
database=raw_input("Input database name:")

db=MySQLdb.connect(databaseID,user,passwd,database)
cursor=db.cursor()

def get_species_name(species_id):
	if species_id==1:
		return "human"
	elif species_id==2:
		return "rhesus"
	elif species_id==3:
		return "chimp"
	elif species_id==5:
		return "mouse"
	elif species_id==6:
		return "guineapig"
	elif species_id==7:
		return "dog"
	elif species_id==8:
		return "treeshrew"
	elif species_id==9:
		return "lemur"
	elif species_id==10:
		return "tarsier"
	elif species_id==11:
		return "marmoset"
	elif species_id==12:
		return "rat"
	elif species_id==13:
		return "rabbit"
	elif species_id==14:
		return "cow"
	elif species_id==15:
		return "elephant"
	elif species_id==16:
		return "hedgehog"
	elif species_id==17:
		return "opossum"
	else:
		return "platypus"

def extract_score(exon_id):
	dict={}
	sql_score="""select species_name,sum(score) from tb_homology_exon where exon_id=%d group by species_name""" % exon_id
	try:
		cursor.execute(sql_score)
		for species_name,sum_score in cursor.fetchall():
			if species_name=="macaque":
				species_name="rhesus"
			elif species_name=="mouselemur":
				species_name="lemur"
			dict[species_name]=sum_score
		return dict
		results=cursor.fetchall()
	except:
		print "Something wrong in extract_score function!!!!!!!!!!!!!"
		
def extract_transposon1(exon_id):
	dict={}
	sql_transposon="""select species_id,poson_start_real,poson_end_real from tb_transposon where exon_id=%d""" % exon_id
	try:
		cursor.execute(sql_transposon)
		for species_id,poson_start_real,poson_end_real in cursor.fetchall():
			species_name=get_species_name(species_id)
			if species_name not in dict:
				dict[species_name]=[(poson_start_real,poson_end_real)]
			else:
				dict[species_name].append((poson_start_real,poson_end_real))
	except:
		print "Failed to extract data from tb_transposon!!!!!!!!!!"
	return dict
	
def extract_transposon2(exon_id,part_num):
	dict={}
	sql_transposon="""select t.species_id,poson_start_real,poson_end_real from tb_transposon as t join tb_homology_exon as h  on t.homology_exon_id=h.homology_exon_id \
	where t.exon_id=%d and h.part_num=%d """ % (exon_id,part_num)
	try:
		cursor.execute(sql_transposon)
		for species_id,poson_start_real,poson_end_real in cursor.fetchall():
			species_name=get_species_name(species_id)
			if species_name not in dict:
				dict[species_name]=[(poson_start_real,poson_end_real)]
			else:
				dict[species_name].append((poson_start_real,poson_end_real))
	except:
		print "Failed to extract data from tb_transposon!!!!!!!!!!"
	return dict
	
def extract_transposon3(exon_id):
	sql_transposon="""select poson_start_real,poson_end_real from tb_transposon where exon_id=%d and species_id=1""" % exon_id
	try:
		cursor.execute(sql_transposon)
		return cursor.fetchall()
	except:
		print "Failed to extract data from tb_transposon!!!!!!!!!!"
		
def sort_transposon(transposon_info):
	transposon_INFO=''
	transposon_list=transposon_info.split(';')
	for transposon in transposon_list:
		startend_list=[]
		start_end_list=[]
		if ',' not in transposon:
			transposon_INFO+=transposon+';'
			continue
		transposon_list_a=transposon.split(',')
		for transposon_a in transposon_list_a:
			start,end=transposon_a.split('-')
			start=int(start)
			end=int(end)
			start_end_list.append(start)
			start_end_list.append(end)
			startend_list.append((start,end))
		check_list=[0]*max(start_end_list)
		for start,end in startend_list:
			check_list[start-1:end]=[1]*(end-start+1)
		flag=0
		tran_info=''
		index=0
		length=len(check_list)
		for base in check_list:
			if base==0 and flag==0:
				index+=1
			elif base==1 and flag==0:
				flag=1
				index+=1
				tran_info+=str(index)+'-'
			elif base==0 and flag==1:
				flag=0
				tran_info+=str(index)+','
				index+=1				
			elif base==1 and flag==1:
				index+=1
				if index>=length:
					tran_info+=str(index)
		tran_info=tran_info.strip(',')
		transposon_INFO+=tran_info+';'
	return transposon_INFO.strip(';')
		


myexon_dict={}		
sql_myexon="""select exon_id,gtf_exon_id,exon_number,gene_id from tb_myexon"""
try:
	cursor.execute(sql_myexon)
	for exon_id,gtf_exon_id,exon_number,gene_id in cursor.fetchall():
		myexon_dict[exon_id]=[gtf_exon_id,exon_number,gene_id]
except:
	print "Failed to extract data from tb_myexon!!!!!!!!!!!!"

cmalign_dict={}       #key---exon_id,value----(part_num,result)
sql_cmalign="""select exon_id,part_num,result from tb_cmalign"""
try:
	cursor.execute(sql_cmalign)
	for exon_id,part_num,result in cursor.fetchall():
		if exon_id not in cmalign_dict:
			cmalign_dict[exon_id]=[(part_num,result)]
		else:
			cmalign_dict[exon_id].append((part_num,result))
except:
	print "Failed to extract data from tb_cmalign!!!!!!!!"
gene_id_start=0
species_view_list=["human","chimp","rhesus","marmoset","tarsier","lemur","treeshrew","rabbit","guineapig","mouse","rat","cow","dog","elephant","hedgehog","opossum","platypus"]
	
for exon_id in myexon_dict:

	cmalign_view_dict={}
	gene_id=myexon_dict[exon_id][2]
	species_cmalign_exon_list=[]
	
	if gene_id!=gene_id_start:
		gene_id_start=gene_id
		sql_homo_gene="""select species_id from tb_homology_gene where gene_id=%d """ % gene_id
		try:
			cursor.execute(sql_homo_gene)
			species_list=["human"]
			for species_id, in cursor.fetchall():
				species_list.append(get_species_name(species_id))
		except:
			print "Failed to extract data from tb_homology_gene!!!!!!!!"
	
			
	
	if species_list:
		if len(cmalign_dict[exon_id])==1:
			cmalign_result=cmalign_dict[exon_id][0][1]
			cmalign_result=cmalign_result.strip('\n')
			cmalign_result=cmalign_result.replace('\t',' ')
			cmalign_list=cmalign_result.split('\n')
			for seq in cmalign_list:
				species_name=seq[:seq.find('-')]
				if (species_name not in species_list) and (species_name!='human') :
					continue
				
				cm_seq=seq[seq.rfind(' ')+1:]
				cm_seq=cm_seq.replace('.','-')
				cm_seq=cm_seq.upper()
				cmalign_view_dict[species_name]=cm_seq
		
			####view_result
			human_length=len(cmalign_view_dict["human"])
			for species_name in species_view_list:
				if species_name not in species_list and species_name!="human":
					cmalign_view_dict[species_name]='-'*human_length
				elif species_name in species_list and species_name not in cmalign_view_dict:
					cmalign_view_dict[species_name]='-'*human_length
			mc=0
			while mc<human_length:
				if cmalign_view_dict["human"][mc]=='-':
					for species_name in cmalign_view_dict:
						if cmalign_view_dict[species_name][mc]!='-':
							break
					else:
						for species_name in cmalign_view_dict:
							cmalign_view_dict[species_name]=cmalign_view_dict[species_name][:mc]+'!'+cmalign_view_dict[species_name][mc+1:]
				mc+=1
			
			for species_name in cmalign_view_dict:
				cmalign_view_dict[species_name]=cmalign_view_dict[species_name].replace('!','')
			view_result=''
			for species_name in species_view_list:
				view_result+=cmalign_view_dict[species_name]+'\n'
			view_result=view_result.strip('\n')
			
			####identity			
			species_score_dict=extract_score(exon_id)    #key---specis_name,value-----sum_score
			for species_name in species_score_dict:
				if species_name not in species_list:
					species_score_dict[species_name]=0.0
			
			identity=''
			L=species_score_dict.values()
			if L:
				score_flag=(max(L)+1)*1.2
			else:
				score_flag=1
				
			if "chimp" not in species_score_dict:
				for species_name in species_view_list:
					if species_name=="human":
						continue
					elif species_name in species_score_dict:
						iden=round(species_score_dict[species_name]/score_flag,2)
						if iden>1:
							iden=1.00
						if len(str(iden))==4 or iden==0:							
							identity+=str(iden)+';' 
						elif len(str(iden))==3:
							identity+=str(iden)+'0'+';'
					else:
						identity+="0.0;"
			else:	
				for species_name in species_view_list:
					if species_name=="human":
						continue
					elif species_name=="chimp":
						if species_score_dict["chimp"]==0.0:
							identity+="0.0;"
						else:
							identity+="1.00;"
					elif species_name in species_score_dict:
						if species_score_dict["chimp"]!=0.0:
							iden=round(species_score_dict[species_name]/(species_score_dict["chimp"]),2)
						else:
							iden=round(species_score_dict[species_name]/score_flag,2)
						if iden>1:
							iden=1.00
						if len(str(iden))==4 or iden==0:							
							identity+=str(iden)+';' 
						elif len(str(iden))==3:
							identity+=str(iden)+'0'+';'
					else:
						identity+="0.0;"
			identity=identity.strip(";")
			
			####transposon
			transposon_dict=extract_transposon1(exon_id)
			transposon=''
			for species_name in species_view_list:
				if species_name not in transposon_dict:
					transposon+='0-0;'
				else:
					transposon_check=cmalign_view_dict[species_name].replace('-','')
					if not transposon_check:
						transposon+='0-0;'
					else:
						m=1
						for poson_start_real,poson_end_real in transposon_dict[species_name]:
							transposon_start_seq1=transposon_check[poson_start_real-1:poson_start_real+7]
							transposon_end_seq1=transposon_check[poson_end_real-8:poson_end_real]
							p1=p2=''
							for base in transposon_start_seq1:
								p1+=(base+'-?-*')
							p1=p1.strip('-?-*')
							for base in transposon_end_seq1:
								p2+=(base+'-?-*')
							p2=p2.strip('-?-*')
							transposon_start_gap=re.search(p1,cmalign_view_dict[species_name])
							transposon_start_seq_gap=transposon_start_gap.group()
							flag_start1=cmalign_view_dict[species_name].find(transposon_start_seq_gap)
							transposon_end_gap=re.search(p2,cmalign_view_dict[species_name])
							transposon_end_seq_gap=transposon_end_gap.group()
							flag_end1=cmalign_view_dict[species_name].find(transposon_end_seq_gap)
							gap_number=transposon_end_seq_gap.count('-')
							transposon_view_start=flag_start1+1
							transposon_view_end=flag_end1+8+gap_number
							if transposon_view_start>0 and transposon_view_end>0 and transposon_view_end>transposon_view_start:
								transposon+=(str(transposon_view_start)+'-'+str(transposon_view_end)+',')
								m=0
						if m:
							transposon+='0-0'
						transposon=transposon.strip(',')
						transposon+=';'
			transposon=sort_transposon(transposon.strip(';'))
			sql_insert="""insert tb_cmalign_view(gene_id,gtf_gene_id,exon_id,identity,view_result,exon_number,transposon) value (%d,%s,%d,%s,%s,%d,%s) """ % (myexon_dict[exon_id][2],repr(myexon_dict[exon_id][0]),exon_id,repr(identity),repr(view_result),myexon_dict[exon_id][1],repr(transposon))
			try:
				cursor.execute(sql_insert)
				db.commit()
			except:
				db.rollback()
				print "Failed to insert tb_cmalign_view!!!!!!!!!"
				
				
		else:
			exon_part_length=0
			cmalign_view_DICT={}
			transposon_view_DICT={}
			for species_name in species_view_list:
				cmalign_view_DICT[species_name]=''
				transposon_view_DICT[species_name]=''
			for part_num,result in cmalign_dict[exon_id]:
				cmalign_view_dict={}
				cmalign_result=result.strip('\n')
				cmalign_result=cmalign_result.replace('\t',' ')
				cmalign_list=cmalign_result.split('\n')
				for seq in cmalign_list:
					species_name=seq[:seq.find('-')]
					if (species_name not in species_list) and (species_name!='human'):
						continue
					cm_seq=seq[seq.rfind(' ')+1:]
					cm_seq=cm_seq.replace('.','-')
					cm_seq=cm_seq.upper()
					cmalign_view_dict[species_name]=cm_seq
				
				human_part_length=len(cmalign_view_dict["human"])
				for species_name in species_view_list:
					if species_name not in species_list and species_name!="human":
						cmalign_view_dict[species_name]='-'*human_part_length
					elif species_name in species_list and (species_name not in cmalign_view_dict):
						cmalign_view_dict[species_name]='-'*human_part_length
						
				mc=0
				while mc<human_part_length:
					if cmalign_view_dict["human"][mc]=='-':
						for species_name in cmalign_view_dict:
							if cmalign_view_dict[species_name][mc]!='-':
								break
						else:
							for species_name in cmalign_view_dict:
								cmalign_view_dict[species_name]=cmalign_view_dict[species_name][:mc]+'!'+cmalign_view_dict[species_name][mc+1:]
					mc+=1
				for species_name in cmalign_view_dict:
					cmalign_view_dict[species_name]=cmalign_view_dict[species_name].replace('!','')	
				for species_name in species_view_list:
					cmalign_view_DICT[species_name]+=cmalign_view_dict[species_name]
				
				####transposon
				transposon_dict=extract_transposon2(exon_id,part_num)
				for species_name in species_view_list:
					if species_name not in transposon_dict:
						continue
					else:
						transposon_check=cmalign_view_dict[species_name].replace('-','')
						if not transposon_check:
							continue
						else:
							for poson_start_real,poson_end_real in transposon_dict[species_name]:
								transposon_start_seq1=transposon_check[poson_start_real-1:poson_start_real+7]
								transposon_end_seq1=transposon_check[poson_end_real-8:poson_end_real]
								p1=p2=''
								for base in transposon_start_seq1:
									p1+=(base+'-?-*')
								p1=p1.strip('-?-*')
								for base in transposon_end_seq1:
									p2+=(base+'-?-*')
								p2=p2.strip('-?-*')
								transposon_start_gap=re.search(p1,cmalign_view_dict[species_name])
								transposon_start_seq_gap=transposon_start_gap.group()
								flag_start1=cmalign_view_dict[species_name].find(transposon_start_seq_gap)
								transposon_end_gap=re.search(p2,cmalign_view_dict[species_name])
								transposon_end_seq_gap=transposon_end_gap.group()
								flag_end1=cmalign_view_dict[species_name].find(transposon_end_seq_gap)
								gap_number=transposon_end_seq_gap.count('-')
								transposon_view_start=flag_start1+1
								transposon_view_end=flag_end1+8+gap_number
								if transposon_view_start>0 and transposon_view_end>0 and transposon_view_end>transposon_view_start:
									transposon_view_DICT[species_name]+=str(exon_part_length+transposon_view_start)+'-'+str(exon_part_length+transposon_view_end)+','
				exon_part_length+=human_part_length
			human_transposon_tuple=extract_transposon3(exon_id)
			
			if human_transposon_tuple:
				transposon_check=cmalign_view_DICT["human"].replace('-','')
				for poson_start_real,poson_end_real in human_transposon_tuple:
					transposon_start_seq1=transposon_check[poson_start_real-1:poson_start_real+7]
					transposon_end_seq1=transposon_check[poson_end_real-8:poson_end_real]
					p1=p2=''
					for base in transposon_start_seq1:
						p1+=(base+'-?-*')
					p1=p1.strip('-?-*')
					for base in transposon_end_seq1:
						p2+=(base+'-?-*')
					p2=p2.strip('-?-*')
					transposon_start_gap=re.search(p1,cmalign_view_DICT["human"])
					transposon_start_seq_gap=transposon_start_gap.group()
					flag_start1=cmalign_view_DICT["human"].find(transposon_start_seq_gap)
					transposon_end_gap=re.search(p2,cmalign_view_DICT["human"])
					transposon_end_seq_gap=transposon_end_gap.group()
					flag_end1=cmalign_view_DICT["human"].find(transposon_end_seq_gap)
					gap_number=transposon_end_seq_gap.count('-')
					transposon_view_start=flag_start1+1
					transposon_view_end=flag_end1+8+gap_number
					if transposon_view_start>0 and transposon_view_end>0 and transposon_view_end>transposon_view_start:
						transposon_view_DICT["human"]+=str(transposon_view_start)+'-'+str(transposon_view_end)+','
			else:
				transposon_view_DICT["human"]="0-0"
					
			transposon=''
			for species_name in species_view_list:
				if transposon_view_DICT[species_name]:
					transposon+=transposon_view_DICT[species_name].strip(',')+';'
				else:
					transposon+="0-0;"
				
				
			####view_result
			view_result=''
			for species_name in species_view_list:
				view_result+=cmalign_view_DICT[species_name]+'\n'
			view_result=view_result.strip('\n')
			
			
			###identity
			species_score_dict=extract_score(exon_id)    #key---specis_name,value-----sum_score
			for species_name in species_score_dict:
				if species_name not in species_list:
					species_score_dict[species_name]=0.0
			
			identity=''	
			L=species_score_dict.values()
			if L:
				score_flag=(max(L)+1)*1.2
			else:
				score_flag=1
				
			if "chimp" not in species_score_dict:
				for species_name in species_view_list:
					if species_name=="human":
						continue
					elif species_name in species_score_dict:
						iden=round(species_score_dict[species_name]/score_flag,2)
						if iden>1:
							iden=1.00
						if len(str(iden))==4 or iden==0:							
							identity+=str(iden)+';' 
						elif len(str(iden))==3:
							identity+=str(iden)+'0'+';'
					else:
						identity+="0.0;"
			else:	
				for species_name in species_view_list:
					if species_name=="human":
						continue
					elif species_name=="chimp":
						if species_score_dict["chimp"]!=0.0:
							identity+="1.00;"
						else:
							identity+="0.0;"
					elif species_name in species_score_dict:
						if species_score_dict["chimp"]!=0.0:
							iden=round(species_score_dict[species_name]/(species_score_dict["chimp"]),2)
						else:
							iden=round(species_score_dict[species_name]/score_flag,2)
						if  iden>1:
							iden=1.00
						if len(str(iden))==4 or iden==0:							
							identity+=str(iden)+';' 
						elif len(str(iden))==3:
							identity+=str(iden)+'0'+';'
					else:
						identity+="0.0;"
			identity=identity.strip(";")
			
			transposon=sort_transposon(transposon.strip(';'))
			sql_insert="""insert tb_cmalign_view(gene_id,gtf_gene_id,exon_id,identity,view_result,exon_number,transposon) value (%d,%s,%d,%s,%s,%d,%s) """ % (myexon_dict[exon_id][2],repr(myexon_dict[exon_id][0]),exon_id,repr(identity),repr(view_result),myexon_dict[exon_id][1],repr(transposon))
			try:
				cursor.execute(sql_insert)
				db.commit()
			except:
				db.rollback()
				print "Failed to insert tb_cmalign_view!!!!!!!!!"

		
		

	
			
				
						
			