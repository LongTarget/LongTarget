import MySQLdb
import fileinput

databaseID=raw_input("Input database ID:")
user=raw_input("Input database user:")
passwd=raw_input("Input database password:")
database=raw_input("Input database name:")

def get_species_id(x):
	if x=='human':
		return 1
	elif x=='rhesus':
		return 2
	elif x=='chimp':
		return 3
	elif x=='gorilla':
		return 4
	elif x=='mouse':
		return 5
	elif x=='guineapig':
		return 6
	elif x=='dog':
		return 7
	elif x=='treeshrew':
		return 8
	elif x=='lemur':
		return 9
	elif x=='tarsier':
		return 10
	elif x=='marmoset':
		return 11
	elif x=='rat':
		return 12
	elif x=='rabbit':
		return 13
	elif x=='cow':
		return 14
	elif x=='elephant':
		return 15
	elif x=='hedgehog':
		return 16
	elif x=='opossum':
		return 17
	elif x=='platypus':
		return 18
	else:
		print "unable to match"
		
exon_id=1325288

for line in fileinput.input("All-genePrediction"):

	line_part_list=line.split( )
	cmalign_UsingSeq_filename=line_part_list[0]+'-'+line_part_list[1]+'.fa'
	all_the_text_list = open(cmalign_UsingSeq_filename).readlines( )    # character string
	x=0
	dict_species_seq={}
	
	while x<len(all_the_text_list):
		species_name=all_the_text_list[x][1:all_the_text_list[x].find('-')]
		dict_species_seq[species_name]=repr(all_the_text_list[x+1])
		x+=2
	row=8
	
	if line_part_list[1].find('.')==-1:
		while row<=112:
			if line_part_list[row] in dict_species_seq:
				sql="""INSERT tb_homology_exon (exon_id,gtf_exon_id_human,start_human,end_human,species_name,start_query,end_query,seqname,strand,start,end,score,species_id,part_num,sequence) VAlUES ("""+str(exon_id)+","+repr(line_part_list[0])+","+line_part_list[6]+","+line_part_list[7]+","+repr(line_part_list[row])+","+line_part_list[row+1]+","+line_part_list[row+2]+","+repr(line_part_list[row+3])+","+repr(line_part_list[row+4])+","+line_part_list[row+5]+","+line_part_list[row+6]+","+line_part_list[row+7]+","+str(get_species_id(line_part_list[row]))+",0"+","+dict_species_seq[line_part_list[row]]+")"
				try:
					cursor.execute(sql)
					db.commit()
				except:
					db.rollback()
			row+=8
		exon_id+=1
	else:
		while row<=112:
			if line_part_list[row] in dict_species_seq:
				sql="""INSERT tb_homology_exon (exon_id,gtf_exon_id_human,start_human,end_human,species_name,start_query,end_query,seqname,strand,start,end,score,species_id,part_num,sequence) VAlUES ("""+str(exon_id)+","+repr(line_part_list[0])+","+line_part_list[6]+","+line_part_list[7]+","+repr(line_part_list[row])+","+line_part_list[row+1]+","+line_part_list[row+2]+","+repr(line_part_list[row+3])+","+repr(line_part_list[row+4])+","+line_part_list[row+5]+","+line_part_list[row+6]+","+line_part_list[row+7]+","+str(get_species_id(line_part_list[row]))+","+line_part_list[1][line_part_list[1].find('.')+1:]+","+dict_species_seq[line_part_list[row]]+")"
				try:
					cursor.execute(sql)
					db.commit()
				except:
					db.rollback()
			row+=8
		if line_part_list[1][line_part_list[1].find('.')+1:]=='1':
			exon_id+=1
			
db.close()