#! /usr/bin/python

import MySQLdb

host=raw_input("host:")
user=raw_input("user:")
passward=raw_input("passward:")
database=raw_input("database:")
db=MySQLdb.connect(host,user,passward,database)
cursor=db.cursor()

f=open("human-lncRNAs-and-orthologues.txt",'r')
for line in f:
	field_list=line.split('\t')
	species_id=int(field_list[1])
	lncRNA=field_list[2]
	TFO_num=field_list[3]
	TFO_seq=field_list[4]
	TTS=field_list[5]
	TTS=TTS.strip('\r\n')
	
	if species_id==1 and TFO_num=="TFO1":
		sql="update tb_gene set TFO1=%s,TTS1=%s where gene_name=%s" %(repr(TFO_seq),repr(TTS),repr(lncRNA))
		try:
			cursor.execute(sql)
			db.commit()
		except:
			db.rollback()
			print "Failed to insert tb_gene when gene_name="+lncRNA
			
	elif species_id==1 and TFO_num=="TFO2":
		sql="update tb_gene set TFO2=%s,TTS2=%s where gene_name=%s" %(repr(TFO_seq),repr(TTS),repr(lncRNA))
		try:
			cursor.execute(sql)
			db.commit()
		except:
			db.rollback()
			print "Failed to insert tb_gene when gene_name="+lncRNA
			
	elif species_id!=1 and TFO_num=="TFO1":
		sql="update tb_homology_gene  as h join tb_gene as g on h.gene_id=g.gene_id set h.TFO1=%s,h.TTS1=%s where g.gene_name=%s and h.species_id=%d" % (repr(TFO_seq),repr(TTS),repr(lncRNA),int(species_id))
		try:
			cursor.execute(sql)
			db.commit()
		except:
			db.rollback()
			print "Failed to insert tb_gene when gene_name="+lncRNA
			
	elif species_id!=1 and TFO_num=="TFO2":
		sql="update tb_homology_gene as h  join tb_gene as g on h.gene_id=g.gene_id set h.TFO2=%s,h.TTS2=%s where g.gene_name=%s and h.species_id=%d" % (repr(TFO_seq),repr(TTS),repr(lncRNA),int(species_id))
		try:
			cursor.execute(sql)
			db.commit()
		except:
			db.rollback()
			print "Failed to insert tb_gene when gene_name="+lncRNA
			
f.close()
db.close()