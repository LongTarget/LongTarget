#! /usr/bin/python

import MySQLdb

databaseID=raw_input("Input database ID:")
user=raw_input("Input database user:")
passwd=raw_input("Input database password:")
database=raw_input("Input database name:")

db=MySQLdb.connect(databaseID,user,passwd,database)
cursor=db.cursor()

homo_gene_id=1
while homo_gene_id<=100276:
	sql_1="select sequence from tb_homology_exon where homology_gene_id=%d" % homo_gene_id
	try:
		cursor.execute(sql_1)
		results=cursor.fetchall()
	except:
		print "Failed to extract data from tb_homology_exon!!!!!!!!!!"
	
	sequence=''
	for t in results:
		sequence+=t[0]
	
	sql_2="update tb_homology_gene set lncRNA_sequence=%s where homology_gene_id=%d" % (repr(sequence),homo_gene_id)
	try:
		cursor.execute(sql_2)
		db.commit()
	except:
		db.rollback()
		print sql_2
		
	homo_gene_id+=1
db.close()