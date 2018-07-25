#!usr/bin/python

import MySQLdb

databaseID=raw_input("Input database ID:")
user=raw_input("Input database user:")
passwd=raw_input("Input database password:")
database=raw_input("Input database name:")

db=MySQLdb.connect(databaseID,user,passwd,database)
cursor=db.cursor()

gene_id_list=range(30265,30277)
for gene_id in gene_id_list:
	sql_1="select start,end,gene_id from tb_myexon where gene_id=%d"  % gene_id
	try:
		cursor.execute(sql_1)
		exon_results=cursor.fetchall()
	except:
		print "Failed to extract data from tb_myexon!!!!!!!!!!!!!!"
		
	exon_dict={}
	for start,end,gene_id in exon_results:
		if gene_id not in exon_dict:
			exon_dict[gene_id]=[(start,end)]
		else:
			exon_dict[gene_id].append((start,end))
			
	for gene_id in exon_dict:
		sql_2="select start,end,strand,sequence from tb_gene where gene_id=%d" % gene_id
		try:
			cursor.execute(sql_2)
			gene_results=cursor.fetchall()
		except:
			print "Failed to extract data from tb_gene!!!!!!!!!!!"
		
		lncRNA_sequence=''
		sequence=gene_results[0][3]
		if gene_results[0][2]=='+':
			for exon_start,exon_end in exon_dict[gene_id]:
				lncRNA_sequence+=sequence[exon_start-gene_results[0][0]+1000:exon_end-gene_results[0][0]+1001]
		else:
			for exon_start,exon_end in exon_dict[gene_id]:
				lncRNA_sequence+=sequence[gene_results[0][1]-exon_start+1000:gene_results[0][1]-exon_end+1001]
		
		sql_3="update tb_gene set lncRNA_sequence=%s where gene_id=%d" % (repr(lncRNA_sequence),gene_id)
		try:
			cursor.execute(sql_3)
			db.commit()
		except:
			db.rollback()
			print "Failed to update tb_gene!!!!!!!!!!!!"
			print sql_3
		
db.close()