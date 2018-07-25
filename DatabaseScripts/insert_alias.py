#! /usr/bin/python

import MySQLdb

databaseID=raw_input("Input database ID:")
user=raw_input("Input database user:")
passwd=raw_input("Input database password:")
database=raw_input("Input database name:")

db=MySQLdb.connect(databaseID,user,passwd,database)
cursor=db.cursor()

sql_alias="select * from tb_alias"
try:
	cursor.execute(sql_alias)
	alias_result=cursor.fetchall()
except:
	print "Failed to extract data from tb_alias!!!!!!!"

alias_dict={}	
for alias_id,gene_id,gene_name,datasource_id in alias_result:
	if gene_id==22778:
		gene_name+=gene_name+' 26'
	if gene_id not in alias_dict:
		alias_dict[gene_id]=gene_name
	else:
		alias_dict[gene_id]+=(','+gene_name)
		
for gene_id in alias_dict:
	sql_insert="insert tb_alias_new(gene_id,alias_name,datasource_id) value(%d,%s,%d)" % (gene_id,repr(alias_dict[gene_id]),0)
	try:
		cursor.execute(sql_insert)
		db.commit()
	except:
		db.rollback()
		print "Failed to insert!!!!!!!!!!!!!"

db.close()