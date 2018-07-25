import MySQLdb
import os

databaseID=raw_input("Input database ID:")
user=raw_input("Input database user:")
passwd=raw_input("Input database password:")
database=raw_input("Input database name:")

gtf_name='to_begin'
Exon_id=1325288
cmalign_id=0

while Exon_id<=1370626:

	sql1="""SELECT gtf_exon_id from tb_myexon WHERE exon_id ="""+str(Exon_id)
	try:
		cursor.execute(sql1)
		Gtf_exon_id=cursor.fetchall()
	except:
		print"Error:unable to fetch data"
	if gtf_name!=Gtf_exon_id[0][0]:
		n=1
		gtf_name=Gtf_exon_id[0][0]
	else:
		n+=1
		
	file_1=Gtf_exon_id[0][0]+"-exon"+str(n)+".cmalign-merge"
	m=1
	
	if os.path.exists(file_1):
		all_the_text = open(file_1).read( )[:open(file_1).read( ).find('GC SS_cons')]   # character string
		sql2="""INSERT tb_cmalign (cmalign_id,exon_id,part_num,result) VALUES ("""+str(cmalign_id)+','+str(Exon_id)+','+str(m-1)+','+repr(all_the_text)+')'
		try:
			cursor.execute(sql2)
			db.commit()
		except:
			db.rollback()
		cmalign_id+=1
	else:
		file_2=Gtf_exon_id[0][0]+"-exon"+str(n)+'.'+str(m)+".cmalign-merge"
		while os.path.exists(file_2):
			all_the_text = open(file_2).read( )[:open(file_2).read( ).find('GC SS_cons')]    # character string
			sql3="""INSERT tb_cmalign (cmalign_id,exon_id,part_num,result) VALUES ("""+str(cmalign_id)+','+str(Exon_id)+','+str(m)+','+repr(all_the_text)+')'
			try:
				cursor.execute(sql3)
				db.commit()
			except:
				db.rollback()
			cmalign_id+=1
			m+=1
			file_2=Gtf_exon_id[0][0]+"-exon"+str(n)+'.'+str(m)+".cmalign-merge"
			
	Exon_id+=1
	
db.close()

			