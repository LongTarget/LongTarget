import MySQLdb
import re

databaseID=raw_input("Input database ID:")
user=raw_input("Input database user:")
passwd=raw_input("Input database password:")
database=raw_input("Input database name:")

db=MySQLdb.connect(databaseID,user,passwd,database)
cursor=db.cursor()

sql_cmalign="select result from tb_cmalign where exon_id=1330999"
try:
	cursor.execute(sql_cmalign)
	cmalign_seq=cursor.fetchall()[0][0]
except:
	print "!!!!!"
	
human_seq=cmalign_seq[:cmalign_seq.find('\n')]
human_seq=human_seq.replace('\t',' ')
human_seq=human_seq[human_seq.find(' ')+1:]

human_seq=human_seq.replace('.','-')
human_check=human_seq.replace('-','')

check_start_seq=human_check[2019105-2018904:2019105-2018904+8]
check_end_seq=human_check[2019105-2018752+1-8:2019105-2018752+1]

p1=p2=''
for base in check_start_seq:
	p1+=(base+"-?-*")
p1=p1.strip("-?-*")
start_gap=re.search(p1,human_seq)
check_start_seq=start_gap.group()
for base in check_end_seq:
	p2+=(base+"-?-*")
p2=p2.strip("-?-*")
end_gap=re.search(p2,human_seq)
check_end_seq=end_gap.group()
poson_start=human_seq.find(check_start_seq)+1
poson_end=human_seq.find(check_end_seq)+8

TFO=str(poson_start)+'-'+str(poson_end)+";0-0"*16

sql_update="update tb_cmalign_view set TFO1=%s where exon_id=1330999" % repr(TFO)
try:
	cursor.execute(sql_update)
	db.commit()
except:
	print "!!!!!!!!!!!!!"