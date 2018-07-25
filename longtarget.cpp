#include "sim.h"
using namespace std;
struct lgInfo
{
  lgInfo(){};
  lgInfo(const string &s1,const string &s2,const string &s3,const string &s4,const string &s5,const string &s6,int s7,const string &s8):lncName(s1),lncSeq(s2),species(s3),dnaChroTag(s4),fileName(s5),dnaSeq(s6),startGenome(s7),resultDir(s8) {};
  string lncName;
  string lncSeq;
  string species;
  string dnaChroTag;
  string fileName;
  string dnaSeq;
  int startGenome;
  string resultDir;
};
struct para
{
  string file1path;
  string file2path;
  string outpath;
  int rule;
  int cutLength;
  int strand;
  int overlapLength;
  int minScore;
  bool detailOutput;
  int ntMin;
  int ntMax;
  float scoreMin;
  float minIdentity;
  float minStability;
  int penaltyT;
  int penaltyC;
  int cDistance;
  int cLength;
};
void cutSequence(string& seq, vector<string>& seqsVec, vector<int>& seqsStartPos, int& cutLength, int& overlapLength,int& cut_num);
void show_help();
void initEnv(int argc,char * const *argv,struct para &paraList);
void LongTarget(struct para &paraList,string rnaSequence,string dnaSequence,vector<struct triplex> &sort_triplex_list);
bool comp(const triplex &a, const triplex &b);
string getStrand(int reverse,int strand);
int same_seq(string &w_str);
void printResult(string &species,struct para paraList,string &lncName,string &dnaFile,vector<struct triplex> &sort_triplex_list,string &chroTag,string &dnaSequence,int start_genome,string &c_tmp_dd,string &c_tmp_length,string &resultDir);
string readDna(string dnaFileName,string &species,string &chroTag,string &startGenome);
string readRna(string rnaFileName,string &lncName);
int main(int argc, char* const* argv)
{
  struct para paraList;
  vector<struct  lgInfo>  lgList;
  initEnv(argc,argv,paraList);
  char c_dd_tmp[10];
  char c_length_tmp[10];
  int c_loop_tmp=0;
  string c_tmp_dd;
  string c_tmp_length;
  sprintf(c_dd_tmp,"%d",paraList.cDistance);
  sprintf(c_length_tmp,"%d",paraList.cLength);
  for(c_loop_tmp=0;c_loop_tmp<strlen(c_dd_tmp);c_loop_tmp++)
  {
    c_tmp_dd+=c_dd_tmp[c_loop_tmp];
  }
  for(c_loop_tmp=0;c_loop_tmp<strlen(c_length_tmp);c_loop_tmp++)
  {
    c_tmp_length+=c_length_tmp[c_loop_tmp];
  }
  string lncName;
  string lncSeq;
  string species;
  string dnaChroTag;
  string fileName;
  string dnaSeq;
  string resultDir;
  string startGenomeTmp;
  long   startGenome;
  dnaSeq=readDna(paraList.file1path,species,dnaChroTag,startGenomeTmp);
  startGenome=atoi(startGenomeTmp.c_str());
  lncSeq=readRna(paraList.file2path,lncName);
  fileName=paraList.file1path.substr(0,paraList.file1path.size()-3);
  lncName.erase(remove(lncName.begin(),lncName.end(),'\r'),lncName.end());
  lncName.erase(remove(lncName.begin(),lncName.end(),'\n'),lncName.end());
  resultDir=paraList.outpath;
  struct lgInfo algInfo;
  algInfo=lgInfo(lncName,lncSeq,species,dnaChroTag,fileName,dnaSeq,startGenome,resultDir);
  lgList.push_back(algInfo);
  int i=0;
  vector<struct triplex> sort_triplex_list;
  LongTarget(paraList,lgList[i].lncSeq,lgList[i].dnaSeq,sort_triplex_list);
  printResult(lgList[i].species,paraList,lgList[i].lncName,lgList[i].fileName,sort_triplex_list,lgList[i].dnaChroTag,lgList[i].dnaSeq,lgList[i].startGenome,c_tmp_dd,c_tmp_length,lgList[i].resultDir);
  cout<<"finished normally"<<endl;
  return 0;
}

void cutSequence(string& seq, vector<string>& seqsVec, vector<int>& seqsStartPos, int& cutLength, int& overlapLength,int &cut_num)
{
	unsigned int pos=0;
  int tmpa=0;
  int tmpb=0;
  seqsVec.clear();seqsStartPos.clear();
	string cutSeq;
  while(pos<seq.size())
	{
		cutSeq = seq.substr(pos,cutLength);
		seqsVec.push_back(cutSeq);
		seqsStartPos.push_back(pos);
		pos += cutLength;
		pos -= overlapLength;
    tmpa++;
	}
  cut_num=tmpa;
}
string readRna(string rnaFileName,string &lncName)
{
  ifstream rnaFile;
  string tmpRNA;
  string tmpStr;
  rnaFile.open(rnaFileName.c_str());
  getline(rnaFile,tmpStr);
  int i=0;
  string tmpInfo;
  for(i=0;i<tmpStr.size();i++)
  {
  	if(tmpStr[i]=='>')
  	{
  		continue;
  	}
  	tmpInfo=tmpInfo+tmpStr[i];
  }
  lncName=tmpInfo;
  cout<<lncName<<endl;
  while(getline(rnaFile,tmpStr))
  {
    tmpRNA=tmpRNA+tmpStr;
  }
  return tmpRNA;
}
string readDna(string dnaFileName,string &species,string &chroTag,string &startGenome)
{
  ifstream dnaFile;
  string tmpDNA;
  string tmpStr;
  dnaFile.open(dnaFileName.c_str());
  getline(dnaFile,tmpStr);
  int i=0;
  int j=0;
  string tmpInfo;
  for(i=0;i<tmpStr.size();i++)
  {
    if(tmpStr[i]=='>')
    {
      continue;
    }
    if(tmpStr[i]=='|' &&j==0)
    {
      species=tmpInfo;
      j++;
      tmpInfo.clear();
      continue;
    }
    if(tmpStr[i]=='|'&&j==1)
    {
      chroTag=tmpInfo;
      j++;
      tmpInfo.clear();
      continue;
    }
    if(tmpStr[i]=='-'&&j==2)
    {
      startGenome=tmpInfo;
      tmpInfo.clear();
      continue;
    }
    tmpInfo=tmpInfo+tmpStr[i];
  }
  cout<<species<<endl;
  cout<<chroTag<<endl;
  cout<<startGenome<<endl;
  while(getline(dnaFile,tmpStr))
  {
    tmpDNA=tmpDNA+tmpStr;
  }
  return tmpDNA;
}

void initEnv(int argc,char * const *argv,struct para &paraList)
{
  const char* optstring = "f:s:r:O:c:m:t:i:S:o:y:z:Y:Z:h:D:E:d";
	struct option long_options[]={
		{"f1",required_argument,NULL,'f'},
		{"f2",required_argument,NULL,'s'},
		{"ni",required_argument,NULL,'y'},
    {"na",required_argument,NULL,'z'},
    {"pc",required_argument,NULL,'Y'},
    {"pt",required_argument,NULL,'Z'},
    {"ds",required_argument,NULL,'D'},
    {"lg",required_argument,NULL,'E'},
    {0,0,0,0}
	};
  paraList.file1path="./";
  paraList.file2path="./";
  paraList.outpath="./";
  paraList.rule=0;
  paraList.cutLength=5000;
  paraList.strand=0;
  paraList.overlapLength=100;
  paraList.minScore=0;
	paraList.detailOutput=false;
  paraList.ntMin=20;
  paraList.ntMax=100000;
  paraList.scoreMin=0.0;
  paraList.minIdentity=60.0;
  paraList.minStability=1.0;
  paraList.penaltyT=-1000;
  paraList.penaltyC=0;
	paraList.cDistance=15;
  paraList.cLength=50;
  int opt;
  if(argc==1)
  {
    show_help();
  }
	while( (opt = getopt_long_only( argc, argv, optstring, long_options, NULL)) != -1 )
	{
		switch(opt)
		{
		case 'f':
			paraList.file1path = optarg;
			break;
		case 's':
			paraList.file2path = optarg;
			break;
		case 'r':
			paraList.rule = atoi(optarg);
			break;
		case 'O':
			paraList.outpath = optarg;
			break;
		case 'c':
			paraList.cutLength = atoi(optarg);
			break;
		case 'm':
			paraList.minScore = atoi(optarg);
			break;
		case 't':
			paraList.strand = atoi(optarg);
			break;
		case 'd':
			paraList.detailOutput = true;
			break;
    case 'i':
      paraList.minIdentity=atoi(optarg);
      break;
    case 'S':
      paraList.minStability=atoi(optarg);
      break;
    case 'y':
      paraList.ntMin=atoi(optarg);
      break;
    case 'z':
      paraList.ntMax=atoi(optarg);
      break;
    case 'Y':
      paraList.penaltyC=atoi(optarg);
      break;
    case 'Z':
      paraList.penaltyT=atoi(optarg);
      break;
    case 'o':
      paraList.overlapLength=atoi(optarg);
      break;
    case 'h':
      show_help();
      break;
    case 'D':
      paraList.cDistance=atoi(optarg);
      break;
    case 'E':
      paraList.cLength=atoi(optarg);
      break;
		}
	}
}
void LongTarget(struct para &paraList,string rnaSequence,string dnaSequence,vector<struct triplex>&sort_triplex_list)
{
	vector< string> dnaSequencesVec;
	vector< int> dnaSequencesStartPos;
  int cut_num=0;
	cutSequence(dnaSequence,dnaSequencesVec,dnaSequencesStartPos,paraList.cutLength,paraList.overlapLength,cut_num);
  vector<struct triplex> triplex_list;
  float Identity=0.0;
  double w_t1,w_t2;
  clock_t time1;
  int minScore=0;
  time1=clock();
	for(int i=0;i<dnaSequencesVec.size();i++)
	{
		long dnaStartPos = dnaSequencesStartPos[i];
	  cout<<"dnaPos="<<dnaStartPos<<endl;
    string seq1=dnaSequencesVec[i];
		if(same_seq(seq1))
    {
      continue;
    }
    if(paraList.strand>=0)
		{
			if(paraList.rule==0)
			{
				for(int j=0;j<6;j++) 
				{
          string seq2=transferString(seq1,0,1,j+1);
					minScore=calc_score(rnaSequence,seq2,dnaStartPos,paraList.rule);
          SIM( rnaSequence, seq2, seq1, dnaStartPos, minScore,5,-4,-12,-4, triplex_list,0,1,j+1,paraList.ntMin,paraList.ntMax,paraList.penaltyT,paraList.penaltyC);
          seq2=transferString(seq1,1,1,j+1);
          reverseSeq(seq2);					
         	minScore=calc_score(rnaSequence,seq2,dnaStartPos,paraList.rule);
          SIM( rnaSequence, seq2, seq1, dnaStartPos, minScore,5,-4,-12,-4, triplex_list,1,1,j+1,paraList.ntMin,paraList.ntMax,paraList.penaltyT,paraList.penaltyC);
        }
			}
			if(paraList.rule>0&&paraList.rule<7)
			{
        string seq2=transferString(seq1,0,1,paraList.rule);
        minScore=calc_score(rnaSequence,seq2,dnaStartPos,paraList.rule);
   			SIM( rnaSequence, seq2, seq1, dnaStartPos, minScore,5,-4,-12,-4, triplex_list,0,1,paraList.rule,paraList.ntMin,paraList.ntMax,paraList.penaltyT,paraList.penaltyC);
        seq2=transferString(seq1,1,1,paraList.rule);
        reverseSeq(seq2);
        minScore=calc_score(rnaSequence,seq2,dnaStartPos,paraList.rule);
        SIM( rnaSequence, seq2, seq1, dnaStartPos, minScore,5,-4,-12,-4, triplex_list,1,1,paraList.rule,paraList.ntMin,paraList.ntMax,paraList.penaltyT,paraList.penaltyC);
      }
		}
		if(paraList.strand<=0)
		{
			if(paraList.rule==0)
			{
				for(int j=0;j<18;j++) 
				{
          string seq2=transferString(seq1,0,-1,j+1);
					minScore=calc_score(rnaSequence,seq2,dnaStartPos,paraList.rule);
          SIM( rnaSequence, seq2, seq1, dnaStartPos, minScore,5,-4,-12,-4, triplex_list,0,-1,j+1,paraList.ntMin,paraList.ntMax,paraList.penaltyT,paraList.penaltyC);
          seq2=transferString(seq1,1,-1,j+1);
          reverseSeq(seq2);
          minScore=calc_score(rnaSequence,seq2,dnaStartPos,paraList.rule);
          SIM( rnaSequence, seq2, seq1, dnaStartPos, minScore,5,-4,-12,-4, triplex_list,1,-1,j+1,paraList.ntMin,paraList.ntMax,paraList.penaltyT,paraList.penaltyC);
        }
			}
			else
			{
        string seq2=transferString(seq1,0,-1,paraList.rule);
        minScore=calc_score(rnaSequence,seq2,dnaStartPos,paraList.rule);
        SIM( rnaSequence, seq2, seq1, dnaStartPos, minScore,5,-4,-12,-4, triplex_list,0,-1,paraList.rule,paraList.ntMin,paraList.ntMax,paraList.penaltyT,paraList.penaltyC);
        seq2=transferString(seq1,1,-1,paraList.rule);
        reverseSeq(seq2);
        minScore=calc_score(rnaSequence,seq2,dnaStartPos,paraList.rule);
        SIM( rnaSequence, seq2, seq1, dnaStartPos, minScore,5,-4,-12,-4, triplex_list,1,-1,paraList.rule,paraList.ntMin,paraList.ntMax,paraList.penaltyT,paraList.penaltyC);
			}
		}
	}
  clock_t time2;
  time2=clock();
  for(int i=0;i<triplex_list.size();i++)
  {
    triplex atr=triplex_list[i];
    if(atr.score>=paraList.scoreMin&&atr.identity>=paraList.minIdentity&&atr.tri_score>=paraList.minStability)
    {
      sort_triplex_list.push_back(atr);
    }
  }
}


void printResult(string &species,struct para paraList,string &lncName,string &dnaFile,vector<struct triplex> &sort_triplex_list,string &chroTag,string &dnaSequence,int start_genome,string &c_tmp_dd,string &c_tmp_length,string &resultDir)
{
  vector<struct tmp_class> w_tmp_class;
  string pre_file2=resultDir+"/"+species+"-"+lncName;
  string pre_file1=dnaFile;
  string outFilePath = pre_file2+"-"+pre_file1+"-TFOsorted";
  ofstream outFile(outFilePath.c_str(),ios::trunc);
  outFile<<"QueryStart\t"<<"QueryEnd\t"<<"StartInSeq\t"<<"EndInSeq\t"<<"Direction\t"<<"StartInGenome\t"<<"EndInGenome\t"<<"MeanStability\t"<<"MeanIdentity(%)\t"<<"Strand\t"<<"Rule\t"<<"Score\t"<<"Nt(bp)\t"<<"Class\t"<<"MidPoint\t"<<"Center\t"<<"TFO sequence"<<endl;
  map<size_t,size_t> class1[6],class1a[6],class1b[6];
  int class_level=5;
  cluster_triplex(paraList.cDistance,paraList.cLength, sort_triplex_list, class1, class1a, class1b, class_level);
  sort(sort_triplex_list.begin(),sort_triplex_list.end(),comp);
  for(int i=0;i<sort_triplex_list.size();i++)
  {
    triplex atr=sort_triplex_list[i];
    if(sort_triplex_list[i].motif==0)
    {
      continue;
    }
    atr.stri_align.erase(remove(atr.stri_align.begin(),atr.stri_align.end(),'-'),atr.stri_align.end());
    if(atr.starj<atr.endj)
      outFile<<atr.stari<<"\t"<<atr.endi<<"\t"<<atr.starj<<"\t"<<atr.endj<<"\t"<<"R\t"<<atr.starj+start_genome-1<<"\t"<<atr.endj+start_genome-1<<"\t"<<atr.tri_score<<"\t"<<atr.identity<<"\t"<<getStrand(atr.reverse,atr.strand)<<"\t"<<atr.rule<<"\t"<<atr.score<<"\t"<<atr.nt<<"\t"<<atr.motif<<"\t"<<atr.middle<<"\t"<<atr.center<<"\t"<<atr.stri_align<<endl;
    else
      outFile<<atr.stari<<"\t"<<atr.endi<<"\t"<<atr.starj<<"\t"<<atr.endj<<"\t"<<"L\t"<<atr.endj+start_genome-1<<"\t"<<atr.starj+start_genome-1<<"\t"<<atr.tri_score<<"\t"<<atr.identity<<"\t"<<getStrand(atr.reverse,atr.strand)<<"\t"<<atr.rule<<"\t"<<atr.score<<"\t"<<atr.nt<<"\t"<<atr.motif<<"\t"<<atr.middle<<"\t"<<atr.center<<"\t"<<atr.stri_align<<endl;
  }
  outFile.close();
  int pr_loop=0;
  for(pr_loop=1;pr_loop<3;pr_loop++)
  {
    print_cluster(pr_loop,class1,start_genome-1,chroTag,dnaSequence.size(),lncName,paraList.cDistance,paraList.cLength,outFilePath,c_tmp_dd,c_tmp_length,w_tmp_class);
  }
  vector<struct tmp_class>tmpClass;
  tmpClass.swap(w_tmp_class);
  for(pr_loop=0;pr_loop<6;pr_loop++)
  {
    class1[pr_loop].clear();
    class1a[pr_loop].clear();
    class1b[pr_loop].clear();
  }
}

bool comp(const triplex &a,const triplex &b)
{
  return a.motif<b.motif;
}
string getStrand(int reverse,int strand)
{
  string Strand;
  if(reverse==0 &&strand==1)
  {
    Strand="ParaPlus";
  }
  else if(reverse==1 &&strand==1)
  {
    Strand="ParaMinus";
  }
  else if(reverse==1 &&strand==-1)
  {
    Strand="AntiMinus";
  }
  else if(reverse==0 &&strand==-1)
  {
    Strand="AntiPlus";
  }
  return Strand;
}

int same_seq(string &w_str)
{
  string A=w_str;
  int i=0;
  int a=0,c=0,g=0,t=0,u=0,n=0;
  for(i=0;i<A.size();i++)
  {
    switch(A[i])
    {
      case 'A':
        a++;
        break;
      case 'C':
        c++;
        break;
      case 'G':
        g++;
        break;
      case 'T':
        t++;
        break;
      case 'U':
        u++;
        break;
      case 'N':
        n++;
        break;
      default:
        cout<<"unknown letter"<<endl;
        break;
    }
  }
  if(a==A.size())
  {
    return 1;
  }
  else if(c==A.size())
  {
    return 1;
  }
  else if(g==A.size())
  {
    return 1;
  }
  else if(t==A.size())
  {
    return 1;
  }
  else if(u==A.size())
  {
    return 1;
  }
  else if(n==A.size())
  {
    return 1;
  }
  else
  {
    return 0;
  }
}

void show_help()
{
  cout<<"This is the help page."<<endl;
  cout<<"options   Parameters      functions"<<endl;
  cout<<"f1   DNA sequence file  used to get the DNA sequence"<<endl;
  cout<<"f2   RNA sequence file  used to get the RNA sequence"<<endl;
  cout<<"r    rules              rules used to construct triplexes.int type.0 is all."<<endl;
  cout<<"O    Output path        if you define this,output result will be in the path.default is pwd"<<endl;
  cout<<"c    Cutlength          Cut sequence's length."<<endl;
  cout<<"m    min_score          Min_score...this option maybe useless.keep it for now."<<endl;
  cout<<"d    detailoutut        if you choose -d option,it will generate a triplex.detail file which describes the sequence-alignment."<<endl;
  cout<<"i    identity           a condition used to pick up triplexes.default is 60.this should be int type such as 60,not 0.6.default is 60."<<endl;
  cout<<"S    stability          a condition like identity,should be float type such as 1.0.default is 1.0."<<endl;
  cout<<"ni   ntmin              triplexes' min length.default is 20."<<endl;
  cout<<"na   ntmax              triplexes' max length.default is 100."<<endl;
  cout<<"pc   penaltyC           penalty about GG.default is 0."<<endl;
  cout<<"pt   penaltyT           penalty about AA.default is -1000."<<endl;
  cout<<"ds   c_dd               distance used by cluster function.default is 15."<<endl;
  cout<<"lg   c_length           triplexes' length threshold used in cluster function.default is 50."<<endl;
  cout<<"all parameters are listed.If you want to run a simple example,type ./LongTarget -f1 DNAseq.fa -f2 RNAseq.fa -r 0 will be OK"<<endl;
  cout<<"any problems or bugs found please send email to us:zhuhao@smu.edu.cn."<<endl;
  exit(1);
}
