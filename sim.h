#include <string.h>
#include <sstream>
#include <fstream>
#include <math.h>
#include <getopt.h>
#include <unistd.h>
#include <map>
#include <algorithm>
#include <ctype.h>
#include <omp.h>
#include <list>
#include<sys/types.h>
#include<dirent.h>
#include<unistd.h>
#include<string>
#include "stats.h"
#define K 50  
#include "rules.h"
using namespace std;
struct triplex 
{
	triplex(){};
	triplex(int n1, int n2, int n3, int n4, int n5, int n6, int n7,int n8, float f1, float f2, float f3, const string &s1, const string &s2,int w1,int w2,int w3,int w4):
	stari(n1), endi(n2), starj(n3), endj(n4), reverse(n5), strand(n6), rule(n7),nt(n8), score(f1), identity(f2), tri_score(f3), stri_align(s1), strj_align(s2), middle(w1),center(w2), motif(w3),neartriplex(w4) {};
  	int stari;
	int endi;
	int starj;
	int endj;
	int reverse;
	int strand;
	int rule;
  	int nt;
	float score;
	float identity;
	float tri_score;
	string stri_align;
	string strj_align;
	int middle;
	int center;
	int motif;
	int neartriplex;
};

typedef struct NODE
{ 
	long  SCORE;
	long  STARI;
	long  STARJ;
	long  ENDI;
	long  ENDJ;
	long  TOP;
	long  BOT;
	long  LEFT;
	long  RIGHT; 
}  vertex, *vertexptr;

struct axis
{
    axis(int n1=0, int n2=0):
        triplexnum(n1), neartriplex(n2) {};              
    int triplexnum;
    int neartriplex;
};

struct tmp_class
{
    tmp_class(){};
    tmp_class(int n1,int n2,int n3,int n4,int n5):genome_start(n1),genome_end(n2),signal_level(n3),peak(n4),row(n5) {};
    int genome_start;
    int genome_end;
    int signal_level;
    int peak;
    int row;
};
float triplex_score(char c1, char c2, int Para)
{
	if (Para>0)
	{
		if     ( c1 == 'T' && c2 == 'T' ) return 3.7;
		else if( c1 == 'A' && c2 == 'G' ) return 2.8;
		else if( c1 == 'C' && c2 == 'G' ) return 2.2;
		else if( c1 == 'C' && c2 == 'T' ) return 2.4;
		else if( c1 == 'C' && c2 == 'C' ) return 4.5;
		else if( c1 == 'G' && c2 == 'T' ) return 2.6;
		else if( c1 == 'G' && c2 == 'C' ) return 2.4;
	}
	else
	{
		if     ( c1 == 'T' && c2 == 'A' ) return 3.0;
		else if( c1 == 'T' && c2 == 'T' ) return 3.5;
		else if( c1 == 'T' && c2 == 'C' ) return 1.0;
		else if( c1 == 'A' && c2 == 'G' ) return 1.0;
		else if( c1 == 'C' && c2 == 'A' ) return 1.0;
		else if( c1 == 'C' && c2 == 'G' ) return 3.0;
		else if( c1 == 'C' && c2 == 'C' ) return 3.0;		
		else if( c1 == 'G' && c2 == 'T' ) return 2.0;		
		else if( c1 == 'G' && c2 == 'C' ) return 1.0;		
	}
	return 0;
}
long addnode(long c, long ci, long cj, long i, long j, vertex  LIST[], long *pnumnode)
{
	short found;			
	long d;
	long  most = 0;
	long  low = 0;
	found = 0;
	for ( d = 0; d < *pnumnode ; d++ )
	{
		most = d;
		if ( LIST[most].STARI == ci && LIST[most].STARJ == cj )
		{
			found = 1;
			break;
		}
	}
	if ( found )
	{
		if ( LIST[most].SCORE < c )
		{
			LIST[most].SCORE = c;
			LIST[most].ENDI = i;
			LIST[most].ENDJ = j;
		}
		if ( LIST[most].TOP > i ) LIST[most].TOP = i;
		if ( LIST[most].BOT < i ) LIST[most].BOT = i;
		if ( LIST[most].LEFT > j ) LIST[most].LEFT = j;
		if ( LIST[most].RIGHT < j ) LIST[most].RIGHT = j;
	}
	else
	{
		if( *pnumnode == K )
		{
			for ( d = 1; d < *pnumnode ; d++ )
				if ( LIST[d].SCORE < LIST[low].SCORE )
					low = d;
			most = low;
		}
		else
			most = (*pnumnode)++;
		LIST[most].SCORE = c;
		LIST[most].STARI = ci;
		LIST[most].STARJ = cj;
		LIST[most].ENDI = i;
		LIST[most].ENDJ = j;
		LIST[most].TOP = LIST[most].BOT = i;
		LIST[most].LEFT = LIST[most].RIGHT = j;
	}
	return 1;
}

int no_cross(vertex  LIST[], long numnode,long  m1,long mm,long n1,long nn,long* prl,long* pcl)
{
	long  cur;
	long i;
	for ( i = 0; i < numnode; i++ )
	{
		cur = i;
		if ( LIST[cur].STARI <= mm && LIST[cur].STARJ <= nn && LIST[cur].BOT >= m1-1 && 
		LIST[cur].RIGHT >= n1-1 && ( LIST[cur].STARI < *prl || LIST[cur].STARJ < *pcl ))
		{
			if ( LIST[cur].STARI < *prl ) *prl = LIST[cur].STARI;
			if ( LIST[cur].STARJ < *pcl ) *pcl = LIST[cur].STARJ;
			break;
		}
	}
	if ( i == numnode )
		return 1;
	else
		return 0;
}

long diff(const char *A, const char *B, long M, long N, long *pI, long *pJ, long tb, long te, long Q, long R, long **psapp, long *plast, long V[][128], list<long> row[], long *CC,long *DD,long *RR,long *SS)
{
	long *sapp;				
	sapp=*psapp;
#define gap(k)  ((k) <= 0 ? 0 : Q+R*(k))	

#define DEL(k)				\
{ (*pI) += k;				\
  if (*plast < 0)				\
    *plast = sapp[-1] -= (k);		\
  else					\
    *plast = *sapp++ = -(k);		\
}
						
#define INS(k)				\
{ (*pJ) += k;				\
  if (*plast < 0)				\
    { sapp[-1] = (k); *sapp++ = *plast; }	\
  else					\
    *plast = *sapp++ = (k);		\
}

						
#define REP 				\
{ *plast = *sapp++ = 0; 			\
}


#define DIAG(ii, jj, x, value)				\
{ for ( tt = 1, it = row[(ii)].begin(); it != row[(ii)].end(); it++ )	\
    if ( *it == (jj) )				\
      { tt = 0; break; }				\
  if ( tt )						\
    x = ( value );					\
}

	long   midi, midj, type;	
	long midc;
	long   i, j;
	long c, e, d, s;
	long t, *va;
	short tt;
	list<long>::iterator it;
	if (N <= 0)
	{
		if (M > 0) DEL(M)
		*psapp=sapp;
		return - gap(M);
	}
	if (M <= 1)
	{
		if (M <= 0)
		{
			INS(N);
			return - gap(N);
		}
		if (tb > te) tb = te;
		midc = - (tb + R + gap(N) );
		midj = 0;
		va = V[A[1]];
		for (j = 1; j <= N; j++)
		{
			for ( tt = 1, it = row[*pI+1].begin(); it != row[*pI+1].end(); it++ )	
				if ( *it == j+(*pJ) )			
					{ tt = 0; break; }		
			if ( tt )			
			{
				c = va[B[j]] - ( gap(j-1) + gap(N-j) );
				if (c > midc)
					{ midc = c; midj = j;}
			}
		}
		if (midj == 0)
			{ INS(N) DEL(1) }
		else
		{
			if (midj > 1) INS(midj-1)
				REP

			(*pI)++; (*pJ)++;
			row[*pI].push_back(*pJ) ;
			if (midj < N) INS(N-midj)
		}
		*psapp=sapp;
		return midc;
    }
	midi = M/2;			
	CC[0] = 0;			
	t = -Q;
	for (j = 1; j <= N; j++)
	{
		CC[j] = t = t-R;
		DD[j] = t-Q;
	}
	t = -tb;
	for (i = 1; i <= midi; i++)
	{
		s = CC[0];
		CC[0] = c = t = t-R;
		e = t-Q;
		va = V[A[i]];
		for (j = 1; j <= N; j++)
		{
			if ((c = c - Q-R) > (e = e - R)) e = c;
			if ((c = CC[j] - Q-R) > (d = DD[j] - R)) d = c;
			DIAG(i+(*pI), j+(*pJ), c, s+va[B[j]])
			if (c < d) c = d;
			if (c < e) c = e;
			s = CC[j];
			CC[j] = c;
			DD[j] = d;
		}
	}
	DD[0] = CC[0];
	RR[N] = 0;		
	t = -Q;			
	for (j = N-1; j >= 0; j--)
	{
		RR[j] = t = t-R;
		SS[j] = t-Q;
	}
	t = -te;
	for (i = M-1; i >= midi; i--)
	{
		s = RR[N];
		RR[N] = c = t = t-R;
		e = t-Q;
		va = V[A[i+1]];
		for (j = N-1; j >= 0; j--)
		{
			if ((c = c - Q-R) > (e = e - R)) e = c;
			if ((c = RR[j] - Q-R) > (d = SS[j] - R)) d = c;
			DIAG(i+1+(*pI), j+1+(*pJ), c, s+va[B[j+1]])
			if (c < d) c = d;
			if (c < e) c = e;
			s = RR[j];
			RR[j] = c;
			SS[j] = d;
		}
	}
	SS[N] = RR[N];
	midc = CC[0]+RR[0];		
	midj = 0;
	type = 1;
	for (j = 0; j <= N; j++)
		if ((c = CC[j] + RR[j]) >= midc)
			if (c > midc || CC[j] != DD[j] && RR[j] == SS[j])
				{ midc = c; midj = j; }
	for (j = N; j >= 0; j--)
		if ((c = DD[j] + SS[j] + Q) > midc)
			{ midc = c; midj = j; type = 2; }

	*psapp=sapp;
	if (type == 1)
	{
		diff(A,B,midi,midj,pI,pJ,tb,Q,Q,R,psapp,plast,V,row,CC,DD,RR,SS);
		diff(A+midi,B+midj,M-midi,N-midj,pI,pJ,Q,te,Q,R,psapp,plast,V,row,CC,DD,RR,SS);
	}
	else
	{
		diff(A,B,midi-1,midj,pI,pJ,tb,0,Q,R,psapp,plast,V,row,CC,DD,RR,SS);
		sapp=*psapp;
		DEL(2);
		*psapp=sapp;
		diff(A+midi+1,B+midj,M-midi-1,N-midj,pI,pJ,0,te,Q,R,psapp,plast,V,row,CC,DD,RR,SS);
	}
	return midc;
}

float display(const char* A, const char* B, long M, long N, long S[], long AP, long BP,string& stri_align,string& strj_align)
{
	long   i, j, f, op, start_i, start_j, match, mis_match;
	string stra, strb;
	float identity;	
	match = mis_match = 0;
	for (i = j = 0; i < M || j < N; ) 
	{
		start_i = i;
		start_j = j;
		while (i < M && j < N && *S == 0) 
		{
			++i;
			++j;
			if (A[i] == B[j])
				++match;
			else
				++mis_match;
			stra+=A[i];
			strb+=B[j];
			S++;
		}
		if (i < M || j < N)
			if ((op = *S++) > 0)
			{
				for(f=0; f<op; f++) {stra+='-';strb+=B[++j];++mis_match;}

			}
			else
			{

				for(f=0; f<-op; f++) {strb+='-';stra+=A[++i];++mis_match;}
			}
	}
	stri_align=stra;
	strj_align=strb;
	identity=(float)(100*match)/(float)(match+mis_match);
	return identity;
}

void SIM(string& strA, string& strB, string& strSrc,long dnaStartPos, long min_score,float parm_M,float parm_I,float parm_O,float parm_E, vector<struct triplex>& triplex_list, long strand, long Para, long rule,int ntMin,int ntMax,int penaltyT,int penaltyC)
{
	short tt;
	long  score; 
	float identity;
	long I, J;			
	long *sapp;
	long  last;
	long  min=0;			
	long V[128][128], Q,R;		
	vertex  LIST[K];			
	long numnode;			
	long  m1, mm, n1, nn;		
	long  rl, cl;			
	short flag;			
	long endi, endj, stari, starj;	
	long count;					
	long  i, j;			
	long  *va;			
	vertex cur; 			
	string stri_align;
	string strj_align;
  	int nt=0;	
  	float final_score=0.0;
	long M, N;
	const char *A, *B;
	string tmpA, tmpB;
	tmpA=' '+strA;
	tmpB=' '+strB;
	A=tmpA.c_str();
	B=tmpB.c_str();
	M=strA.size();
	N=strB.size();
	long CC[N + 1], DD[N + 1];							
	long RR[N + 1], SS[N + 1], EE[N + 1], FF[N + 1];	
	long HH[M + 1], WW[M + 1];							
	long II[M + 1], JJ[M + 1], XX[M + 1], YY[M + 1];	
	long S[N+M+2];										
	list<long> row[M+1];								
	list<long>::iterator it;
	V['A']['A'] = V['C']['C'] = V['G']['G'] = V['T']['T'] = 10 *parm_M;
	V['A']['G'] = V['G']['A'] = V['C']['T'] = V['T']['C'] = 10 *parm_I;
	V['A']['C'] = V['A']['T'] = V['C']['A'] = V['C']['G'] =
	V['G']['C'] = V['G']['T'] = V['T']['A'] = V['T']['G'] = 10 *parm_I;
	Q = -10 * parm_O;
	R = -10 * parm_E;
	numnode = 0;

#define DIAG(ii, jj, x, value)				\
{ for ( tt = 1, it = row[(ii)].begin(); it != row[(ii)].end(); it++ )	\
    if ( *it == (jj) )				\
      { tt = 0; break; }				\
  if ( tt )						\
    x = ( value );					\
}


#define ORDER(ss1, xx1, yy1, ss2, xx2, yy2)		\
{ if ( ss1 < ss2 )					\
    { ss1 = ss2; xx1 = xx2; yy1 = yy2; }		\
  else							\
    if ( ss1 == ss2 )					\
      { if ( xx1 < xx2 )				\
	  { xx1 = xx2; yy1 = yy2; }			\
	else						\
	  if ( xx1 == xx2 && yy1 < yy2 )		\
	    yy1 = yy2;					\
      }							\
}


	long  c;		
	long  f;			
	long  d;			
	long  p;			
	long  ci, cj;		
	long  di, dj;		
	long  fi, fj;		
	long  pi, pj;		
	short  cflag, rflag;		
	long  limit;				
	for ( j = 1; j <= N ; j++ )
	{
		CC[j] = 0;
		RR[j] = 0;
		EE[j] = j;
		DD[j] = - (Q);
		SS[j] = 0;
		FF[j] = j;
	}
	for ( i = 1; i <= M; i++) 
	{
		c = 0;				
		f = - (Q);
		ci = fi = i;
		va = V[A[i]];
		p = 0;
		pi = i - 1;
		cj = fj = pj = 0;
		for ( j = 1 ; j <= N ; j++ )  
		{  
			f = f - R;
			c = c - Q-R;
			ORDER(f, fi, fj, c, ci, cj)
			c = CC[j] - Q-R; 
			ci = RR[j];
			cj = EE[j];
			d = DD[j] - R;
			di = SS[j];
			dj = FF[j];
			ORDER(d, di, dj, c, ci, cj)
			c = 0;
			DIAG(i, j, c, p+va[B[j]])		
			if ( c <= 0 )
			{ c = 0; ci = i; cj = j; }
			else
			{ ci = pi; cj = pj; }
			ORDER(c, ci, cj, d, di, dj)
			ORDER(c, ci, cj, f, fi, fj)
			p = CC[j];
			CC[j] = c;
			pi = RR[j];
			pj = EE[j];
			RR[j] = ci;
			EE[j] = cj;
			DD[j] = d;
			SS[j] = di;
			FF[j] = dj;
			if ( c > min_score )	
				addnode(c, ci, cj, i, j, LIST, &numnode);
		}
	}
	for ( count = numnode - 1; count >= 0 ; count-- )
	{
		
		for ( j = 0, i = 1; i < numnode ; i++ )
			if ( LIST[i].SCORE > LIST[j].SCORE )
				j = i;
		memcpy(&cur, LIST+j, sizeof(vertex));
		if ( j != --numnode )
		{
			memcpy(LIST+j,LIST+numnode, sizeof(vertex));
			memcpy(LIST+numnode, &cur, sizeof(vertex));
		}
		score = cur.SCORE;
		stari = ++cur.STARI;
		starj = ++cur.STARJ;
		endi = cur.ENDI;
		endj = cur.ENDJ;
		m1 = cur.TOP;
		mm = cur.BOT;
		n1 = cur.LEFT;
		nn = cur.RIGHT;
		rl = endi - stari + 1;
		cl = endj - starj + 1;
		I = stari - 1;
		J = starj - 1;
		sapp = S;
		last = 0;
    	nt=endi-stari+1;
		diff(&A[stari]-1, &B[starj]-1,rl,cl,&I,&J,Q,Q,Q,R,&sapp,&last,V,row,CC,DD,RR,SS);
		if( score/10.0 <= min_score )
			break;
		struct triplex atriplex;
		float tri_score = 0.0;
		char  prechar=0, curchar=0;
		identity=display(&A[stari]-1,&B[starj]-1,rl,cl,S,stari,starj,stri_align,strj_align);
		nt=strj_align.size();
    if(strand==0&&(nt>=ntMin&&nt<=ntMax))
		{
			string seqtmp=strSrc;
			if(Para>0)  complement(seqtmp);
			j = 0;
			float hashvalue=0, prescore=0;
			for( i = 0; i < strj_align.size() ; i++ )
			{
				if(strj_align[i]=='-')
				{
					curchar='-';
					hashvalue=triplex_score(curchar, stri_align[i],Para);
				}
				else
				{
					curchar=seqtmp[starj+j-1];
					hashvalue=triplex_score(curchar, stri_align[i], Para);
					j++;
				}
				if( (curchar==prechar) && curchar=='A')
				{
					tri_score=tri_score-prescore+penaltyT;
					hashvalue=penaltyT;
				} 
				if( (curchar==prechar) && curchar=='G') 
				{
					tri_score=tri_score-prescore+penaltyC;
					hashvalue=penaltyC;
				}
				prescore=hashvalue;
				prechar=curchar;
				tri_score+=hashvalue;
			}
      		score/=10;
      		final_score=(float)score/nt;
      		tri_score/=nt;
			atriplex=triplex(stari,endi,starj+dnaStartPos,endj+dnaStartPos,strand,Para,rule,nt,final_score,identity,tri_score,stri_align,strj_align,0,0,0,0);
		}
		else if(strand==1&&(nt>=ntMin&&nt<=ntMax))
		{
			string seqtmp=strSrc;
			if(Para<0) complement(seqtmp);
			j = 0;
			float hashvalue=0, prescore=0;
			for( i = 0; i < strj_align.size() ; i++ )
			{
				if(strj_align[i]=='-')
				{
					curchar='-';
					hashvalue=triplex_score(curchar, stri_align[i], Para);
				}
				else
				{
					curchar=seqtmp[N-starj-j];
					hashvalue=triplex_score(curchar, stri_align[i], Para);
					j++;
				}
				if( (curchar==prechar) && curchar=='A')
				{
					tri_score=tri_score-prescore-1000;
					hashvalue=-1000;
				} 
				if( (curchar==prechar) && curchar=='G') 
				{
					tri_score=tri_score-prescore;
					hashvalue=0;
				}
				prescore=hashvalue;
				prechar=curchar;
				tri_score+=hashvalue;
			}
      		score/=10;
      		final_score=(float)score/nt;
      		tri_score/=nt;
			atriplex=triplex(stari,endi,N-starj+dnaStartPos+1,N-endj+dnaStartPos+1,strand,Para,rule,nt,final_score,identity,tri_score,stri_align,strj_align,0,0,0,0);
		}
		if(nt>=ntMin)
      		triplex_list.push_back(atriplex);
		if ( count )
		{ 
			flag = 0;
			for ( j = nn; j >= n1 ; j-- )
			{
				CC[j] = 0;
				EE[j] = j;
				DD[j] = - (Q);
				FF[j] = j;
				RR[j] = SS[j] = mm + 1;
			}
			for ( i = mm; i >= m1; i-- )
			{
				c = p = 0;
				f = - (Q);
				ci = fi = i;
				pi = i + 1;
				cj = fj = pj = nn + 1;
				va = V[A[i]];
				limit = n1;
				for ( j = nn; j >= limit ; j-- )  
				{
					f = f - R;
					c = c - Q-R;
					ORDER(f, fi, fj, c, ci, cj)
					c = CC[j] - Q-R; 
					ci = RR[j];
					cj = EE[j];
					d = DD[j] - R;
					di = SS[j];
					dj = FF[j];
					ORDER(d, di, dj, c, ci, cj)
					c = 0;
					DIAG(i, j, c, p+va[B[j]])		
					if ( c <= 0 )
					{ c = 0; ci = i; cj = j; }
					else
					{ ci = pi; cj = pj; }
					ORDER(c, ci, cj, d, di, dj)
					ORDER(c, ci, cj, f, fi, fj)
					p = CC[j];
					CC[j] = c;
					pi = RR[j];
					pj = EE[j];
					RR[j] = ci;
					EE[j] = cj;
					DD[j] = d;
					SS[j] = di;
					FF[j] = dj;
					if ( c > min )
						flag = 1;
				}
				HH[i] = CC[n1];
				II[i] = RR[n1];
				JJ[i] = EE[n1];
				WW[i] = f;
				XX[i] = fi;
				YY[i] = fj;
			}
			for ( rl = m1, cl = n1; ; )
			{
				for ( rflag = cflag = 1; ( rflag && m1 > 1 ) || ( cflag && n1 > 1 ) ;  )
				{
					if ( rflag && m1 > 1 )	
					{
						rflag = 0;
						m1--;
						c = p = 0;
						f = - (Q);
						ci = fi = m1;
						pi = m1 + 1;
						cj = fj = pj = nn + 1;
						va = V[A[m1]];
						for ( j = nn; j >= n1 ; j-- )  
						{
							f = f - R;
							c = c - Q-R;
							ORDER(f, fi, fj, c, ci, cj)
							c = CC[j] - Q-R; 
							ci = RR[j];
							cj = EE[j];
							d = DD[j] - R;
							di = SS[j];
							dj = FF[j];
							ORDER(d, di, dj, c, ci, cj)
							c = 0;
							DIAG(m1, j, c, p+va[B[j]])		
							if ( c <= 0 )
							{ c = 0; ci = m1; cj = j; }
							else
							{ ci = pi; cj = pj; }
							ORDER(c, ci, cj, d, di, dj)
							ORDER(c, ci, cj, f, fi, fj)
							p = CC[j];
							CC[j] = c;
							pi = RR[j];
							pj = EE[j];
							RR[j] = ci;
							EE[j] = cj;
							DD[j] = d;
							SS[j] = di;
							FF[j] = dj;
							if ( c > min )
								flag = 1;
							if ( ! rflag && ( ci > rl && cj > cl || di > rl && dj > cl || fi > rl && fj > cl ) )
								rflag = 1;
						}
						HH[m1] = CC[n1];
						II[m1] = RR[n1];
						JJ[m1] = EE[n1];
						WW[m1] = f;
						XX[m1] = fi;
						YY[m1] = fj;
						if ( ! cflag && ( ci > rl && cj > cl || di > rl && dj > cl || fi > rl && fj > cl ) )
							cflag = 1;
					}
					if ( cflag && n1 > 1 )
					{
						cflag = 0;
						n1--;
						c = 0;
						f = - (Q);
						cj = fj = n1;
						va = V[B[n1]];

						p = 0;
						ci = fi = pi = mm + 1;
						pj = n1 + 1;
						limit = mm;
						for ( i = limit; i >= m1 ; i-- )  
						{
							f = f - R;
							c = c - Q-R;
							ORDER(f, fi, fj, c, ci, cj)
							c = HH[i] - Q-R; 
							ci = II[i];
							cj = JJ[i];
							d = WW[i] - R;
							di = XX[i];
							dj = YY[i];
							ORDER(d, di, dj, c, ci, cj)
							c = 0;
							DIAG(i, n1, c, p+va[A[i]])
							if ( c <= 0 )
								{ c = 0; ci = i; cj = n1; }
							else
								{ ci = pi; cj = pj; }
							ORDER(c, ci, cj, d, di, dj)
							ORDER(c, ci, cj, f, fi, fj)
							p = HH[i];
							HH[i] = c;
							pi = II[i];
							pj = JJ[i];
							II[i] = ci;
							JJ[i] = cj;
							WW[i] = d;
							XX[i] = di;
							YY[i] = dj;
							if ( c > min )
								flag = 1;
							if ( ! cflag && ( ci > rl && cj > cl || di > rl && dj > cl || fi > rl && fj > cl ) )
								cflag = 1;
						}
						CC[n1] = HH[m1];
						RR[n1] = II[m1];
						EE[n1] = JJ[m1];
						DD[n1] = f;
						SS[n1] = fi;
						FF[n1] = fj;
						if ( ! rflag && ( ci > rl && cj > cl || di > rl && dj > cl || fi > rl && fj > cl ) )
							rflag = 1;
					}
				}
				if ( m1 == 1 && n1 == 1 || no_cross(LIST,numnode,m1,mm,n1,nn,&rl,&cl) )	 break;
			}
			m1--;
			n1--;
			if ( flag )
			{
				for ( j = n1 + 1; j <= nn ; j++ )
				{
					CC[j] = 0;
					RR[j] = m1;
					EE[j] = j;
					DD[j] = - (Q);
					SS[j] = m1;
					FF[j] = j;
				}
				for ( i = m1 + 1; i <= mm; i++) 
				{
					c = 0;			
					f = - (Q);
					ci = fi = i;
					va = V[A[i]];

					p = 0;
					pi = i - 1;
					cj = fj = pj = n1;
					limit = n1 + 1;
					for ( j = limit ; j <= nn ; j++ )  
					{
						f = f - R;
						c = c - Q-R;
						ORDER(f, fi, fj, c, ci, cj)
						c = CC[j] - Q-R; 
						ci = RR[j];
						cj = EE[j];
						d = DD[j] - R;
						di = SS[j];
						dj = FF[j];
						ORDER(d, di, dj, c, ci, cj)
						c = 0;
						DIAG(i, j, c, p+va[B[j]])		
						if ( c <= 0 )
							{ c = 0; ci = i; cj = j; }
						else
							{ ci = pi; cj = pj; }
						ORDER(c, ci, cj, d, di, dj)
						ORDER(c, ci, cj, f, fi, fj)
						p = CC[j];
						CC[j] = c;
						pi = RR[j];
						pj = EE[j];
						RR[j] = ci;
						EE[j] = cj;
						DD[j] = d;
						SS[j] = di;
						FF[j] = dj;
						if ( c > min )	
							min = addnode(c, ci, cj, i, j, LIST, &numnode);
					}
				}
			}
		}
	}
}

void cluster_triplex(int dd,int length,vector<struct triplex>& triplex_list,map<size_t,size_t> class1[],map<size_t,size_t> class1a[],map<size_t,size_t> class1b[],int class_level)
{
    int i,j;
    int find=0;
    map<size_t, struct axis> axis_map;
    int max_neartriplexnum=0,max_pos=0;
    int middle=0;
    int count=0;
    for(vector<struct triplex>::iterator it=triplex_list.begin();it!=triplex_list.end();it++)
    {
      if(it->nt>length)
      {
        count++;
        middle=(int)((it->stari+it->endi)/2);
        it->middle=middle;
        it->motif=0;
        axis_map[middle].triplexnum++;

        for (i=-dd;i<=dd;i++)
        {
            if(i>0)
            {
            axis_map[middle+i].neartriplex = axis_map[middle+i].neartriplex + (dd-i);
            }
            else if(i<0)
            {
              axis_map[middle+i].neartriplex=axis_map[middle+i].neartriplex+(dd+i);
            }           
            else
            { 
            }
            if(axis_map[middle].triplexnum>0)
            {
            if (axis_map[middle+i].neartriplex >max_neartriplexnum)
              {
                max_neartriplexnum=axis_map[middle+i].neartriplex;
                max_pos=middle+i;
                find=1;
              }
            }
        }
        it->neartriplex=axis_map[middle].neartriplex;
      }
    }
    int theclass=1;
    while( find )
    {
        for(i=max_pos-dd; i<=max_pos+dd; i++)
        {
            for(vector<struct triplex>::iterator it=triplex_list.begin();it!=triplex_list.end();it++)
            {
              if(it->middle==i && it->motif==0)
                {
                    it->motif = theclass;
                    it->center=max_pos;
                    if (theclass > class_level) 
                    {
                      continue;
                    }
                    if(it->endj>it->starj)
                        for(j=it->starj;j<it->endj;j++)
                        {
                          class1[theclass][j]++;
                          class1a[theclass][j]++;
                        }
                    else
                        for(j=it->endj;j<it->starj;j++)
                        {
                          class1[theclass][j]++;
                          class1b[theclass][j]--;
                        }
                }
            }
            axis_map.erase(i);
        }
        max_neartriplexnum=0;
        find=0;
        for(map<size_t, struct axis>::iterator it=axis_map.begin(); it!=axis_map.end(); it++)
        {
            if(it->second.neartriplex>=max_neartriplexnum&&it->second.triplexnum>0)
            {
                max_neartriplexnum=it->second.neartriplex;
                max_pos=it->first;
                find=1;
            }
        }
        ++theclass;
    }
}

void print_cluster(int c_level,map<size_t,size_t> class1[],int start_genome,string &chro_info,int dna_size,string &rna_name,int distance,int length,string &outFilePath,string &c_tmp_dd,string &c_tmp_length,vector<struct tmp_class> &w_tmp_class)
{
  struct tmp_class a_tmp_class;
  char c_level_tmp[3];
  sprintf(c_level_tmp,"%d",c_level);
  string c_tmp_level;
  int c_level_loop=0;
  for(c_level_loop=0;c_level_loop<strlen(c_level_tmp);c_level_loop++)
  {
    c_tmp_level+=c_level_tmp[c_level_loop];
  }
  string class_name=outFilePath.substr(0,outFilePath.size()-10)+"-TFOclass"+c_tmp_level+"-"+c_tmp_dd+"-"+c_tmp_length;
	ofstream outfile(class_name.c_str(),ios::trunc);
  int map_tmp0=0,map_tmp1=0,map_tmp2=0,map_tmp3=0,map_count=0,map_count1=0;
  int map_first1=0,map_second1=0;
  int map_first0=0,map_second0=0;
  int if_map1=0,if_map2=0,if_map3=0,if_map4=0;
  int if_map_flag=0;
  outfile<<"browser position "<<chro_info<<":"<<start_genome<<"-"<<start_genome+dna_size<<endl;
  outfile<<"browser hide all"<<endl;
  outfile<<"browser pack refGene encodeRegions"<<endl;
  outfile<<"browser full altGraph"<<endl;
  outfile<<"# 300 base wide bar graph, ausoScale is on by default == graphing"<<endl;
  outfile<<"# limits will dynamically change to always show full range of data"<<endl;
  outfile<<"# in viewing window, priority = 20 position this as the second graph"<<endl;
  outfile<<"# Note, zero-relative, half-open coordinate system in use for bedGraph format"<<endl;
  outfile<<"track type=bedGraph name='"<<rna_name<<" TTS ("<<c_level<<")' description='"<<distance<<"-"<<length<<"' visibility=full color=200,100,0 altColor=0,100,200 priority=20"<<endl;
  int final_genome=0;
  for(map<size_t,size_t>::iterator it=class1[c_level].begin();it!=class1[c_level].end();it++)
  {
    final_genome=it->first+start_genome;
  }
  for(map<size_t,size_t>::iterator it=class1[c_level].begin();it!=class1[c_level].end();it++)
  {
    map_tmp1=it->first;
    map_tmp2=it->second;
    map_first0=it->first;
    map_second0=it->second;
    if((map_count==0)&&(map_count1==0))
    {
      map_tmp3=it->first;
      map_tmp2=it->second;
      map_count1++;
    }
    if(it->first==final_genome)
    {
      break;
    }
    it++;
    map_first1=it->first;
    map_second1=it->second;
    if(abs((long)(it->first-map_tmp1))==1&&(it->second==map_tmp2))
    {
      map_tmp1=it->first;
      map_tmp2=it->second;
      it++;
      while(abs((long)(it->first-map_tmp1))==1&&(it->second==map_tmp2))
      {
        map_tmp1=it->first;
        map_tmp2=it->second;
        it++;
      }
      if(map_count==0)
      {
      a_tmp_class=tmp_class(map_first0+start_genome-2,map_tmp1+start_genome,map_tmp2,0,0);
      w_tmp_class.push_back(a_tmp_class);
      map_count++;
      }
      else
      {
        a_tmp_class=tmp_class(map_first0+start_genome-1,map_tmp1+start_genome,map_tmp2,0,0);
        w_tmp_class.push_back(a_tmp_class);
      }
    }
    if(abs((long)(it->first-map_tmp1))!=1&&(it->second==map_tmp2))
    {
      a_tmp_class=tmp_class(map_tmp1+start_genome,it->first+start_genome-1,0,0,0);
      w_tmp_class.push_back(a_tmp_class);
      a_tmp_class=tmp_class(it->first+start_genome-1,it->first+start_genome,it->second,0,0);
      w_tmp_class.push_back(a_tmp_class);
    }
    if(abs((long)(it->first-map_tmp1))==1&&(it->second!=map_tmp2))
    {
      a_tmp_class=tmp_class(map_tmp1+start_genome,it->first+start_genome,it->second,0,0);
      w_tmp_class.push_back(a_tmp_class);
    }
    if(abs((long)(it->first-map_tmp1))!=1&&(it->second!=map_tmp2))
    {
      a_tmp_class=tmp_class(map_tmp1+start_genome,it->first+start_genome-1,0,0,0);
      w_tmp_class.push_back(a_tmp_class);
      a_tmp_class=tmp_class(it->first+start_genome-1,it->first+start_genome,it->second,0,0);
      w_tmp_class.push_back(a_tmp_class);
    }
    it--;
  }
  int w_class_loop=0;
  for(w_class_loop=0;w_class_loop<w_tmp_class.size();w_class_loop++)
  {
    tmp_class btc=w_tmp_class[w_class_loop];
    tmp_class ctc=w_tmp_class[w_class_loop+1];
    if(btc.genome_start==final_genome)
    {
      break;
    }
    if(w_class_loop+1==w_tmp_class.size())
    {
    }
    if(btc.genome_start==ctc.genome_start)
    {
      if(1)
      {
        outfile<<chro_info<<"\t"<<btc.genome_start<<"\t"<<ctc.genome_end<<"\t"<<ctc.signal_level<<endl;
      }
      w_class_loop+=1;
    }
    else
    {
      outfile<<chro_info<<"\t"<<btc.genome_start<<"\t"<<btc.genome_end<<"\t"<<btc.signal_level<<endl;

    }
  }
}

