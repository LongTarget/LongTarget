#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<ctype.h>
#include<emmintrin.h>
#include<string>
#include<string.h>
#include<sys/time.h>
#include<iostream>
#include<time.h>
#define BIGNUM 1000000
#define MAXTST 1500
#define MAXLIB 10000
#define EL 125
#define ES 126
#define MIN_RES 1000
#define NA 123
#define MAXSQ 60
#define AA 16807
#define MM 2147483647
#define QQ 127773	
#define RW 2836
#define PI_SQRT6 1.28254983016186409554
#define TINY 1.0e-6
#define MAX_NIT 100
#define MAX_LNCRNA 1000000
using namespace std;
double first_deriv_cen(double lambda, int *sptr, int *n1, int start, int stop, double sumlenL, double cenL, double sumlenH, double cenH);
double second_deriv_cen(double lambda, int *sptr, int *n1, int start, int stop, double sumlenL, double cenL, double sumlenH, double cenH);
double *mle_cen  (int *, int,int *, int, double, double , double );
void st_sort (int *v, int n)
{
  int gap, i, j;
  int tmp;
  double dtmp;
  int w_tmp;
  w_tmp=0;
  for (gap = 1; gap < n/3; gap = 3*gap +1) ;
  for (; gap > 0; gap = (gap-1)/3) 
  {
    for (i = gap; i < n; i++) 
    {
      for (j = i - gap; j >= 0; j -= gap) 
      {
       if (v[j] <= v[j + gap]) break;
       tmp = v[j];
       v[j] = v[j + gap];
       v[j + gap] = tmp;
      }
    }
  }
}

double *mle_cen(int *sptr, int n_len,int *n1, int m_len, double fc,double Lambda, double K_tmp)
{

  double sumlenL, sumlenH, cenL, cenH;
  double sum_s, sum2_s, mean_s, var_s, dtmp;
  int start, stop;
  int i, nf;
	double *wtmpa;
	wtmpa=(double *)malloc(2*sizeof(double));
  int nit = 0;
	double wtmpb=0.0;
  double deriv, deriv2, lambda, old_lambda, sum = 0.0;
  int w_tmp;
  w_tmp=0;
  nf = (fc/2.0) * n_len;
  start = nf;
  stop = n_len - nf;
  st_sort(sptr,n_len);
  sum_s = sum2_s = 0.0;
  for (i=start; i<stop; i++) 
  {
    sum_s += sptr[i];
  }
  dtmp = (double)(stop-start);
  mean_s = sum_s/dtmp;
  for (i=start; i<stop; i++) 
  {
    sum2_s += sptr[i] * sptr[i];
  }
  var_s = sum2_s/(dtmp-1.0);
  sumlenL = sumlenH = 0.0;
  for (i=0; i<start; i++) 
  {
    sumlenL += (double)n1[i];
  }
  for (i=stop; i<n_len; i++) 
  {
    sumlenH += (double)n1[i];
  }

  if (nf > 0) 
  {
    cenL = (double)sptr[start];
    cenH = (double)sptr[stop];
  }
  else 
  {
    cenL = (double)sptr[start]/2.0;
    cenH = (double)sptr[start]*2.0;
  }
  if (cenL >= cenH) 
  {
    printf("cenL is larger than cenH!mle_cen is wrong!\n");
    return NULL;
  }


  lambda = PI_SQRT6/sqrt(var_s);
  if (lambda > 1.0)
  {
    fprintf(stderr," Lambda initial estimate error: lambda: %6.4g; var_s: %6.4g\n",lambda,var_s);
    lambda = 0.2;
  }

  do 
  {
    deriv =   first_deriv_cen(lambda, sptr,n1, start, stop,sumlenL, cenL, sumlenH, cenH);
    deriv2 = second_deriv_cen(lambda, sptr,n1, start, stop,sumlenL, cenL, sumlenH, cenH); 
    old_lambda = lambda;
    if (lambda - deriv/deriv2 > 0.0) lambda = lambda - deriv/deriv2;
    else lambda = lambda/2.0;
    nit++;
  } while (fabs((lambda - old_lambda)/lambda) > TINY && nit < MAX_NIT);



  if (nit >= MAX_NIT) return NULL;
  
  for(i = start; i < stop ; i++) 
  {
    sum += (double) n1[i] * exp(- lambda * (double)sptr[i]);
  }
  wtmpa[0] = lambda;
  wtmpa[1]= (double)n_len/((double)(m_len)*(sum+sumlenL*exp(-lambda*cenL)-sumlenH*exp(-lambda*cenH)));
  return wtmpa;
}
double
first_deriv_cen(double lambda, int *sptr,int *n1, int start, int stop,double sumlenL, double cenL, double sumlenH, double cenH) 
{
  int i;
  double sum = 0.0, sum1 = 0.0, sum2 = 0.0;
  double s, l, es;

  for(i = start ; i < stop ; i++) 
  {
    s = (double)sptr[i];
    l = (double)n1[i];
    es = exp(-lambda * s );
    sum += s;
    sum2 += l * es;
    sum1 += s * l * es;
  }
  sum1 += sumlenL*cenL*exp(-lambda*cenL) - sumlenH*cenH*exp(-lambda*cenH);
  sum2 += sumlenL*exp(-lambda*cenL) - sumlenH*exp(-lambda*cenH);
  return (1.0 / lambda) - (sum /(double)(stop-start)) + (sum1 / sum2);
}

double
second_deriv_cen(double lambda, int *sptr,int *n1, int start, int stop,double sumlenL, double cenL, double sumlenH, double cenH) 
{

  double sum1 = 0.0, sum2 = 0.0, sum3 = 0.0;
  double s, l, es;
  int i;

  for(i = start ; i < stop ; i++) 
  {
    s = (double)sptr[i];
    l = (double)n1[i];
    es = exp(-lambda * s);
    sum2 += l * es;
    sum1 += l * s * es;
    sum3 += l * s * s * es;
  }
  sum1 += sumlenL*cenL*exp(-lambda*cenL) - sumlenH*cenH*exp(-lambda*cenH);
  sum2 += sumlenL*exp(-lambda * cenL) -  sumlenH*exp(-lambda * cenH);
  sum3 += sumlenL*cenL*cenL * exp(-lambda * cenL) -
    sumlenH*cenH*cenH * exp(-lambda * cenH);
  return ((sum1 * sum1) / (sum2 * sum2)) - (sum3 / sum2) - (1.0 / (lambda * lambda));
}
void findmax_score(int *a,int *b,int n)
{
  int i=0;
  int c[500];
  int wtmp;
  for(i=0;i<n;i++)
  {
    if(a[i]>b[i]||a[i]==b[i])
    {
      c[i]=a[i];
    }
    else
    {
      c[i]=b[i];
    }
  }
}
 int nascii[128]={
	EL,NA,NA,NA,NA,NA,NA,NA,NA,NA,EL,NA,NA,EL,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,ES,NA,NA,16,NA,NA,
	NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,ES,NA,NA,ES,NA,
	NA, 1,15, 2,12,NA,NA, 3,13,NA,NA,11,NA, 8,16,NA,
	 6, 7, 6,10, 4, 5,14, 9,17, 7,NA,NA,NA,NA,NA,NA,
	NA, 1,15, 2,12,NA,NA, 3,13,NA,NA,11,NA, 8,16,NA,
	 6, 7, 6,10, 4, 5,14, 9,17, 7,NA,NA,NA,NA,NA,NA};

int npam[450] = {
	 5,					
	-4, 5,						
	-4,-4, 5,				
	-4,-4,-4, 5,				
	-4,-4,-4, 5, 5,				
	 2,-1, 2,-1,-1, 2,				
	-1, 2,-1, 2, 2,-2, 2,				
	 2, 2,-1,-1,-1,-1,-1, 2,		
	 2,-1,-1, 2, 2, 1, 1, 1, 2,		
	-1, 2, 2,-1,-1, 1, 1, 1,-1, 2,			
	-1,-1, 2, 2, 2, 1, 1,-1, 1, 1, 2,		
	 1,-2, 1, 1, 1, 1,-1,-1, 1,-1, 1, 1,		
	 1, 1,-2, 1, 1,-1, 1, 1, 1,-1,-1,-1, 1,		
	 1, 1, 1,-2,-2, 1,-1, 1,-1, 1,-1,-1,-1, 1,	
	-2, 1, 1, 1, 1,-1, 1,-1,-1, 1, 1,-1,-1,-1, 1,	
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1, 
	-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1}; 
struct pstruct{
	int maxlen;
	int pam2[MAXSQ][MAXSQ];
	int dnaseq;
	int pam_h;
	int pam_l;
	int pamoff;
	int have_pam2;
	int nsq;
};
struct m_rand_struct {
	int seed;
};
struct f_struct {
  int 							max_res;
  unsigned char      bias;
  unsigned char *    byte_score;
  void *             workspace;
  int                alphabet_size;
  void *             word_score_memory;
  unsigned short *    word_score;
  void *             byte_score_memory;
  void *             workspace_memory;
  int                try_8bit;
  int                done_8bit;
  int                done_16bit;
} ftr;
struct m_rand_struct *my_srand(int set)	
{
  struct timeval t;
  int n;
  struct m_rand_struct *my_rand_state;

  if ((my_rand_state = (struct m_rand_struct *)calloc(1, sizeof(struct m_rand_struct)))==NULL) 
  {
    fprintf(stderr," *** [my_srand] cannot allocate random state ***\n");
    exit(1);
  } 
  gettimeofday(&t,NULL);
  n = t.tv_usec % 65535;
  if ((n % 2)==0) n++;
  if (set > 0) {  my_rand_state->seed = set;}
  else {my_rand_state->seed = n;}
  my_rand_state->seed=33;
  return my_rand_state;
}
unsigned int my_nrand(int n, struct m_rand_struct *my_rand_state)
{
  unsigned int rn;
  int lo, hi, test;
  hi = my_rand_state->seed / QQ;
  lo = my_rand_state->seed % QQ;
  test = AA * lo - RW * hi;
  if (test > 0) { my_rand_state->seed = test;}
  else {my_rand_state->seed = test + MM;}
  rn = my_rand_state->seed;
  return rn%n;
}
void shuffle(unsigned char *from, unsigned char *to, int n, struct m_rand_struct *rand_state)
{
  int i,j; unsigned char tmp;

  if (from != to) 
  {
    memcpy((void *)to,(void *)from,n);
  }

  for (i=n; i>0; i--) 
  {
    j = my_nrand(i, rand_state);
    tmp = to[j];
    to[j] = to[i-1];
    to[i-1] = tmp;
  }
  to[n] = 0;
}
unsigned char *cg_str(char *str,int *ascii,int *n)
{
	int i=0,j=0;
	unsigned	char *a;
	i=strlen(str);
	a=(unsigned char *)malloc(i*sizeof(unsigned char));
	for(j=0;j<i;j++)
	{
    switch(ascii[str[j]])
    {
      case 1:
      a[j]='\001';break;
      case 2:
      a[j]='\002';break;
      case 3:
      a[j]='\003';break;
      case 4:
      a[j]='\004';break;
      case 5:
      a[j]='\005';break;
      case 16:
      a[j]='\020';break;
      default:
      a[j]='\020';break;
    }
	}
  *n=i;
	return a;
}
void alloc_pam (int d1, int d2, struct pstruct *ppst)
{
  int     i, *d2p;
  char err_str[128];
	ppst->dnaseq=0;
	ppst->pam_h=-1;
	ppst->pam_l=1;
	ppst->pamoff=0;
  ppst->have_pam2 = 1;
}
void
init_pam2 (struct pstruct *ppst,int *aascii)
{
  int     i, j, k, nsq, sa_t;
  int ix_j, ix_l, ix_i, p_i, p_j;
	char pam_sq[]="\0ACGTURYMWSKDHVBNX";
  nsq = ppst->nsq;
	int pam_sq_n=17;
  ppst->pam2[0][0] = -BIGNUM;
  ppst->pam_h = -1; ppst->pam_l = 1;
	int *pam;
	pam=npam;
  k = 0;
	int *pascii;
	pascii=aascii;
  sa_t = nascii['X'];
  for(i=0;i<MAXSQ;i++)
  {
    for(j=0;j<MAXSQ;j++)
    {
      ppst->pam2[i][j]=0;
    }
  }
  for (i = 1; i < sa_t; i++) 
  {
    p_i = pascii[pam_sq[i]];
    ppst->pam2[0][p_i] = ppst->pam2[p_i][0] = -BIGNUM;
    for (j = 1; j <= i; j++) 
    {
      p_j = pascii[pam_sq[j]];
      ppst->pam2[p_j][p_i] = ppst->pam2[p_i][p_j] = pam[k++] - ppst->pamoff;
      if (ppst->pam_l > ppst->pam2[p_i][p_j]) ppst->pam_l = ppst->pam2[p_i][p_j];
      if (ppst->pam_h < ppst->pam2[p_i][p_j]) ppst->pam_h = ppst->pam2[p_i][p_j];
    }
  }
  for (i = sa_t+1; i < pam_sq_n; i++) 
  {
    p_i = pascii[pam_sq[i]];
    ppst->pam2[0][p_i] = ppst->pam2[p_i][0] = -BIGNUM;
  } 
}
void init_work(unsigned char *aa0,int n0,struct pstruct *ppst,struct f_struct **f_arg)
{
  struct f_struct *f_str;
  int e,f,i,j;
  unsigned char *pc;
	unsigned short *ps;
  int ip=0;
  int n_count;
  int col_len;
  int bias,data,nsq=34,overflow;
  f_str=(struct f_struct *)calloc(1,sizeof(struct f_struct));
  f_str->workspace_memory  = (void *)malloc(3*16*(MAXTST+MAXLIB+32)+256);
  f_str->workspace  = (void *) ((((size_t) f_str->workspace_memory) + 255) & (~0xff));
  f_str->word_score_memory = (void *)malloc((n0 + 32) * sizeof (short) * (nsq + 1) + 256);
  f_str->byte_score_memory = (void *)malloc((n0 + 32) * sizeof (char) * (nsq + 1) + 256);
  f_str->word_score = (unsigned short *) ((((size_t) f_str->word_score_memory) + 255) & (~0xff));
  f_str->byte_score = (unsigned char *) ((((size_t) f_str->byte_score_memory) + 255) & (~0xff));
  overflow = 0;
  bias = 127;
  for (i = 1; i < nsq ; i++) 
	{
    for (j = 1; j < nsq ; j++) 
		{
      data = ppst->pam2[i][j];
	    if (data < -128) 
	    {
	      fprintf(stderr,"*** ERROR *** data out of range: %d[%d][%d,%d]\n",data, ip, i, j);
	    }
      if (data < bias) 
			{
        bias = data;
      }
    }
  }
	ps = f_str->word_score;
  col_len = (n0 + 7) / 8;
  n_count = (n0 + 7) & 0xfffffff8;
  for (f = 0; f < n_count; ++f) 
	{
    *ps++ = 0;
  }
  for (f = 1; f < nsq ; f++) 
	{
    for (e = 0; e < col_len; e++) 
		{
      for (i = e; i < n_count; i += col_len) 
			{
        if (i >= n0) 
				{
          data = 0;
        } 
        else 
				{
          data = ppst->pam2[aa0[i]][f];
        }
        *ps++ = (unsigned short)data;
      }
    }
  }
	pc = f_str->byte_score;
  col_len = (n0 + 15) / 16;
  n_count = (n0 + 15) & 0xfffffff0;
  for (f = 0; f < n_count; ++f) 
	{
    *pc++ = 0;
  }
  for (f = 1; f < nsq ; f++) 
	{
    for (e = 0; e < col_len; e++) 
		{
      for (i = e; i < n_count; i += col_len) 
			{
        if (i >= n0) 
				{
          data = -bias;
        }
        else 
				{
          data = ppst->pam2[aa0[i]][f] - bias;
        }
        if (data > 255) 
				{
          printf("Fatal error. data: %d bias: %d, position: %d/%d, Score out of range for 8-bit SSE2 datatype.\n",data, bias, f, e);
          exit(1);
        }
        *pc++ = (unsigned char)data;
      }
    }
  }
  f_str->bias = (unsigned char) (-bias);
  f_str->try_8bit = (overflow == 0) ? 1 : 0;
  f_str->done_8bit  = 0;
  f_str->done_16bit = 0;
  f_str->max_res =3*n0/2>MIN_RES?3*n0/2:MIN_RES;
  *f_arg = f_str;
}
void close_work_f_str(struct f_struct **f_arg)
{
  struct f_struct *f_str;
  f_str=*f_arg;
  if(f_str!=NULL)
  {
    free(f_str->workspace_memory);
    free(f_str->word_score_memory);
    free(f_str->byte_score_memory);
  }
  free(f_str);
  *f_arg=NULL;
}
int smith_waterman_sse2_word(const unsigned char *     query_sequence,
                         unsigned short *    query_profile_word,
                         const int                 query_length,
                         const unsigned char *     db_sequence,
                         const int                 db_length,
                         unsigned short      gap_open,
                         unsigned short      gap_extend,
                         struct f_struct *   f_str)
{
  int     i, j, k;
  short   score;
  int     cmp;
  int     iter = (query_length + 7) / 8;
  __m128i *p;
  __m128i *workspace = (__m128i *) f_str->workspace;
  __m128i E, F, H;
  __m128i v_maxscore;
  __m128i v_gapopen;
  __m128i v_gapextend;
  __m128i v_min;
  __m128i v_minimums;
  __m128i v_temp;
  __m128i *pHLoad, *pHStore;
  __m128i *pE;
  __m128i *pScore;
  v_gapopen = _mm_setzero_si128();	
  v_gapopen = _mm_insert_epi16 (v_gapopen, gap_open, 0);
  v_gapopen = _mm_shufflelo_epi16 (v_gapopen, 0);
  v_gapopen = _mm_shuffle_epi32 (v_gapopen, 0);
  v_gapextend = _mm_setzero_si128();
  v_gapextend = _mm_insert_epi16 (v_gapextend, gap_extend, 0);
  v_gapextend = _mm_shufflelo_epi16 (v_gapextend, 0);
  v_gapextend = _mm_shuffle_epi32 (v_gapextend, 0);
  v_maxscore = _mm_setzero_si128();	
  v_maxscore = _mm_cmpeq_epi16 (v_maxscore, v_maxscore);
  v_maxscore = _mm_slli_epi16 (v_maxscore, 15);
  v_minimums = _mm_shuffle_epi32 (v_maxscore, 0);
  v_min = _mm_shuffle_epi32 (v_maxscore, 0);
  v_min = _mm_srli_si128 (v_min, 14);
  k = 2 * iter;
  p = workspace;
  for (i = 0; i < k; i++)
  {
      _mm_store_si128 (p++, v_maxscore);
  }
  pE = workspace;
  pHStore = pE + iter;
  pHLoad = pHStore + iter;
  for (i = 0; i < db_length; ++i)
  {
    pScore = (__m128i *) query_profile_word + db_sequence[i] * iter;
    F = _mm_setzero_si128();	
    F = _mm_cmpeq_epi16 (F, F);
    F = _mm_slli_epi16 (F, 15);
    H = _mm_load_si128 (pHStore + iter - 1);
    H = _mm_slli_si128 (H, 2);
    H = _mm_or_si128 (H, v_min);
    p = pHLoad;
    pHLoad = pHStore;
    pHStore = p;
    for (j = 0; j < iter; j++)
    {
      E = _mm_load_si128 (pE + j);
      H = _mm_adds_epi16 (H, *pScore++);
      v_maxscore = _mm_max_epi16 (v_maxscore, H);
      H = _mm_max_epi16 (H, E);
      H = _mm_max_epi16 (H, F);
      _mm_store_si128 (pHStore + j, H);
      H = _mm_subs_epi16 (H, v_gapopen);
      E = _mm_subs_epi16 (E, v_gapextend);
      E = _mm_max_epi16 (E, H);
      F = _mm_subs_epi16 (F, v_gapextend);
      F = _mm_max_epi16 (F, H);
      _mm_store_si128 (pE + j, E);
      H = _mm_load_si128 (pHLoad + j);
    }
    j = 0;
    H = _mm_load_si128 (pHStore + j);
    F = _mm_slli_si128 (F, 2);
    F = _mm_or_si128 (F, v_min);
    v_temp = _mm_subs_epi16 (H, v_gapopen);
    v_temp = _mm_cmpgt_epi16 (F, v_temp);
    cmp  = _mm_movemask_epi8 (v_temp);
    while (cmp != 0x0000) 
    {
      E = _mm_load_si128 (pE + j);
      H = _mm_max_epi16 (H, F);
      _mm_store_si128 (pHStore + j, H);
      H = _mm_subs_epi16 (H, v_gapopen);
      E = _mm_max_epi16 (E, H);
      _mm_store_si128 (pE + j, E);
      F = _mm_subs_epi16 (F, v_gapextend);
      j++;
      if (j >= iter)
      {
        j = 0;
        F = _mm_slli_si128 (F, 2);
        F = _mm_or_si128 (F, v_min);
      }
      H = _mm_load_si128 (pHStore + j);
      v_temp = _mm_subs_epi16 (H, v_gapopen);
      v_temp = _mm_cmpgt_epi16 (F, v_temp);
      cmp  = _mm_movemask_epi8 (v_temp);
    }
  }
  v_temp = _mm_srli_si128 (v_maxscore, 8);
  v_maxscore = _mm_max_epi16 (v_maxscore, v_temp);
  v_temp = _mm_srli_si128 (v_maxscore, 4);
  v_maxscore = _mm_max_epi16 (v_maxscore, v_temp);
  v_temp = _mm_srli_si128 (v_maxscore, 2);
  v_maxscore = _mm_max_epi16 (v_maxscore, v_temp);
  score = _mm_extract_epi16 (v_maxscore, 0);
  return score + 32768;
}

int
smith_waterman_sse2_byte(const unsigned char *     query_sequence,
                         unsigned char *     query_profile_byte,
                         const int                 query_length,
                         const unsigned char *     db_sequence,
                         const int                 db_length,
                         unsigned char       bias,
                         unsigned char       gap_open,
                         unsigned char       gap_extend,
                         struct f_struct *   f_str)
{
  int     i, j, k;
  int     score;
  int     dup;
  int     cmp;
  int     iter = (query_length + 15) / 16; 
  __m128i *p;
  __m128i *workspace = (__m128i *) f_str->workspace;
  __m128i E, F, H;
  __m128i v_maxscore;
  __m128i v_bias;
  __m128i v_gapopen;
  __m128i v_gapextend;
  __m128i v_temp;
  __m128i v_zero;
  __m128i *pHLoad, *pHStore;
  __m128i *pE;
  __m128i *pScore;
  dup    = ((short) bias << 8) | bias;
  v_bias = _mm_setzero_si128();
  v_bias = _mm_insert_epi16 (v_bias, dup, 0);
  v_bias = _mm_shufflelo_epi16 (v_bias, 0);
  v_bias = _mm_shuffle_epi32 (v_bias, 0);
  dup  = ((short) gap_open << 8) | gap_open;
  v_gapopen = _mm_setzero_si128();
  v_gapopen = _mm_insert_epi16 (v_gapopen, dup, 0);
  v_gapopen = _mm_shufflelo_epi16 (v_gapopen, 0);
  v_gapopen = _mm_shuffle_epi32 (v_gapopen, 0);
  dup  = ((short) gap_extend << 8) | gap_extend;
  v_gapextend = _mm_setzero_si128();
  v_gapextend = _mm_insert_epi16 (v_gapextend, dup, 0);
  v_gapextend = _mm_shufflelo_epi16 (v_gapextend, 0);
  v_gapextend = _mm_shuffle_epi32 (v_gapextend, 0);
  v_maxscore = _mm_setzero_si128();	
  v_zero = _mm_setzero_si128();	
  k = iter * 2;
  p = workspace;
  for (i = 0; i < k; i++)
  {
      _mm_store_si128 (p++, v_maxscore);
  }
  pE = workspace;
  pHStore = pE + iter;
  pHLoad = pHStore + iter;
  for (i = 0; i < db_length; ++i)
  {
    pScore = (__m128i *) query_profile_byte + db_sequence[i] * iter;
    F = _mm_setzero_si128();
    H = _mm_load_si128 (pHStore + iter - 1);
    H = _mm_slli_si128 (H, 1);
    p = pHLoad;
    pHLoad = pHStore;
    pHStore = p;
    for (j = 0; j < iter; j++)
    {
      E = _mm_load_si128 (pE + j);
      H = _mm_adds_epu8 (H, *pScore++);
      H = _mm_subs_epu8 (H, v_bias);
      v_maxscore = _mm_max_epu8 (v_maxscore, H);
      H = _mm_max_epu8 (H, E);
      H = _mm_max_epu8 (H, F);
      _mm_store_si128 (pHStore + j, H);
      H = _mm_subs_epu8 (H, v_gapopen);
      E = _mm_subs_epu8 (E, v_gapextend);
      E = _mm_max_epu8 (E, H);
      F = _mm_subs_epu8 (F, v_gapextend);
      F = _mm_max_epu8 (F, H);
      _mm_store_si128 (pE + j, E);
      H = _mm_load_si128 (pHLoad + j);
    }
    j = 0;
    H = _mm_load_si128 (pHStore + j);
    F = _mm_slli_si128 (F, 1);
    v_temp = _mm_subs_epu8 (H, v_gapopen);
    v_temp = _mm_subs_epu8 (F, v_temp);
    v_temp = _mm_cmpeq_epi8 (v_temp, v_zero);
    cmp  = _mm_movemask_epi8 (v_temp);
    while (cmp != 0xffff) 
    {
      E = _mm_load_si128 (pE + j);
      H = _mm_max_epu8 (H, F);
      _mm_store_si128 (pHStore + j, H);
      H = _mm_subs_epu8 (H, v_gapopen);
      E = _mm_max_epu8 (E, H);
      _mm_store_si128 (pE + j, E);
      F = _mm_subs_epu8 (F, v_gapextend);
      j++;
      if (j >= iter)
      {
        j = 0;
        F = _mm_slli_si128 (F, 1);
      }
      H = _mm_load_si128 (pHStore + j);
      v_temp = _mm_subs_epu8 (H, v_gapopen);
      v_temp = _mm_subs_epu8 (F, v_temp);
      v_temp = _mm_cmpeq_epi8 (v_temp, v_zero);
      cmp  = _mm_movemask_epi8 (v_temp);
    }
  }
  v_temp = _mm_srli_si128 (v_maxscore, 8);
  v_maxscore = _mm_max_epu8 (v_maxscore, v_temp);
  v_temp = _mm_srli_si128 (v_maxscore, 4);
  v_maxscore = _mm_max_epu8 (v_maxscore, v_temp);
  v_temp = _mm_srli_si128 (v_maxscore, 2);
  v_maxscore = _mm_max_epu8 (v_maxscore, v_temp);
  v_temp = _mm_srli_si128 (v_maxscore, 1);
  v_maxscore = _mm_max_epu8 (v_maxscore, v_temp);
  score = _mm_extract_epi16 (v_maxscore, 0);
  score = score & 0x00ff;
  if (score + bias >= 255)
  {
      score = 255;
  }
  return score;
}

int calc_score(string &strA,string &strB,int dnaStartPos,int rule)
{
  struct pstruct pst;
  struct f_struct f_str;
  struct f_struct *f_str_all[2];
  f_str_all[0]=f_str_all[1]=NULL;
  unsigned char *aa0_save,*aa1_save,*aa0_rc_save,*aa1_shuf_save;
  int i=0,n0=0,n1=0;
  int score=0,rc_score=0;
  int w_char=0;
  char aa0[MAX_LNCRNA];
  char aa1[MAX_LNCRNA];
  char wtmp1,wtmp2;
  char aa0_rc[MAX_LNCRNA];
  int shuf_cnt=0;
  int shuf_max=1002;
  int shuf_save_score[shuf_max];
  int max_shuf_score[500];
  int aa1_len[500];
  int shuf_score[500];
  int shuf_rc_score[500];
  int mle_thresh=0;
  long mle_thresh_tmp=0;
  struct m_rand_struct *rand_state;
  double *mle_rst;
  double lambda_tmp=0.0;
  double K_tmp=0.0;
  aa1_shuf_save=(unsigned char *)malloc(MAX_LNCRNA*sizeof(unsigned char));
  for(w_char=0;w_char<strA.size();w_char++)
  {
    aa0[w_char]=strA[w_char];
  }
  aa0[w_char]='\0';
  for(w_char=0;w_char<strB.size();w_char++)
  {
    aa1[w_char]=strB[w_char];
  }
  aa1[w_char]='\0';
  for(w_char=0;w_char<strlen(aa0);w_char++)
  {
    switch(aa0[strlen(aa0)-w_char-1])
    {
      case 'A':
      aa0_rc[w_char]='T';break;
      case 'T':
      aa0_rc[w_char]='A';break;
      case 'C':
      aa0_rc[w_char]='G';break;
      case 'G':
      aa0_rc[w_char]='C';break;
      case 'U':
      aa0_rc[w_char]='A';break;
      case 'N':
      aa0_rc[w_char]='N';break;
      default:
      aa0_rc[w_char]='N';break;
    }
  }
  aa0_rc[w_char]='\0';
  aa0_save=cg_str(aa0,nascii,&n0);
  aa1_save=cg_str(aa1,nascii,&n1);
  aa0_rc_save=cg_str(aa0_rc,nascii,&n0);
  alloc_pam(MAXSQ,MAXSQ,&pst);
  init_pam2(&pst,nascii);
  init_work(aa0_save,n0,&pst,&f_str_all[0]);
  init_work(aa0_rc_save,n0,&pst,&f_str_all[1]);
  score=smith_waterman_sse2_byte(aa0_save,f_str_all[0]->byte_score,n0,aa1_save,n1,f_str_all[0]->bias,'\020','\004',f_str_all[0]);
  if(score>=255)
  {
    score=smith_waterman_sse2_word(aa0_save,f_str_all[0]->word_score,n0,aa1_save,n1,'\020','\004',f_str_all[0]);
  }
  rc_score=smith_waterman_sse2_byte(aa0_rc_save,f_str_all[1]->byte_score,n0,aa1_save,n1,f_str_all[1]->bias,'\020','\004',f_str_all[1]);
  if(rc_score>=255)
  {
    rc_score=smith_waterman_sse2_word(aa0_rc_save,f_str_all[1]->word_score,n0,aa1_save,n1,'\020','\004',f_str_all[1]);
  }
  int check_max_score=0;
  int check_num=0;
  int tmp_scorea=0,tmp_scoreb=0;
  check_max_score=score>rc_score?score:rc_score;
  rand_state=my_srand(0);
  for(shuf_cnt=0;shuf_cnt<shuf_max;shuf_cnt++)
  {
    shuffle(aa1_save,aa1_shuf_save,n1,rand_state);
    if(shuf_cnt%2==0)
    {
    tmp_scorea=smith_waterman_sse2_byte(aa0_save,f_str_all[0]->byte_score,n0,aa1_shuf_save,n1,f_str_all[0]->bias,'\020','\004',f_str_all[0]);
    if(tmp_scorea>=255)
    {
      tmp_scorea=smith_waterman_sse2_word(aa0_save,f_str_all[0]->word_score,n0,aa1_shuf_save,n1,'\020','\004',f_str_all[0]);
    }
    shuf_save_score[shuf_cnt]=tmp_scorea;
    }
    else
    {
    tmp_scoreb=smith_waterman_sse2_byte(aa0_rc_save,f_str_all[1]->byte_score,n0,aa1_shuf_save,n1,f_str_all[1]->bias,'\020','\004',f_str_all[1]);
    if(tmp_scoreb>=255)
    {
      tmp_scoreb=smith_waterman_sse2_word(aa0_rc_save,f_str_all[1]->word_score,n0,aa1_shuf_save,n1,'\020','\004',f_str_all[1]);
    }
    shuf_save_score[shuf_cnt]=tmp_scoreb;
    }
    if(shuf_cnt<500)
    {
    aa1_len[shuf_cnt]=n1;
    }
  }
  int tmp_num=0;
  int last_score=0;
  for(shuf_cnt=0;shuf_cnt<shuf_max;shuf_cnt+=2)
  {
    if(tmp_num<500)
    {
      max_shuf_score[tmp_num]=shuf_save_score[shuf_cnt]>shuf_save_score[shuf_cnt+1]?shuf_save_score[shuf_cnt]:shuf_save_score[shuf_cnt+1];
      tmp_num++;
    }
    else
    {
      last_score=shuf_save_score[shuf_cnt]>shuf_save_score[shuf_cnt+1]?shuf_save_score[shuf_cnt]:shuf_save_score[shuf_cnt+1];
    }
  }
  max_shuf_score[150]=last_score;
  mle_rst=mle_cen(max_shuf_score,500,aa1_len,n0,0.0,lambda_tmp,K_tmp);
  if(mle_rst==NULL)
  {
    return 0;
  }
  lambda_tmp=mle_rst[0];
       K_tmp=mle_rst[1];
  mle_thresh=(int)((log((double)K_tmp*n1*n0)-log((double)10))/lambda_tmp+0.5);
  close_work_f_str(&f_str_all[0]);
  close_work_f_str(&f_str_all[1]);
  free(mle_rst);
  free(rand_state);
  free(aa0_save);
  free(aa1_save);
  free(aa0_rc_save);
  free(aa1_shuf_save);
  return mle_thresh;
}
