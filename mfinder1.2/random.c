/************************************************************************
*
*  File name: random.c
*
*  Description: Random number generator functions
*               Based on sources from Numerical Recepies
*  
*  Copyright © 2002-2004 Weizmann Institute of Science, 
*			   76100 Rehovot Israel, All rights reserved 
*
*************************************************************************/#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "random.h"





#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

#define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

#define IB1 1
#define IB2 2
#define IB5 16
#define IB18 131072

#define NR_END 1


/*************************** Global variables ****************************/
//The random seed
//should be initialized to negative value
//in the beginning of each run 
long Rand_seed; 



/******************************* Functions *******************************/
float ran2(long *idum)
{
	int j;
	long k;
	static long idum2=123456789;
	static long iy=0;
	static long iv[NTAB];
	float temp;

	if (*idum <= 0) {
		if (-(*idum) < 1) *idum=1;
		else *idum = -(*idum);
		idum2=(*idum);
		for (j=NTAB+7;j>=0;j--) {
			k=(*idum)/IQ1;
			*idum=IA1*(*idum-k*IQ1)-k*IR1;
			if (*idum < 0) *idum += IM1;
			if (j < NTAB) iv[j] = *idum;
		}
		iy=iv[0];
	}
	k=(*idum)/IQ1;
	*idum=IA1*(*idum-k*IQ1)-k*IR1;
	if (*idum < 0) *idum += IM1;
	k=idum2/IQ2;
	idum2=IA2*(idum2-k*IQ2)-k*IR2;
	if (idum2 < 0) idum2 += IM2;
	j=iy/NDIV;
	iy=iv[j]-idum2;
	iv[j] = *idum;
	if (iy < 1) iy += IMM1;
	if ((temp=AM*iy) > RNMX) return RNMX;
	else return temp;
}

float ran3(long *idum)
/* Returnes a uniform random deviate between 0.0 and 1.0.
   Set idum to any negative value to initialize or reinitialize
   the sequence */
{
	static int inext,inextp;
	static long ma[56];
	static int iff=0;
	long mj,mk;
	int i,ii,k;

	if (*idum < 0 || iff == 0) {
		iff=1;
		mj=MSEED-(*idum < 0 ? -*idum : *idum);
		mj %= MBIG;
		ma[55]=mj;
		mk=1;
		for (i=1;i<=54;i++) {
			ii=(21*i) % 55;
			ma[ii]=mk;
			mk=mj-mk;
			if (mk < MZ) mk += MBIG;
			mj=ma[ii];
		}
		for (k=1;k<=4;k++)
			for (i=1;i<=55;i++) {
				ma[i] -= ma[1+(i+30) % 55];
				if (ma[i] < MZ) ma[i] += MBIG;
			}
		inext=0;
		inextp=31;
		*idum=1;
	}
	if (++inext == 56) inext=1;
	if (++inextp == 56) inextp=1;
	mj=ma[inext]-ma[inextp];
	if (mj < MZ) mj += MBIG;
	ma[inext]=mj;
	return mj*FAC;
}

int irbit1(unsigned long *iseed)
/* Returnes as an integer a random bit, based on the 18
   low-significance bits in iseed (which is modified for the next
   call) */
{
	unsigned long newbit;	// The accumulated XOR's

	newbit = (*iseed & IB18) >> 17		// get bit 18.
		^ (*iseed & IB5) >> 4			// XOR with bit 5
		^ (*iseed & IB2) >> 1			// XOR with bit 2
		^ (*iseed & IB1);				// XOR with bit 1
	*iseed=(*iseed << 1) | newbit;		// leftshift the seed and put
	return (int) newbit;		// the result of the XOR's in its bit 1
}


//return negative random seed
//this is needed to ran2
void
init_random_seed()
{
	long neg_mask=0x80000000;
	Rand_seed= (long)time(NULL) | neg_mask;
	//Rand_seed= 0xfff00000;
}

// return random integer in the interval 1 to max_val
int
get_rand(int max_val)
{
	float tmp;
	long r;
	
	tmp=ran2(&Rand_seed)*(float)(max_val);

	
	r = ((long)(tmp))+1;
	//just to be sure on the limits (if r2 returns 0.000 or 1.000 exactly 
	if(r==0)
		r=1;
	if(r>max_val)
		r=max_val;
	return (int)(r);
}


double
get_rand_double()
{
	float frnd;

	frnd=ran2(&Rand_seed);
	return((double)frnd);
}




int scale_rand(int num_nodes,double *cumpd)
/* produces a random integer in the range 1:num_nodes according
   to a given  distribution */
{
//	p=min(find(rand<cPk));
	int i;
	double k;

	k=get_rand_double();
	i=1;
	while (cumpd[i]<k)
	{
		i++;
	}
	
	return(i);
	
	
}


/*
* get as input vec of size
*/
void cumsum(int len,double *vec,double *cumvec)
/* returnes the cummulative sum of vec */
{
	int i;
	
	for (i=0;i<len;i++)
		cumvec[i]=0.0;
	for (i=1;i<len;i++)
		cumvec[i]=cumvec[i-1]+vec[i];
}


/*
* init rand prob distribution structure
* arguments:
* rnd_pfd -refernce to struct (will be allocated inside) 
* size-length of vec (if vec [1..n] then length=n)
* pd_vec - vector of probabalities [1..n]
*/
void
init_rnd_pd(Rnd_pd **rnd_pd_p,int size,double *pd_vec)
{
	int i;
	//double n_factor;
	Rnd_pd *rnd_pd;

	rnd_pd=(Rnd_pd*)calloc(1,sizeof(Rnd_pd));
	rnd_pd->size=size;
	//the vector includes two added cells :[0]=0 and [size+1]=1
	rnd_pd->cumpd_vec=(double*)calloc(sizeof(double),size+2);
	for (i=0;i<=size+1;i++)
		rnd_pd->cumpd_vec[i]=0.0;
	for (i=1;i<=size;i++)
		rnd_pd->cumpd_vec[i]=rnd_pd->cumpd_vec[i-1]+pd_vec[i];
	//normalize by dividing by cumpd_vec[size]
	//n_factor=rnd_pd->cumpd_vec[size];

	//for (i=1;i<=size;i++)
		//rnd_pd->cumpd_vec[i]/=n_factor;
	//out of the loop in order to make sure its value is exactly 1
	rnd_pd->cumpd_vec[size]=1.0; //should be anyway but in case that is epsilon away
	rnd_pd->cumpd_vec[size+1]=1.0;
	*rnd_pd_p=rnd_pd;
}

void
free_mem_rnd_pd(Rnd_pd **rnd_pd_p)
{
	if((*rnd_pd_p)!=NULL){
		if((*rnd_pd_p)->cumpd_vec!=NULL)		
			free((*rnd_pd_p)->cumpd_vec);
		free(*rnd_pd_p);
	}
}


/*
* get random number according to some probability distribution
* arguments:
* rnd_pd - Rnd_pd structure already initialized
* returned value:
* the random number (int)
*/
int
get_rand_pd(Rnd_pd *rnd_pd)
{
  int i;
	double k;

	k=get_rand_double();

	i=1;
	while (rnd_pd->cumpd_vec[i]<k)
	{
		i++;
	}
#if 0	
	/* 
	*faster by binary search
	*/
	l=0;
	r=rnd_pd->size+1;
	i=rnd_pd->size/2;
	while(( (rnd_pd->cumpd_vec[i]<=k) || (rnd_pd->cumpd_vec[i-1]>=k) ) && (l<r)) {
		if(k>rnd_pd->cumpd_vec[i]){
			l=i;
			i=(l+r)/2;
		} else {
			r=i;
			i=(l+r)/2;
		}
#if 0
		t++;
		if(t>10) {
			printf("k=%.4f l=%d r=%d i=%d\n",k,l,r,i);
			break;
			}
#endif
	}
#endif
#if 0
	if(i!=j)
		printf("i should be %d, binary search found i=%d\n",j,i);
#endif
	return(i);
}

	


