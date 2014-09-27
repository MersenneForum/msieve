
/************************************************************************/
/*																		*/
/*	Program to factor numbers using the Quadratic Sieve					*/
/*	Designed for TINY composites; at most 2^60 i.e.double precision 	*/
/*	with RADIX = 2^30.													*/
/*  for now, do not use LP variation; fulls only; may change later		*/
/*	N.B.  if N is greater than (say) (750M)^2 ~ 5.625e+17,  then N is	*/
/*	very unlikely to split into two primes both less than 2^30. So		*/
/*	therefore I don't bother! This gets handled by the lattice siever.	*/
/*	750M may need adjusting.											*/
/*																		*/
/*	Note:  This is designed to work for double precision N only.		*/
/*	sometimes, when N is close to being triple precision, we can only	*/
/*  select 1 as a multiplier and this might be VERY BAD. If so, this	*/
/*  routine which typically takes 2-3 millisec can take 10's of millsec	*/
/*	OR run out of candidate polynomials with A = D^2 being single prec.	*/
/*	if so, we quit.														*/
/*  The code allows 8 more rows than columns. So 1 time out of 256 we	*/
/*  might fail.															*/
/*																		*/
/************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "mpdim.h"
#include "primes.h"
 
#define	DEBUG		0					/* debugging switches			*/
#define DEBUG1		0
#define	DEBUGM		0
#define DEBUG_DIV	0
#define	PROFILE		0
#define	DISPLAY		0
#define DEBUG_NORM	0
#define DISPLAY_GCD	0

#define	SHOW_CANDIDATE	0

#define	TIME		0
#define FBASE		120					/* Max Factor base size			*/
#define SMALLP		29					/* Small prime variation limit	*/
#define	NUMSMALLP	10					/* number of small primes		*/
#define MAX_SIEVE	(3*1024)			/* maximum sieve length			*/
#define	SIEVE_INIT	(1 << 25) + (1 << 17) + (1 << 9) + 2
#define	MASK		0x80808080
#define	LOG			.20	
#define	XTRA		7	
#define	RATE		1600000.0 
#define	LISTSIZE	50
#define	MAXSOL		25
#define	NFACTORS	60
#define	MATSIZE		3					/* ceiling[NFACTORS/32]			*/
#define	MAXFULLS	NFACTORS+XTRA
#define	TOLERANCE	1.00

/************************************************************************/
/*																		*/
/*				G L O B A L      V A R I A B L E S						*/
/*																		*/
/************************************************************************/ 									
 
int factors[FBASE],							/* factor base				*/
    critical,								/* critic size for single p	*/
    pmax;									/* largest factor in base	*/

int sqrtkn[FBASE],							/* SQRT(kN) mod p			*/
    kN_sqrt[MPDIM],							/* SQRT(kn)					*/
    A,										/* Polynomial coefficients	*/
    B,										/* one set of A,B,C,D,DINV	*/
    C,										/* 							*/
    D,										/* SQRT(A) mod p			*/
    DINV[MPDIM],							/* 1/(2D) mod kN			*/
    roots1[FBASE],							/* roots of q(x) = 0 mod p	*/
    roots2[FBASE],							/* roots of q(x) = 0 mod p	*/
	OLD_D;									/* restart value for D		*/

int faclist[FBASE+10][LISTSIZE];			/*	stores prime factors	*/

int indices[FBASE],							/* points in primes array	*/
    prime_pointer;

int	matrix[FBASE+10][3],
	lhs[FBASE+10][MPDIM];
	
int numfull,								/* full factorizations		*/
    ini,									/* loop variables			*/
    il;															  	
 
int nfactors,								/* size of factor base		*/
	numsmallp,
    solution_count,							/* counts # of dependencies	*/
    crossover1,								/*  q(x) becomes negative	*/
    crossover2;
 

double tolerance;		/* tolerance for large prime: pmax^tolerance	*/
 
char test_val,			/* test value for full factorization			*/
     sieve[MAX_SIEVE],	/* byte array for sieving						*/
     logp[FBASE];		/* scaled logs of factor base					*/

double sieve_tot, discrim_tot, select_tot, trial_tot, div_only, roots_tot, scan_tot, norm_tot;
 
static int smallp[] = {1,128,243,125,343,121,169,289,361,529,841,961};

static int log_primes[] = 

{3, 5, 8, 10, 12, 13, 14, 15, 16, 17, 17, 18, 19, 19, 19, 20,
20, 21, 21, 21, 21, 22, 22, 22, 23, 23, 23, 23, 23, 24, 24, 24,
25, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, 27, 27,
27, 27, 27, 27, 27, 28, 28, 28, 28, 28, 28, 28, 28, 28, 29, 29,
29, 29, 29, 29, 29, 29, 29, 29, 30, 30, 30, 30, 30, 30, 30, 30,
30, 30, 30, 30, 30, 30, 31, 31, 31, 31, 31, 31, 31, 31, 31, 31,
31, 31, 31, 31, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32, 32,
32, 32, 32, 32, 32, 32, 32, 32, 32, 33, 33, 33, 33, 33, 33, 33,
33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 33, 34, 34, 34,
34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34, 34,
34, 34, 34, 34, 34, 34, 34, 35, 35, 35, 35, 35, 35, 35, 35, 35,
35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35, 35,
35, 35, 35, 35, 35, 36, 36, 36};


static int fsquared[] = 

{4, 9, 25, 49, 121, 169, 289, 361, 529, 841, 961, 1369, 1681, 1849, 2209, 2809, 3481, 3721, 4489,
5041, 5329, 6241, 6889, 7921, 9409, 10201, 10609, 11449, 11881, 12769, 16129, 17161, 18769, 19321,
22201, 22801, 24649, 26569, 27889, 29929, 32041, 32761, 36481, 37249, 38809, 39601, 44521, 49729,
51529, 52441, 54289, 57121, 58081, 63001, 66049, 69169, 72361, 73441, 76729, 78961, 80089, 85849,
94249, 96721, 97969, 100489, 109561, 113569, 120409, 121801, 124609, 128881, 134689, 139129, 143641,
146689, 151321, 157609, 160801, 167281, 175561, 177241, 185761, 187489, 192721, 196249, 201601,
208849, 212521, 214369, 218089, 229441, 237169, 241081, 249001, 253009, 259081, 271441, 273529,
292681, 299209, 310249, 316969, 323761, 326041, 332929, 344569, 351649, 358801, 361201, 368449,
375769, 380689, 383161, 398161, 410881, 413449, 418609, 426409, 434281, 436921, 452929, 458329,
466489, 477481, 491401, 502681, 516961, 528529, 537289, 546121, 552049, 564001, 573049, 579121,
591361, 597529, 619369, 635209, 654481, 657721, 674041, 677329, 683929, 687241, 703921, 727609,
734449, 737881, 744769, 769129, 776161, 779689, 786769, 822649, 829921, 844561, 863041, 877969,
885481, 896809, 908209, 935089, 942841, 954529, 966289, 982081, 994009, 1018081, 1026169, 1038361,
1042441, 1062961, 1067089, 1079521, 1100401, 1104601, 1125721, 1129969, 1142761, 1181569, 1190281,
1194649, 1203409, 1216609, 1229881, 1247689, 1261129, 1274641, 1324801, 1329409, 1352569, 1371241,
1394761, 1408969, 1423249, 1442401, 1471369, 1481089, 1495729};

static unsigned int bits[] = 
	{ 020000000000, 010000000000, 004000000000, 002000000000,
	  001000000000, 000400000000, 000200000000, 000100000000,
	  000040000000, 000020000000, 000010000000, 000004000000,
	  000002000000, 000001000000, 000000400000, 000000200000,
	  000000100000, 000000040000, 000000020000, 000000010000,
	  000000004000, 000000002000, 000000001000, 000000000400,
	  000000000200, 000000000100, 000000000040, 000000000020,
	  000000000010, 000000000004, 000000000002, 000000000001 };

void write_num(), get_mult(), mul_single(),mpsqrt(), comp_fact(),
	 qs_get_logs(),discrim_roots(),get_fact(),reverse_startpts(),
	 read_params(), qs_sieve_init(),qs_siever(),single_subt(),
	 get_roots(),sieve_scan(),out_put(),single_add(),
	 knuth(),output(),single_modinv_dp(),reduce_matrix(),single_power(),
	 init_vars(), check_solution();

extern int div_single(),mod_mult_asm(),mpcmp(),smod(),prp2(),next_prime(),spower(),
		   single_modinv(),gcd();

extern void add(),subt(),mult(),divv(),square(),mpcopy(),modinverse(),mpower(),mpctoi(),
			mul_single_prec_asm(),square_dble_prec_asm();

int scan_it(),new_modinv(),msqrt(),mpresidue(),isresidue(), prp_single(), refer_index(),
	compute_gcd(), dump_matrix(), dump_row(), finish(), dump_id(), select_coeffs(),
    jacobi_symbol(), get_next_prime(), vsmall_power();

extern __int64 prectime();

#define setbit(row, index)		row[index >> 5] ^= bits[index & 31]
#define bittest(row)		  ((row[current_word] & bits[current_bit]) != 0)


typedef unsigned __int64 uint64;
int dc, dc0;
 
/************************************************************************/
/*																		*/
/*							M A I N										*/ 
/*																		*/
/************************************************************************/
 
small_qs(N, fac1, fac2)
int N[], *fac1, *fac2;
 
{   /* start of small_qs */

register int i,status;
int k,kN[MPDIM];
char get_mainlog();
double t1;

init_vars(); 
solution_count = 0;
numfull = 		 0;
nfactors =		NFACTORS;
tolerance =		TOLERANCE;
sieve_tot = 0.0;
select_tot = 0.0;
trial_tot = 0.0;
div_only = 0.0;
roots_tot = 0.0;
scan_tot = 0.0;
dc = 0;
dc0 = 0;

get_mult(&k, N);
if (DEBUG1) (void) printf("multiplier = %d\n",k);
mul_single(N,k,kN);
if (DISPLAY) {(void) printf("qs trying: "); write_num(kN, numstr);} 
if (DEBUG1) write_num(kN, numstr);

mpsqrt(kN,kN_sqrt);								/* [SQRT(kN)]			*/
if (TIME) t1 =  (double) prectime();
comp_fact(NFACTORS,factors,kN);					/* compute factor base	*/
if (TIME) (void) printf("comp_fact time = %g\n", ((double)prectime() - t1)/RATE);

numsmallp = 0;
for (i=1; i<=NUMSMALLP; i++) if (factors[i] <= SMALLP) numsmallp++;

qs_get_logs(NFACTORS,factors,logp);				/* and their logs		*/
if (TIME) t1 = (double) prectime();
discrim_roots(kN, NFACTORS,sqrtkn);				/* SQRT(kN) mod p 		*/
if (TIME) (void) printf("discrim root time = %g\n", ((double)prectime() - t1)/RATE);

pmax = factors[NFACTORS];					   /* get largest factor	*/
prime_pointer = indices[NFACTORS];

if (DEBUG1) (void) printf("pmax = %d\n",pmax);

OLD_D = primes[prime_pointer];
if (DEBUG1) (void) printf("OLD_D = %d\n", OLD_D);
while ((OLD_D & 3) == 1)					/* search for p = 3 mod 4   */
	{
    prime_pointer++;
	OLD_D = primes[prime_pointer];
	if (DEBUG1) (void) printf("new old_d = %d\n",OLD_D);
	}
if (DEBUG1) (void) printf("OLD_D = %d\n", OLD_D);

/*	residues are fully factored when their sieved values exceed			*/
/*	.5*log(kN)-tolerance*log(pmax). Compute this value					*/

test_val = get_mainlog(kN, tolerance, pmax);
if (DEBUG) (void) printf("sieve comparison value = %d\n",test_val);

/*	For each time we use a different polynomial: Ax^2 + Bx + C			*/
/*	with A = D^2. Select these coefficients								*/
/*	We must have D == -1 Mod 4, and (D/kN) = 1 and D about equal		*/
/*	to sqrt(kN/2)/MAX_SIEVE so residues will be small.					*/
/*  For N ~ 2^62 and sieve length 5K,  D ~ 1K							*/

select_coeffs(kN);

/*	for each polynomial q(x) we must find its roots mod p for each		*/
/*	p in the factor base. There are 2 for each prime p					*/

get_roots(kN);							/* roots mod p[i]				*/

while(1)
	{
	qs_sieve_init(); 
	qs_siever(kN, 1);

	if (numfull >= NFACTORS + XTRA)
		{
		if (numfull > MAXFULLS) numfull = MAXFULLS;
		if (PROFILE)
			{
			(void) printf("Sieve total = %lf\n",sieve_tot);
			(void) printf("select total = %lf\n",select_tot);
			(void) printf("trial_tot = %lf\n",trial_tot);
			(void) printf("div_only = %lf\n",div_only);
			(void) printf("roots tot = %lf\n",roots_tot);
			(void) printf("scan tot = %lf\n",scan_tot);
			(void) printf("norm tot = %lf\n",norm_tot);
			(void) printf("dc0, dc = %d %d\n",dc0, dc);
			}
		return(finish(N, fac1, fac2));
		}

/*			Now fill the sieve in the other direction (for -x)			*/
/*					And sieve the array again							*/
	qs_sieve_init();
	reverse_startpts();
	qs_siever(kN, -1);

 
/*		Check to see if we have enough factorizations (exit)			*/
 
	if (numfull >= NFACTORS + XTRA)
		{
		if (numfull > MAXFULLS) numfull = MAXFULLS;
		if (PROFILE)
			{
			(void) printf("Sieve total = %lf\n",sieve_tot);
			(void) printf("select total = %lf\n",select_tot);
			(void) printf("trial_tot = %lf\n",trial_tot);
			(void) printf("div_only = %lf\n",div_only);
			(void) printf("roots tot = %lf\n",roots_tot);
			(void) printf("scan tot = %lf\n",scan_tot);
			(void) printf("norm tot = %lf\n",norm_tot);
			(void) printf("dc0, dc = %d %d\n",dc0, dc);
			}
		return(finish(N, fac1, fac2));
		}
 
	OLD_D = D;								/* set up more polynomials	*/
	status = select_coeffs(kN);				/* get more coefficients	*/
	if (status == -1)						/* can't get more; try to	*/
		{									/* finish if we can or quit	*/
		if (numfull >= NFACTORS-2)  return(finish(N, fac1, fac2));
		else return(-1);
		}
	get_roots(kN);							/* and roots of the polys	*/

	}   /* end of qs while(1) loop */

return(0);

}   /* end of small_qs */

/************************************************************************/
/*																		*/
/*		Routine to initialize the sieve array							*/
/*																		*/
/************************************************************************/

void qs_sieve_init()

{   /* start of qs_sieve_init */
int i,count,*sieve_init;


//  memset faster?

sieve_init = (int*) sieve;
count = MAX_SIEVE >> 6;
if (factors[1] == 2)
   for (i=0; i<count; i++)
	{
	*(sieve_init++) = SIEVE_INIT; *(sieve_init++) = SIEVE_INIT;
	*(sieve_init++) = SIEVE_INIT; *(sieve_init++) = SIEVE_INIT;
	*(sieve_init++) = SIEVE_INIT; *(sieve_init++) = SIEVE_INIT;
	*(sieve_init++) = SIEVE_INIT; *(sieve_init++) = SIEVE_INIT;
	*(sieve_init++) = SIEVE_INIT; *(sieve_init++) = SIEVE_INIT;
	*(sieve_init++) = SIEVE_INIT; *(sieve_init++) = SIEVE_INIT;
	*(sieve_init++) = SIEVE_INIT; *(sieve_init++) = SIEVE_INIT;
	*(sieve_init++) = SIEVE_INIT; *(sieve_init++) = SIEVE_INIT;
	}
else
   for (i=0; i<count; i++)
	{
	*(sieve_init++) = 0; *(sieve_init++) = 0;
	*(sieve_init++) = 0; *(sieve_init++) = 0;
	*(sieve_init++) = 0; *(sieve_init++) = 0;
	*(sieve_init++) = 0; *(sieve_init++) = 0;
	*(sieve_init++) = 0; *(sieve_init++) = 0;
	*(sieve_init++) = 0; *(sieve_init++) = 0;
	*(sieve_init++) = 0; *(sieve_init++) = 0;
	*(sieve_init++) = 0; *(sieve_init++) = 0;
	}

}   /* end of qs_sieve_init */

/************************************************************************/
/*																		*/
/*		SIEVER															*/
/*  new version;  do both roots in same loop                            */
/*																		*/
/************************************************************************/

void qs_siever(kN,sign)
int sign;
int kN[];

{   /* start of siever */
int i, pos, pos1, prime;
register char logval;
int fourp, stop;
char *sieve1, *sieve2, *sieve3;
double t1;

if (TIME) t1 = (double)prectime();

for (i=1; i<=numsmallp; i++)		/* first root of q(x) = 0 mod p	*/
   	{							 	/* small primes first			*/
   	pos = min(roots1[i],roots2[i]);
	pos1 = max(roots1[i], roots2[i]); /* 2nd root                   */
   	logval = logp[i];
   	prime = smallp[indices[i]+1];
   	fourp = prime << 2;
   	stop = MAX_SIEVE - fourp + prime; 
   	sieve1 = sieve + prime;
   	sieve2 = sieve1 + prime;
   	sieve3 = sieve2 + prime;
   	while (pos1 < stop)
		{	
		sieve[pos] += logval;  sieve[pos1] += logval;
		sieve1[pos] += logval; sieve1[pos1] += logval;
		sieve2[pos] += logval; sieve2[pos1] += logval;
		sieve3[pos] += logval; sieve3[pos1] += logval;
		pos += fourp;
		pos1 += fourp;
		}
   	while (pos1 < MAX_SIEVE)
		{
		sieve[pos] += logval; sieve[pos1] += logval;
		pos += prime;
		pos1 += prime;  
		}
	}

for (; i<=NFACTORS; i++)				/* now the larger primes	*/
   	{
	pos = min(roots1[i],roots2[i]);
	pos1 = max(roots1[i], roots2[i]);
   	logval = logp[i];
   	prime = factors[i];

	fourp = prime << 2;
	stop = MAX_SIEVE - fourp + prime;
	sieve1 = sieve + prime;
	sieve2 = sieve1 + prime;
	sieve3 = sieve2 + prime;
	while (pos1 < stop)
	   	{
		sieve[pos] += logval;  sieve[pos1] += logval;
		sieve1[pos] += logval; sieve1[pos1] += logval;
		sieve2[pos] += logval; sieve2[pos1] += logval;
		sieve3[pos] += logval; sieve3[pos1] += logval;
		pos += fourp;
	   	pos1 += fourp;
	   	}
	while (pos1 < MAX_SIEVE)
		{
		sieve[pos] += logval; sieve[pos1] += logval;
		pos += prime;
		pos1 += prime;
		}
	}

if (TIME) 
	{
	t1 = ((double)prectime() - t1)/RATE;
	sieve_tot += t1;
	(void) printf("siever time = %g\n",t1); 
	}

sieve_scan(kN,sign);

}   /* end of siever */
		
/************************************************************************/
/*																		*/
/*	routine to scan the sieve array after it has been sieved			*/
/*	we identify completely factored residues by large remaining			*/
/*	values and produce their exact factorization via trial division		*/
/*	then output the exponent vector to the data file along with the		*/
/*	value of q(x).														*/
/*																		*/
/*	The value of 'sign' is 1 or -1 and indicates whether the sieve		*/
/*	was done over a positive or negative interval						*/
/*	NOTE!  fix e[] array --> mechanism of storing successful factors	*/
/*		and output routine;												*/
/*																		*/
/************************************************************************/

void sieve_scan(kN,sign)
int sign;
int kN[];

{   /* start of sieve_scan */
register int index,limit,count;
int signbit,rem,j;
int temp1[MPDIM],temp3[MPDIM],qx[MPDIM],temp[MPDIM];
int number[MPDIM],list[LISTSIZE];
int success[1000];
double t1,t2, t3, t4;
unsigned __int64 candidate;

if (TIME) t1 = (double)prectime();
if (TIME) t3 = (double)prectime();

count = scan_it(success);			/* find success locations		*/
scan_tot += ((double)prectime() - t3)/RATE;


for (j=0; j<LISTSIZE; j++) list[j] = 0;
if (DEBUG_NORM) (void) printf("after scan, count = %d\n",count);

if (count == 0) return;

for (limit=0; limit<count; limit++)				/* sieve scan loop	*/
   {
	if (TIME) t4 = (double)prectime(); 
   index = success[limit];
   if (DEBUG) (void) printf("hit at %d\n",index);
   signbit = 0;
   SIZE(temp) = SINGLE;
   LASTDIGIT(temp) = A;
   mul_single(temp, index, temp1);    // need  c = a*b  in assembler; slightly faster
   if (sign == 1)
	  single_add(temp1, B, temp3);	 // replace 
   else
	  single_subt(temp1, B, temp3);  

   if (DEBUG_NORM) printf("A,B, index = %d %d %d\n",A,B, index);
   if (DEBUG1) { (void) printf("%d * %d + %d = ",A,index, B); write_num(temp3,numstr); }

   mult(DINV, temp3, temp1);					/* (Ai + B)*D^-1	*/

   if (DEBUG_NORM) {printf("(Ai+B)/D = "); write_num(temp1,numstr);  }

   divv(temp1,kN,temp3,qx);						/* qx^2=Ai^2 + 2Bi + C	*/
   												/* MOD kN				*/
   if (DEBUG_NORM) { printf("qx = "); write_num(qx,numstr); }

/*			We now have qx*qx == Ai^2 + 2Bi + C MOD kN					*/
/*			Compute qx*qx and factor it via trial division				*/
 
   square_dble_prec_asm(qx,temp3);				/*	qx*qx				*/
   divv(temp3,kN,temp,number);					/* number = Q*Q mod kN	*/
   j = index*sign;
   if (j >= crossover1 && j <= crossover2)
		{										/* Poly is negative		*/
   		subt(kN,number,number);
   		signbit = 1;							/* add -1 as a factor	*/
   		if (DEBUG1) { printf("crossover = "); write_num(number,numstr); }
		}
   if (SHOW_CANDIDATE) {(void) printf("factoring: "); write_num(number,numstr);}

/*				OK , now lets trial divide								*/

if (TIME)
	{
	norm_tot += ((double)prectime() - t4)/RATE;
	t2 = (double)prectime();
	}

candidate = (uint64)FIRST(number)*(uint64)1073741824  + (uint64)LASTDIGIT(number);


if (SHOW_CANDIDATE) (void) printf("candidate = %I64d\n",candidate);

list[0] = 0;
for (j=1; j<=NFACTORS; j++) 
	{
	rem = (int) (candidate % factors[j]); 
	while (rem == 0)
			{
			if (DEBUG_DIV) (void) printf("divide at: %d\n",factors[j]);
			list[0]++;
			list[list[0]] = j; 
			candidate /= factors[j];  
			if (DEBUG_DIV) printf("new candidate = %I64d\n",candidate);
			if ((uint64)fsquared[indices[j]] > candidate)		/*  at most one factor remains */
				{ 
				if (DEBUG)  (void) printf("RESIDUE = %ld\n",candidate);
   				if (candidate <= (uint64)pmax)					/* full factorization	       */
					{ 
					index = refer_index((int)candidate); 
					list[0]++; 
					list[list[0]] = index;
					output(signbit,qx,NFACTORS,list);
					numfull++;
					goto outer_loop;
   					}
				}
			rem = (int) (candidate % factors[j]);
			}

	}								/* end of while j loop				*/
outer_loop:;

  if (TIME) div_only += ((double)prectime() - t2)/RATE;
  }									/* end of sieve scanning loop		*/

if (TIME) 
	{
	t1 =  ( (double)prectime() - t1)/RATE;
	trial_tot += t1;
	(void) printf("sieve scan/trial time = %g\n",t1);
	}

}   /* end of sieve_scan */
 
 
/************************************************************************/
/*																		*/
/*	routine to find roots of ax^2 + bx + c = 0 mod p					*/
/*																		*/
/************************************************************************/
 
void get_roots(kN)
int kN[];
 
{   /* start of get_roots */
 
register int  j, p, r1, r2;
int  cmodp, bmodp, amodp, A_inv;
int temp[MPDIM], temp1[MPDIM];
double t1;

if (TIME) t1 = (double)prectime();

for (j=1; j<=NFACTORS; j++)
	{
	p = factors[j];
	if (j <= numsmallp) p = smallp[indices[j]+1];
	amodp = A % p;								/* reduce A,		*/
	bmodp = B % p;								/* B	  mod p		*/
	if (DEBUG1) (void) printf("A mod, B mod p = %d %d %d\n",amodp,bmodp,p);

/*  unless amodp == 0, in which case there is only					*/
/*	one root: X such that bmodp X + C == 0 mod p					*/

	if (amodp == 0)
	   { 
	   if (DEBUG1) (void) printf("single root: A, A mod p = 0 %d %d\n",A, p);
	   SIZE(temp1) = SINGLE;
	   LASTDIGIT(temp1) = B;
	   square_dble_prec_asm(temp1,temp);
	   subt(kN,temp,temp);
	   div_single(temp, A, temp1);				/* C = (kN-B^2)/A	*/
	   cmodp = div_single(temp1, p, temp);
	   A_inv = single_modinv(bmodp,p);
 
	   r1 = mod_mult_asm(A_inv,cmodp,p);
	   r2 = r1;
	   }
	else
		{		 
		A_inv = single_modinv(amodp,p);  
		r1 = ((-bmodp - sqrtkn[j])*A_inv) % p;     //  modmultadd_asm?
		if (r1 < 0) r1 += p;
		r2 = ((sqrtkn[j] - bmodp)*A_inv) % p;
		if (r2 < 0) r2 += p;
		}

	roots1[j] = r1;
	roots2[j] = r2;
	if (DEBUG)
		{
		(void) printf("root1[%d] = %d\n",p,r1);
		(void) printf("root2[%d] = %d\n",p,r2);
		fflush(stdout);
		}
	}
if (TIME) 
	{
	t1 =  ((double)prectime() - t1)/RATE;
	roots_tot += t1;
	(void) printf("get roots time = %g\n",t1);
	}
}   /* end of get_roots */
 

/************************************************************************/
/*																		*/
/*			routine to find roots of kN mod p							*/
/*																		*/
/************************************************************************/
 
void discrim_roots(kN, nfactors,sqrtkn)
int nfactors,
    sqrtkn[],
	kN[];

{   /* start of discrim_roots */
 
register int i,modp,test;
int temp[MPDIM];
 
for (i=1; i<=NFACTORS; i++)
	{
	test = factors[i];
	if (i <= numsmallp) test = smallp[indices[i]+1];
	modp = div_single(kN,test,temp);				/* get kN mod p		*/
	sqrtkn[i] = msqrt(modp,test);					/* sqrt(kN) 		*/
	}

}   /* end of discrim_roots */
 
 
/************************************************************************/
/*																		*/
/*			routine to select coefficients for the polynomials			*/
/*																		*/					   	
/************************************************************************/
 
select_coeffs(ktN)
int ktN[];
 
{   /* start of select_coeffs */
register int done,j;
int temp2,temp3,result1,result2,testv;
int bigtemp1[MPDIM], bigtemp4[MPDIM],result0[MPDIM];
double t1;

#define DEBUG_SELECT	0

if (TIME) t1 = (double)prectime();
	 				
testv = get_next_prime();		    	    							   	

/*	find a prime near temp1 which equals 3 mod 4 and which has		*/
/*	kN as a quadratic residue.                             	        */
/*  note: simple optimization; delete 1 mod 4 entries in table!     */
/*  only above 800 (say). Need ALL primes up to maxp                */
// Note: compare timing on this routine against old one.
 
done = 0;
while (done == 0)
	{
	while ((testv & 3) == 1) testv = get_next_prime();
	if (DEBUG_SELECT) (void) printf("testv = %d\n",testv);
	j = div_single(ktN, testv, bigtemp1);
	if (jacobi_symbol(j,testv) == 1)		
		{ 
		done = 1;								/* we found one		*/
		D = testv;
		A = D*D;
		}
	else testv = get_next_prime();
	if (testv > SQRT_RAD) return(-1);	/* testv gets too big; quit	*/	
	}

/*		set DINV = 1/(D) mod kN (we'll need it later)				*/

if (DEBUG_SELECT) (void) printf("D, A = %d %d\n", D, A);
OLD_D = testv;

single_modinv_dp(D,ktN,DINV);
if (DEBUG_SELECT) { (void) printf("DINV = "); write_num(DINV, numstr); }

if (DEBUG_SELECT)
	{
    mul_single(DINV, D, bigtemp1);
    divv(bigtemp1, ktN, bigtemp4, result0);
	printf("D * DINV = "); write_num(result0, numstr);
	}

/*		now select B and C subject to B^2 - AC = kN					*/

temp2 = (D+1) >> 2;
j = div_single(ktN, D, bigtemp1);
result1 = vsmall_power(j, temp2, D);

if (spower(ktN,temp2,D) != result1) printf("FAIL!\n"); 
 
if (DEBUG_SELECT) (void) printf("result1 = %d\n", result1);

/*		solve the equation 2*result1*X = (kN-result1^2)/D mod D		*/

temp3 = single_modinv((result1 << 1),D); 
temp2 = result1*result1; 

if (DEBUG_SELECT) (void) printf("result1^2, temp3 = %d %d\n",temp2,temp3);

SIZE(result0) = SINGLE;
LASTDIGIT(result0) = temp2;
subt(ktN, result0, bigtemp1);

j = div_single(bigtemp1,D,bigtemp4);		/* (kN-result1^2)/D	*/
if (DEBUG_SELECT) {(void) printf("bigtemp4 = "); write_num(bigtemp4, numstr);}

if (DEBUG_SELECT) (void) printf("j = %d\n",j);

mul_single(bigtemp4,temp3,result0);			/*	X				*/
result2 = div_single(result0,D,bigtemp1);	/* reduce mod D		*/

if (DEBUG)  (void) printf("X = %d\n", result2); 
 
/*			now B = X*D + result1								*/

B = (result2 * D + result1) % A;

if (DEBUG) (void) printf("B is: %d\n",B);  
 
/*		Find out where this polynomial becomes negative					*/

single_add(kN_sqrt,B,bigtemp4);					/* B[i] + SQRT(kN)		*/

(void) div_single(bigtemp4,A,result0);			/* divide by 2A			*/
crossover1 = -LASTDIGIT(result0);

single_subt(kN_sqrt,B,bigtemp4);				/* SQRT(kN) - B[i]		*/
(void) div_single(bigtemp4,A,result0);			/* divide by 2A			*/
crossover2 = LASTDIGIT(result0);

if (DEBUG_SELECT) (void) printf("cross points: %d %d\n",crossover1,crossover2);

if (TIME) 
	{
	t1 =  ( (double)prectime() - t1)/RATE;
	select_tot += t1;
//	(void) printf("select time = %g\n",t1);
	}

return(0);

}   /* end of select_coeffs */

/************************************************************************/
/*																		*/
/*		routine to get scaled logarithms of the primes					*/
/*																		*/
/************************************************************************/
 
void qs_get_logs(nfactors,factors,logp)
int nfactors;
int factors[];
char logp[];
 
{   /* start of get_logs */

register int i;
 
for (i=1; i<=NFACTORS; i++)
	{
	logp[i] = log_primes[indices[i]];
	}

}   /* end of get_logs */
 
/************************************************************************/
/*																		*/
/*		routine to get scaled log of test value for factorization		*/
/*																		*/
/*		the average value of q(x) is M*sqrt(kN)/2. compute log of		*/
/*		this value less large prime tolerance							*/
/*																		*/
/************************************************************************/
 
char get_mainlog(kN,tolerance,pmax)
double tolerance;
int pmax,
    kN[];
 
{   /* start of get_mainlog */

double t1;
register char t2;
 
/*				first compute log(kN) 									*/
 	   	
t1 = log((double)RADIX) + log((double)FIRST(kN));
 
/* 				now compute the test value								*/
 
t2 = (char) ((t1/2.0 - tolerance*log((double)pmax))/LOG) +
     (char) (log((double)MAX_SIEVE/2.0)/LOG) + 1;

return(t2);
 
}   /* end of get_mainlog */
 
/************************************************************************/
/*																		*/
/*		routine to save a factored vector								*/
/*																		*/
/************************************************************************/
 
void output(signbit,num,nfactors,list)
int num[],list[],signbit,nfactors;

 
{   /* start of output */
int i,count;

if (DEBUG1) 
	{
	(void) printf("outputting: %d %d\n",numfull,list[0]);
	(void) printf("lhs = "); write_num(num,numstr);
	(void) printf("rhs = %d ",signbit);
	for (i=1; i<=list[0]; i++)
		(void) printf("%d ",list[i]);
	(void) printf("\n");
	for (i=1; i<=list[0]; i++)
		(void) printf("%d ",factors[list[i]]);
	(void) printf("\n===============\n");
	}
mpcopy(num, lhs[numfull]);

count = list[0];
faclist[numfull][0] = count;
for (i=1; i<=count; i++) faclist[numfull][i] = list[i];

if (signbit == 1) setbit(matrix[numfull], 0);

for (i=1; i<=count; i++) setbit(matrix[numfull], list[i]); 

if (DEBUG1) dump_row(numfull);
if (DEBUG1) (void) printf("================================================================\n");
 
}   /* end of output */
 
/************************************************************************/
/*																		*/
/*			routine to compute the factor base							*/
/*																		*/
/************************************************************************/
 
void comp_fact(nfactors,factors,base)
int nfactors,factors[],base[];
 
{   /* start of comp_fact */
register int ii, next, prime;
int t, temp[MPDIM];
 
factors[0] = 1; indices[0] = -1;
factors[1] = 2; indices[1] = 0;
ii = 2;
next = 1;
while (ii <= nfactors)
	{
	prime = primes[next];
	t = div_single(base, prime, temp);
	if (jacobi_symbol(t,prime) == 1)
//	if (isresidue(base,prime))
		{
		factors[ii] = prime;  
		indices[ii] = next;
		ii++; 
		}
	next++;
	}
}   /* end of routine */

/************************************************************************/
/*																		*/
/*		routine to see if a prime is a quadratic residue of base		*/
/*																		*/
/************************************************************************/
 
isresidue(base,prime)
int base[],prime;
 
{   /* start of isresidue */
 /* hopeless slow;  reduce base mod p then use single precision quad recip.*/
register int k,test;
int temp[MPDIM];
 
k = (prime-1) >> 1;
if (k == 1)
	test = div_single(base,prime,temp);
else
	test = spower(base,k,prime);
 
if (test == 0 || test == 1) return(1);
else return(0);
 
}   /* end of isresidue */
 

/************************************************************************/
/*																		*/
/*		routine to back reference an index for exponent vector			*/
/*																		*/
/************************************************************************/
 
refer_index(value)
int value;
 
{   /* start of refer_index */
 
register int i;
 
for (i=NFACTORS; i>=1; i--)
	{
	if (value == factors[i]) 
		{
		return(i);
		}
	} 
}   /* end of refer index */ 


/************************************************************************/
/*																		*/
/*		Routine to scan sieve array; This returns an array holding		*/
/*		locations for possible successes.								*/
/*																		*/
/************************************************************************/
 
scan_it(success)
int success[];

{   /* start of scan_it */
int i,j,k,count,check;
int *ptr, addval, tval;

count = 0;
ptr = (int *) sieve;
addval = 128 - test_val;
tval = (addval << 24) + (addval << 16) + (addval << 8) + addval;

for (i=0; i<(MAX_SIEVE >> 2); i+=4)
   {
   check = (ptr[i] + tval) | (ptr[i+1] + tval) |
	   (ptr[i+2] + tval) | (ptr[i+3] + tval);
   if ((check & MASK)  != 0)
		{
		for (k=i; k<i+4; k++)
			{
			if (((ptr[k] + tval) & MASK) == 0) continue;
			else for (j=4*k; j<4*k+4; j++)
				{
				if (sieve[j] > test_val)
					{
					success[count] = j;
					count++;
					}
				}
			}
		}
   }

if (DEBUG) (void) printf("hit count = %d\n",count);
return(count);

}   /* end of scan_it */


/************************************************************************/
/*																		*/
/*		Reset start points when sieve reverses direction				*/
/*																		*/
/************************************************************************/

void reverse_startpts()

{   /* start of reverse_startpts */
int i;

for (i=1; i<=NFACTORS; i++)
   {
   if (factors[i] <= SMALLP)
		{
		roots1[i] = smallp[indices[i]+1] - roots1[i];
		roots2[i] = smallp[indices[i]+1] - roots2[i];
		}
   else 
		{
		roots1[i] = factors[i] - roots1[i];
		roots2[i] = factors[i] - roots2[i];
		}
	}

}   /* end of reverse_startpts */

/************************************************************************/
/*																		*/
/*		Test single precision number for primality						*/
/*																		*/
/************************************************************************/

prp_single(N)
int N;

{
int i,bound;

bound = (int)sqrt((double)N);

i = 1;
while (primes[i] <= bound)  
	{
	if (N % primes[i] == 0) return(1);
	i++;
	}
return(0);
}


/**************/
/* single_add */   // too slow;  replace with direct addition
/**************/

void single_add(num, x, result)
int num[], x, result[];

{
int temp[MPDIM];

SIZE(temp) = SINGLE;
LASTDIGIT(temp) = x;
add(num, temp, result);
}

/***************/
/* single_subt */
/***************/

void single_subt(num, x, result)
int num[], x, result[];

{
int temp[MPDIM];

SIZE(temp) = SINGLE;
LASTDIGIT(temp) = x;
subt(num, temp, result);
}


/************************************************************************/
/*																		*/
/*		Finish up														*/
/*																		*/
/************************************************************************/

finish(kN,fac1, fac2)
int kN[], *fac1, *fac2;
{
int solutions[MAXSOL][FBASE],numsol;


if (DEBUG1) {(void) printf("calling finish on: "); write_num(kN,numstr);}

reduce_matrix(solutions,&numsol);

if (DEBUG1) {(void) printf("after reduce kN = "); write_num(kN,numstr);}

return(compute_gcd(kN, solutions, numsol, fac1, fac2));

}


/************************************************************************/
/*																		*/
/*		solve matrix													*/
/*																		*/
/************************************************************************/

void reduce_matrix(solutions, count)
int solutions[MAXSOL][FBASE];
int *count;

{
int rightmost[96],ident[96][3];
int i,j,k,current_col,ident_size,current_word,current_bit;
int ncolumns,temp,nrows,tcount;
double t1;

if (TIME) t1 = (double)prectime();
nrows = numfull;
ncolumns = MATSIZE;						/* bit matrix: # of columns	*/
ident_size = MATSIZE;					/* size of identity matrix	*/
solution_count = 0;

if (DEBUG1) (void) printf("BEFORE reduce\n");
if (DEBUG1) dump_matrix(numfull);

for (i=0; i<nrows; i++)
	{
	ident[i][0] = 0;
	ident[i][1] = 0;
	ident[i][2] = 0;
	rightmost[i] = -1;		/* and the rightmost bit pointers	*/
	setbit(ident[i],i);
	}
 
current_col = NFACTORS;				/* search for rightmost 1	*/

while (current_col >= 0)
   {
   current_word = current_col >> 5;
   current_bit = current_col & 31;
   for (i=0; i<nrows; i++)
	{
	if (bittest(matrix[i]))
	   {
	   if (rightmost[i] > current_col) continue;
	   rightmost[i] = current_col;	/* set rightmost for this row	*/
	   for (k=i+1; k<nrows; k++)	/* rightmost 1 found		*/
			{
			if (bittest(matrix[k]))
				{
				if (rightmost[k] > current_col) continue;
				for (j=0; j<current_word+1; j++)  
					matrix[k][j] ^= matrix[i][j];

				for (j=0; j<ident_size; j++)
					ident[k][j] ^= ident[i][j];
				}
			}
	   break;
	   }
	}
   current_col--;
   }   /* end of outer loop */
if (DEBUG1) (void) printf("after reduce\n");
if (DEBUG1) dump_matrix(numfull);
 
/*	 now search for a row of zeros in matrix (solution)		*/

for (i=nrows-1; i>0; i--)
	{
	for (j=0; j<ncolumns; j++)
		if (matrix[i][j] != 0) goto endloop;
 
/*	found a row of zeros, extract solution from ident matrix	*/
 
	temp = 0;
	for (j=0; j<nrows; j++)
		{ 
		current_word = j >> 5;
		current_bit = j & 31;
		if (bittest(ident[i]))  temp++;
		}
	if (temp < 2) continue;
	solutions[solution_count][0] = temp;
	tcount= 1;
	for (j=0; j<nrows; j++)
		{
		current_word = j >> 5;
		current_bit = j & 31;
		if (bittest(ident[i]))
		    {
		 if (DEBUG1) { (void) printf("doing row %d\n",j); fflush(stdout);}
			solutions[solution_count][tcount] = j;
			tcount++;  
			if (DEBUG1) { (void) printf("tcount = %d\n",tcount); fflush(stdout);}
		    }
		}
	if (DEBUG1) (void) printf("solution = ");
	if (DEBUG1) for (j=0; j<=temp; j++) printf("%d ",solutions[solution_count][j]); 
	if (DEBUG1) (void) printf("\n");
	if (DEBUG1) check_solution(solutions, solution_count);
	solution_count++;
	if (solution_count > 20) break;
endloop:;
	}

*count = solution_count;
if (TIME) 
	{
	t1 =  ( (double)prectime() - t1)/RATE;
	(void) printf("solution time = %g\n",t1);
	(void) printf("final solution count = %d\n", *count);
	}
}


/************************************************************************/
/*																		*/
/*		compute final GCD												*/
/*																		*/
/************************************************************************/

compute_gcd(kN, solutions, count, fac1, fac2)
int solutions[MAXSOL][FBASE];
int count, *fac1, *fac2;
int kN[];
{
int left[MPDIM], right[MPDIM], temp[MPDIM], temp1[MPDIM];
int answer[MPDIM], index, index1, fcount;
int accum[FBASE];
int quo[MPDIM],rem[MPDIM];
int i,j,k,setsize;
double t1;


if (TIME) t1 = (double)prectime();

if (DEBUG1) (void) printf("computing GCD %d %d\n",fac1, fac2);

for (i=0; i<count; i++)
	{
	SIZE(left) = SINGLE;
	SIZE(right) = SINGLE;
	LASTDIGIT(left) = 1;
	LASTDIGIT(right) = 1;
//	left[2] = 0;
//	right[2] = 0;
	setsize = solutions[i][0];   
	if (DEBUG1) (void) printf("setsize = %d\n",setsize);
	for (j=1; j<=setsize; j++)
		{
		if (DEBUG1) {(void) printf("multiplying: "); write_num(lhs[solutions[i][j]], numstr);}
		mult(left, lhs[solutions[i][j]], temp); 
		if (DEBUG1) {(void) printf("done mult\n"); fflush(stdout);}
		if (DEBUG1) {(void) printf("kN = %d %d",fac1,fac2); write_num(kN,numstr);}
		divv(temp, kN, temp1, left);
		if (DEBUG1) {(void) printf("product: "); fflush(stdout); write_num(left, numstr);}
		}
	if (DEBUG1) {(void) printf("final LHS = %d %d",fac1,fac2); write_num(left, numstr); fflush(stdout);}
	for (j=0; j<=NFACTORS; j++) accum[j] = 0;
	for (j=1; j<=setsize; j++)
		{
		index = solutions[i][j];
		if (DEBUG1) (void) printf("index = %d\n",index);
		fcount = faclist[index][0];
		if (DEBUG1) (void) printf("#factors in row = %d\n",fcount);
		for (k=1; k<=fcount; k++)
			{
			index1 = faclist[index][k];
			  if (DEBUG1) (void) printf("factor index = %d\n",index1);
			accum[index1] ++;
			}
		}

if (DEBUG1) (void) printf("initial accum = %d %d\n",fac1, fac2);
if (DEBUG1) { for (j=0; j<=NFACTORS; j++) (void) printf("%d ",accum[j]); }
if (DEBUG1) (void) printf("\n");

	for (j=1; j<=NFACTORS; j++)
		{
//		if ((accum[j] & 1)  == 1) { printf("parity err at %d\n",j); write_num(kN,numstr); exit(0); }
		accum[j] = accum[j] >> 1;
		}

if (DEBUG1) (void) printf("accum = %d %d\n",fac1, fac2);
if (DEBUG1) { for (j=0; j<=NFACTORS; j++) printf("%d ",accum[j]);   }
if (DEBUG1) {(void) printf("\n"); fflush(stdout);}

	for (j=1; j<=NFACTORS; j++)
		{
		if (accum[j] == 0) continue;
		if (accum[j] == 1) mul_single(right, factors[j], temp);
		else if (accum[j] == 2) mul_single(right, factors[j]*factors[j], temp);
		else if (accum[j] == 3) mul_single(right, factors[j]*factors[j]*factors[j], temp);
		else
			{
			single_power(factors[j], accum[j], kN, answer);
			if (DEBUG1) {(void) printf("j: %d,  %d ^ %d mod kN = ", 
								j,factors[j], accum[j]); write_num(answer,numstr); fflush(stdout);}
			mult(right, answer, temp);
			}
		divv(temp, kN, temp1, right);
		}

if (DEBUG1)
	{
	(void) printf("checking:\n");
	mult(left,left,temp);
	divv(temp,kN,quo,rem);
	(void) printf("left^2 = "); write_num(rem,numstr);
	mult(right,right,temp);
	divv(temp,kN,quo,rem);
	(void) printf("right^2 = "); write_num(rem,numstr);
	(void) printf("DOING sol #%d %d %d\n",i,fac1, fac2); fflush(stdout);
	}
	add(left, right, temp);  
	gcd(temp, kN, answer);  
	if (DISPLAY_GCD) {(void) printf("GCD = "); write_num(answer,numstr); fflush(stdout);}
	if (mpcmp(kN, answer) == 0) continue;
	if (SIZE(answer) > SINGLE || LASTDIGIT(answer) != 1)  //  both need to be single prec!
		{ 
		if (SIZE(answer) > SINGLE) return(0);
		*fac1 = LASTDIGIT(answer);
		div_single(kN, *fac1, temp);
		if (SIZE(temp) > SINGLE) return(0);
		*fac2 = LASTDIGIT(temp);
		if (DEBUG1) (void) printf("qs returning 1\n");
		if (TIME) 
			{
			t1 =  ( (double)prectime() - t1)/RATE;
			(void) printf("GCD time = %g\n",t1);
			}
		return(1);
		}
	}
if (TIME) 
	{
	t1 =  ( (double)prectime() - t1)/RATE;
	(void) printf("GCD fail  time = %g\n",t1);
	}
return(-1);  // fail;  no non-trivial GCD appeared.  rare

}



/************************************************************************/
/*																		*/
/*	routine to compute b^e MOD m for integers b,e, and MP m				*/
/*																		*/
/************************************************************************/
 
void single_power(b,e,m,a)
int b,e,m[],a[];
 
{   /* start of single_power */

int junkc[MPDIM],junk4[MPDIM];
int bits[MPBITS+1],ebit;
register int i,j;

ebit = e;
i = 0;
SIZE(a) = SINGLE;
LASTDIGIT(a) = b;
for (j=1; j<=MPBITS; j++)
	{
	i++;
	bits[i] = ebit & 1;
	ebit = ebit >> 1;
	if (ebit == 0) { i--; break; }
	}
for (; i>0; i--)
	{
	square_dble_prec_asm(a,junkc);
	divv(junkc,m,junk4,a);
	if (bits[i] != 0)
		{
		mul_single(a,b,junkc);
		divv(junkc,m,junk4,a);
		}
	}
 
}   /* end of single_power */


/************************************************************************/
/*																		*/
/*	routine to initialize variables										*/
/*																		*/
/************************************************************************/

void init_vars()

{
int i;

for (i=0; i<FBASE+10; i++)
   {
   matrix[i][0] = 0;
   matrix[i][1] = 0;
   matrix[i][2] = 0;
   faclist[i][0] = 0;
   }
for (i=0; i<FBASE; i++) indices[i] = 0;

}


/************************************************************************/
/*																		*/
/*		routine to compute multiplier									*/
/*																		*/
/************************************************************************/
 
void get_mult(k,n)
int *k,n[];
 
{   /* start of get_mult */
 
static int test[] = { 1,3,5,7,11,13 };
 
int temp1[MPDIM];
register int i,mod4, s = 6;
int fbase[11];
double result,max;


mod4 = LASTDIGIT(n) & 3;			/* get n mod 4		*/
 
*k = 1;
max = 1.0; 
for (i=0; i<s; i++)
	{
	mul_single(n,test[i],temp1);
	if (SIZE(temp1) > DOUBLE) return;	// is this really needed? YES!
	comp_fact(10,fbase,temp1);
	knuth(test[i],10,fbase,&result);
	if (DEBUG1) printf("KNUTH value = %d %d %f\n",i, test[i], result);
	if (result > max)
		{
		max = result;
		*k = test[i];
		}
	}
 
}   /* end of get_mult */

/************************************************************************/
/*																		*/
/*		routine to evaluate a multiplier								*/
/*	uses a modified version of the Knuth-Schroeppel function			*/
/*	the modification is that the probabilities are different for		*/
/*	QS than they are for the continued fraction algorithm				*/
/*																		*/
/*	Q(x) is divisible by 2 iff kN = 1 mod 8 and is divisible			*/
/*	on average 2 times.													*/
/*	If p | multiplier then Q(x) is divisible by p with probability		*/
/*	1/p																	*/
/*	else Q(x) = 0 mod p with probability 2/p							*/
/*																		*/
/************************************************************************/
 
void knuth(mult,numfac,factors,sum)
int mult,numfac,factors[];
double *sum;
 
{   /* start of knuth */
double t1;
register int i;

#define	LOG2		.693147181  
*sum = 0.0;
 
for (i=1; i<=numfac; i++)
	{
	t1 = (double)factors[i];
	if (factors[i] == 2)		 	/* Q(x)=0 mod 2		*/
		   *sum += 2.0*LOG2;		/* iff kN = 1 mod 8	*/
	else if (mult % factors[i] == 0)
		*sum += 1.0/t1 *log(t1);
	else
		*sum += 2.0/t1 *log(t1);
	}

*sum = *sum - .5*log((double)mult);
 
}   /* end of knuth */

/************************************************************************/
/*                                                                      */
/*	Routine to compute Legendre symbol			                        */
/*                                                                      */
/************************************************************************/

jacobi_symbol(a,n)
int a,n;
{   /* start of jacobi_symbol */
int j, temp, res;

j = 1;

   
while (a != 0)
   {
   while ((a & 1) == 0)        // a is even
      {
      a >>= 1;
      res = n & 7;
      if (res == 3 || res == 5) j = -j;
      }
   temp = a;
   a = n;
   n = temp;

   if ((a & 3) ==3 && (n & 3) ==3) j = -j;

   a = a % n;
   }

if (n==1) return (j);
else return (0);
}   /* end of jacobi_symbol */	


/************************************************************************/
/*                                                                      */
/*	Routine to compute a^n mod p for small integers                     */
/*                                                                      */
/************************************************************************/

vsmall_power(base, exp, modulus)

{    /* start of vsmall_power */
int base1, i, j, bits[16];

base1 = base;
i = 0;
for (j=1; j<16; j++)
   {
   i++;
   bits[i] = exp & 1;
   exp >>= 1;
   if (exp == 0) {i--; break;}
   }

for (; i>0; i--)
   {
   base1 = (base1 * base1) % modulus;
   if (bits[i] != 0)  base1 = (base * base1) % modulus;
   }

return(base1);
}   /* start of vsmall_power */

/************************************************************************/
/*                                                                      */
/*  Macro to fetch next entry from primes[] table                       */
/*  Note:  prime_pointer is global                                      */
/*                                                                      */
/************************************************************************/

int __inline get_next_prime()

{   /* start of get_next prime */
int p;

p = primes[prime_pointer];
prime_pointer++;
return(p);
}


/*************************/
/* check matrix solution */
/*************************/

void check_solution(solutions,which)
int which;
int solutions[MAXSOL][FBASE];
{
int i,j, acc[FBASE];
int count,row;

for (i=0; i<=NFACTORS; i++) acc[i] = 0;
count = solutions[which][0];

for (i=1; i<=count; i++)
   {
   row = solutions[which][i];
   if (row < 10) printf("fr%d:  ",row);
   else printf("fr%d: ",row);
   for (j=1; j<=faclist[row][0]; j++) (void) printf("%d ",faclist[row][j]);
   (void) printf("\n");
   for (j=1; j<=faclist[row][0]; j++)
		{
		acc[faclist[row][j]]++;
		}
	}

(void) printf("checking solution %d\n",which);
for (i=1; i<=NFACTORS; i++) if ((acc[i] & 1) == 1) (void) printf("parity bad at %d: %d\n",i,acc[i]);
for (i=1; i<=NFACTORS; i++) (void) printf("%d ", acc[i]);
(void) printf("\n");
}

/*******************/
/* dump the matrix */
/*******************/

dump_matrix(count)
int count;
{
int i,j;
int current_word, current_bit;

(void) printf("dumping matrix; numfull = %d\n",numfull);
for (i=0; i<count; i++)
	{
	(void) printf("%d : ",i);
	for (j=0; j<=NFACTORS; j++)
		{
		current_word = j >> 5;
		current_bit = j & 31;
		(void) printf("%d", bittest(matrix[i]));
		}
	printf(" -1 \n");
	}

fflush(stdout);
return(1);
}


/*******************/
/* dump a mat row  */
/*******************/

dump_row(row)
int row;
{
int j;
int current_word, current_bit;

if (row < 10) (void) printf("r%d:  ",row);
else (void) printf("r%d: ",row);
for (j=0; j<=NFACTORS; j++)
		{
		current_word = j >> 5;
		current_bit = j & 31;
		(void) printf("%d", bittest(matrix[row]));
		}
(void) printf("\n");


fflush(stdout);

return(1);
}

/************/
/*  dump id	*/
/************/

dump_id(ident,count)
int count;
int ident[96][3];
{
int i,j;
int current_word, current_bit;

(void) printf("dumping id\n");
(void) printf("numfull = %d\n",numfull);
for (i=0; i<count; i++)
	{
	(void) printf("%d: ",i);
	for (j=0; j<=NFACTORS; j++)
		{
		current_word = j >> 5;
		current_bit = j & 31;
		(void) printf("%d", bittest(ident[i]));
		}
	(void) printf(" -1 \n");
	}

fflush(stdout);
return(1);
}
