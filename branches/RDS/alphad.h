
#define MPDIM			150				/* size of MP arrays		*/
#define MPCDIM			300				/* char size of MP number	*/
#define	MAX_DEGREE		7				/* maximum degree of field	*/
#define RADIX			1073741824		/* MP number base			*/
#define	HALF_RAD		46340			/* [SQRT(RADIX)]			*/
#define FBASE			1850001			/* Max Factor base size		*/
#define	BIGNUM			127				/* largest number for 1 byte*/
#define	MASK			0x80808080		/* top bit in each byte lit	*/
#define	RESET			0x81818181		/* -127 put in each byte	*/
#define SMALLP			29				/* Small prime variation lim*/
#define NUMSMALLP		10				/* Number of small primes	*/
#define AREA			7
#define MAX_SIEVE		(AREA*1024)		/* maximum sieve length		*/
#define	NCOLUMNS		(AREA*1024)		/* # of columns to sieve	*/
#define	LOG2			.69314718		/* Natural Log of 2.0		*/
#define	LOG3			1.0986122		/* Natural Log of 3.0		*/
#define SIEVE_INIT		16843009		/* 1 in each byte of a word	*/
#define	STORE_LIMIT		MAX_SIEVE		/* storage limits for factor*/
#define	PRIME_LIMIT		11				/* by resieve				*/
#define SUCCESS_SIZE	(AREA/2*MAX_SIEVE)	
#define	FLIST_SIZE		MAX_SIEVE
#define	BIG				MAX_SIEVE+1		/* do not sieve proj root	*/
#define	NUM_PROJ		16				/* # projective primes		*/
#define	NSIEVEP			8192			/* # line sieve primes		*/
#define	GOOD			255				/* marker for good hits		*/
#define	MAX_LINES		2000			/* max #lines in range file	*/
#define	STAT_FREQ		1000			/* how often display stats	*/
#define	CLOCKRATE		1500000.		/* 1.5 GHz --> millisec		*/
#define	MAX_RED_BITS	18				/* limit on V1, V2			*/

#define	COLUMN_STORE	1				/* Storage locations in		*/
#define	LOC_STORE		0				/* success array			*/

#define	DMAX_SIEVE		(double)MAX_SIEVE
#define	DCOLUMNS		(double)NCOLUMNS				
		
#define	TRUE			1
#define	FALSE			0

#define	FUZZY			3				/* linesieve: FUZZY*MAX_SIEVE*/

#define TSIZE			MAX_SIEVE*80	/* storage for resieve		 */
#define	SHIFT			6				/* for hashing;				 */

/*		Algorithm Parameters									     */

#define	SLEEPTIME	3000			/* time to wait: fopens millisec	*/

#define	ALG_REFINE_FREQ	100			/* how often refine alg threshhold	*/
#define INT_REFINE_FREQ	100			/* "   "   int threshhold			*/
 
 
/*		Multi-precise number Macros									*/

#define	SIZE(x)		x[0]			/* Size of x in base RADIX		*/
#define FIRST(x)	x[SIZE(x)-1]	/* First digit of x				*/
#define	SECOND(n)	n[SIZE(n)-2]	/*  2nd limb					*/
#define LAST(x)		x[1]			/* Last digit of x				*/
#define SINGLE		2				/* single precision size of MP	*/
#define	DOUBLE		3				/* size of dble precision int	*/

#define COPY(a,b)	for (il=0; il<SIZE(a); il++) b[il] = a[il]
#define COPY1(a,b,ind)	for (ind=0; ind<SIZE(a); ind++) b[ind] = a[ind]

#define	PUREFIELD	1				/* f(x) = x^d - k				*/
#define	GENFIELD	1
#define	SIGN(x)		(x < 0 ? -1 : 1)
#define	LOGR		log((double)RADIX)


/*				Bit set/test Macros									*/

#define	bitset(i)	bitarray[(i) >> 3] += pow2[(i) & 7]
#define bittest(i)	bitarray[(i) >> 3] & pow2[(i) & 7]
 
/*	Defines for different output cases				*/
 
#define	IL1	1
#define	IL2	2
#define AL1	4
#define	AL2	8

#define	FF	0
#define	FP	4
#define	PF	1
#define	PPF	3
#define	PP	5
#define	PPQ	7
#define	QPP	13
#define FPP	12
#define	PPPP 15

typedef	struct
	{
	int list[10];
	int count;
	int alg_rem, alg_rem2;
	int a,b;
	} factorization;

/************************************************************************/
/*																		*/
/*	Floor & ceiling macros; surprisingly these are 9x faster than		*/
/*	ceil() and floor()  on the Pentium!	!!!!!							*/
/*	note: wrong for exact integers										*/
/*																		*/
/************************************************************************/

#define	iceil(a)  ((a) <= 0.0 ? (int)(a) : (int)(a) + 1)
#define	ifloor(a) ((a) >= 0.0 ? (int)(a) : (int)(a) - 1)

/************************************************************************/
/*																		*/
/*	Macro to compute lower bound on f as a function of e				*/
/*																		*/
/************************************************************************/

//#define lower_bound(e) (e < lower_midpt ? iceil(lm1*e + lc1) : iceil(lm2*e + lc2))
#define lower1(e) (e < lower_midpt ? iceil(lm1*e + lc1) : iceil(lm2*e + lc2))
#define lower_bound(e) (e < lower_midpt ? round_double(lm1*e + lc1) : round_double(lm2*e + lc2))
#define lower_boundf(e) (e < lower_midpt ? ceil(lm1*e + lc1) : ceil(lm2*e + lc2))
#define lower1_bound(e) (ceil(lm1*e + lc1))
#define	lower2_bound(e) (ceil(lm2*e + lc2))

/************************************************************************/
/*																		*/
/*	Macro  to compute upper bound on f as a function of e.				*/
/*																		*/
/************************************************************************/

#define upper_bound(e) (e < upper_midpt ? ifloor(um1*e + uc1) : ifloor(um2*e + uc2))
#define	upper1_bound(e)	(floor(um1*e + uc1))
#define	upper2_bound(e)	(floor(um2*e + uc2))

//#define	min(a,b) ((a) < (b) ? (a) : (b))   already in stdlib.h
//#define   max(a,b) ((a) > (b) ? (a) : (b))

/************************************************************************/
/*																		*/
/*	Defines needed for sieve by bucket                                  */
/*																		*/
/************************************************************************/

#define SCALE					180000
#define MEM_FOR_POINT_UPDATES   (400*1024*1024)      
#define SORT_CHUNK_SIZE         (255*1024)      
#define POINTS_PER_SORT_CHUNK   (SORT_CHUNK_SIZE / sizeof (struct point_data))
#define NUM_SORT_CHUNKS         (MEM_FOR_POINT_UPDATES / SORT_CHUNK_SIZE + 1)
#define LAST_POINT_PTR			(((struct point_data *)((char *)point_array + MEM_FOR_POINT_UPDATES))-1)

// By combining multiple d values into one sort key, in do_delayed_points
// we read more cache lines in each sort chunk.  This reduces TLB misses
// and gives the hardware prefetcher a chance to do us some good.  It also
// saves some memory in sorting each chunk.  The downside is in
// do_delayed_points we'll use more L2 cache to hold several sieve
// lines.

#define SORT_KEY_GROUPING   5   /* sort key has 8 d values */

#define NUM_SORT_KEYS      ((MAX_SIEVE>>SORT_KEY_GROUPING)+1)
#define MAKE_SORT_KEY(d)   ((d)>>SORT_KEY_GROUPING)
