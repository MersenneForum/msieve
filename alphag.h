FILE *input;						/* data input					*/
FILE *rest;							/* holds restart results		*/
FILE *fbase;						/* holds factor base etc.		*/
FILE *output;						/* holds a,b data pairs			*/
FILE *count;						/* counts results				*/
FILE *outfile;						/* holds output data.			*/
FILE *rejfile;						/* holds rejected q's			*/
 
int  int_fbase[FBASE],				/* integer factor base			*/
     int_squares[8000],				/* and their small squares		*/
     alg_squares[8000],
     alg_fbase[FBASE],				/* algebraic factor base		*/
     qqueue[55],					/* queue for squfof routine		*/
     qpoint,						/* and pointer for queue		*/
     critical,						/* critical size for single p	*/
     alg_pmax,						/* largest prime on alg side	*/
     int_pmax;						/* largest prime on int side	*/
 						   	
int  N[MPDIM],						/* number to be factored		*/
     b,								/* coefficient of linear form	*/
     alg_roots[FBASE],				/* roots of f(x) = 0 mod p		*/
     int_roots[FBASE],				/* roots of a + bM = 0 mod p	*/
     int_startpts[NSIEVEP],			/* starting locations for sieve	*/
     alg_startpts[NSIEVEP],			/* ""  on alg side				*/
	 int_reset_start[NSIEVEP],		/* saved for resieving			*/
	 alg_reset_start[NSIEVEP],		/* i.e. factor by resieve		*/
     int_smallroot[NUMSMALLP],		/* small prime power roots		*/
     alg_smallroot[2*NUMSMALLP],	/* of a +bM						*/
     int_faclist[TSIZE][12],		/* integer factors: resieved	*/
	 alg_faclist[TSIZE][12],		/* algebraic factors: resieved	*/
     normpart[MPDIM],				/* part of norm : c*b^d			*/
     poly_coef[MAX_DEGREE+1],		/* polynomial coefficients		*/
	 num_alg_success,				/* #algebraic successes			*/
	 num_total_success,				/* #successes/both sides		*/
	 int_pmax_squared[MPDIM],		/* largest prime^2				*/
	 alg_pmax_squared[MPDIM],		/* largest prime^2				*/
     int_polyc[MAX_DEGREE+1][MPDIM];/* lhs polynomial coeffs		*/

short
	 success[SUCCESS_SIZE][2];		/* success locations in sieve	*/

int  v1[2],							/* q reduced lattice; 1st vector*/
     v2[2],							/* 2nd vector					*/	
     special_q,						/* current value of special q	*/
     special_q_root,				/* and its root					*/
	 int_count,						/* count of total saved primes	*/
	 alg_count;						/* by resieving					*/

int  special_q_smalliroot[NUMSMALLP],  /* special q start points	*/
     special_q_smallaroot[NUMSMALLP],  /* for 1st column			*/
     special_q_int_roots[NSIEVEP],
     special_q_alg_roots[NSIEVEP];

int	proj_int_primes[NUM_PROJ],		/* holds proj prime indices		*/
	proj_alg_primes[NUM_PROJ];

char proj_int_logs[NUM_PROJ],		/* and their logarithms			*/
	 proj_alg_logs[NUM_PROJ];

int int_root_ok[NUMSMALLP],			/* indicates when no root exists*/
    alg_root_ok[NUMSMALLP];			/* for powers of small primes	*/

int alg_prime_power_root[NUMSMALLP],	/* and the actual roots		*/
	power_root_count;

int	int_v1[FBASE][2],				/* vectors for sub_p reduced	*/
    int_v2[FBASE][2],				/* lattices						*/
    alg_v1[FBASE][2],				/* LOTS of memory!				*/
	alg_v2[FBASE][2];
 
int  numtot,						/* total factorizations			*/
     numpf,							/* partial factorizations		*/
     numfp,
     numpp,
     numff,
     numppf,
     numppq,
	 numqpp,
	 numfpp,
	 numpppp,
     fdegree,						/* degree of the field			*/
     int_degree,					/* degree of LHS polynomial		*/
     poly_constant,					/* constant for f(x)			*/
     field_base,					/* b of N = b^k + s				*/
     lhs_degree,					/* d of b^d = alpha mod N		*/
     starting_q,					/* starting q index for this run*/
     ending_q,						/* ending q index for this run	*/
     starting_a,					/* starting a for sieve interval*/
     genpoly,						/* flag for general polynomial	*/
     restart;						/* restart option				*/
 
int size_int_base,					/* size of factor base: lhs		*/
    size_alg_base,					/* size of factor base: rhs		*/
    lhs_value[MPDIM],				/* base^exponent for lhs		*/
	int_cutoff,						/* minimal indices for resieve	*/
	alg_cutoff,
	int_8hits,
	int_4hits,						/* factor base partitions		*/
	int_3hits,						/* according to #hits in sieve	*/
	int_2hits,						/* interval						*/
	int_1hits,
	int_linesieve,
	alg_8hits,
	alg_4hits,
	alg_3hits,
	alg_2hits,
	alg_1hits,
	alg_linesieve,
	BPRIME_LIMIT,					/* largest big prime; read from file	*/
    SPLIT_LIMIT,					/* this^2. size limit for cofactors		*/
    nsieve;							/* number of intervals to sieve			*/

double tolerance,					/* tol for large prime: pmax^tolerance	*/
       log_b,						/* log base 3 of current value of b		*/
       log_pmax,					/* log(pmax)/LOG3 * tolerance			*/
       log_field_base,				/* log of the polynomial const (base 3)	*/
       log_PART1,					/* log of the lhs for aM2+bM1, logM1	*/
	   log_PART2;					/* log of the lhs for aM2+bM1, logM2	*/

double emin,						/* coordinates and slopes of lines that	*/
	   emax,						/* define vector sieving boundaries		*/
	   lower_midpt,
	   upper_midpt,
	   lm1, lc1,
	   lm2, lc2,
	   um1, uc1,
	   um2, uc2;
 
unsigned char 
	 test_int_val,  				/* test value for full factor:lhs		*/
     test_alg_val,					/* test value for full factor:rhs		*/
     log_special_q,					/* log(special_q)						*/
     numstr[MPCDIM],				/* string for multi-precise  output		*/
     *sieve,						/* byte array for line sieving			*/
     int_logp[FBASE],				/* scaled logs of factor base: int side	*/
     alg_logp[FBASE],				/* scaled logs of factor base: alg side	*/
     s3465[3465],					/* table s3465 of QR's for squfof		*/
     sieve_array[NCOLUMNS+1][MAX_SIEVE+1];	/* global sieve array			*/

static int smallp[] = { 256,243,625,343,121,169,289,361,529,841,961 };

static int alg_smallp[] = {0,0,256,243,0,625,0,343,0,0,0,121,0,169,0,0,0,
						   289,0,361,0,0,0,529,0,0,0,0,0,841,0,961};

static int pow2[8] = {1, 2, 4, 8, 16, 32, 64, 128};

char label[256];							/* identifying label		*/
unsigned char intlogsmall[NUMSMALLP];		/* logs of small primes,	*/
											/* not their powers.		*/

double degree,								/* (double)field_degree		*/
	fconst,									/* polynomial constant		*/
	pcoef0,									/* coeffs if not pure field	*/
	polyc_float[MAX_DEGREE+1],
	pcoefd,
	log_pmax_rhs,							/* scaled logs of pmax for	*/
	log_pmax_lhs,							/* left and right sides		*/
	lhs_tol;								/* left tolerance (rhs *.95)*/

double LOGB;								/* logarithm base			*/
 
int	field_flag,								/* indicates if 'pure' field*/
    active_project_id;
 
int	int_bad,								/* int side squfof too big	*/
	alg_bad,								/* alg side squfof too big	*/
	int_prime,								/* int side prime too big	*/
	alg_prime,								/* alg side prime too big	*/
	int_cofactor,							/* comp cofactor too big	*/
	alg_cofactor,							/* ditto for rhs			*/
	total_hits,								/* count number of hits		*/
	squfof_fail,							/* squfof statistics		*/
	overflow_count,
	squfof_calls,
    badsplit,
	rho_calls,
	relation_count,
	not_coprime,
	total_q_successes;
 
int	int_factor_list[128],					/* stores primes 			*/
	alg_factor_list[128],
	bpsquared[MPDIM],						/* BPRIME_LIMIT^2			*/
	coeffmodp[NUMSMALLP];					/* proj root/b coeff mod p	*/

int	dump_cutoff;							/* ouput bound large primes	*/
 
double vector_setup_time,					/* timing totals			*/
	   int_linesieve_time,
	   alg_linesieve_time,
	   even_int_linesieve_time,
	   even_alg_linesieve_time,
       special_q_time,
	   int_vector_time,
	   alg_vector_time,
	   int_resieve_time,
	   alg_resieve_time,
	   trial_int_time,
	   trial_alg_time,
	   find_startpt_time,
	   alg_scan_time,
	   squfof_time,
	   latred_time,
	   find_success_time,
	   global_time,
	   total_ptime,
	   invtime,
	   total_q_time,
	   clockrate;

int alg_upper_area[500000],					/* for debugging only		*/
	alg_lower_area[500000],
	int_upper_area[500000],
	int_lower_area[500000];	   