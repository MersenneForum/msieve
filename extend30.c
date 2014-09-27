/************************************************************************/
/*                                                                      */
/*               L A T T I C E   S I E V E                              */
/*                                                                      */
/*   Program to factor numbers using the Number Field Sieve             */
/*   algorithm with 2 large primes on each side (4 in all)              */
/*                                                                      */
/*                                                                      */
/*   This code handles the sieving and factoring of residues.           */
/*   The algorithm has a + bM = phi(a + b alpha)   mod N.               */
/*   where alpha is a complex root of the extension polynomial, and M   */
/*   is a root of that polynomial mod N.                                */
/*   We refer to the left hand side as the integer side, and the        */
/*   right as the algebraic side.                                       */
/*                                                                      */
/*   NOTE: code intended to be read by editor with tabs = 4 spaces      */
/*   This is the WINDOZE version!  [not Unix!]                          */
/*                                                                      */
/*   Bob Silverman  V1: 10/04;                                          */
/*   Additions by George Woltman V3: 7/05                               */
/*   Extensive optimizations 7/05, 8/05                                 */
/*                                                                      */
/************************************************************************/

/************************************************************************/
/*   Note:   special code for projective roots; We deal with            */
/* them as follows:  In find_start_points if a denominator is 0, then   */
/* we have a projective root. Keep a list of these primes. This list    */
/* will be small since sum(1/p) ~ loglog(maxprime). For each proj prime */
/* set its startpt to beyond the sieve interval. Set its root to 0 so it*/
/* does not change as j <-j+1. Thus it will never get sieved. In sieve  */
/* init, add the value of log(proj prime) to the init value when prime  */
/* divides j because it divides for all i for this j. Then, in trial_div*/
/* deal with this prime separately.                                     */
/************************************************************************/

/************************************************************************/
/*   Note:  we sieve a region from -MAX_SIEVE < a < MAX_SIEVE and       */
/*   1 <= b <= NCOLUMNS;  first we do the upper half plane (a > 0)      */
/*   then the lower, keeping track of the sign.    In the q-lattice     */
/*  for primes that have at least one hit [p < MAX_SIEVE], we line      */
/*  sieve; for bigger primes, we sieve by vectors                       */
/*  We refer to the GLOBAL lattice as the (a,b) lattice = (cV1 + dV2)   */
/*  and the q-lattice as the (c,d) lattice.  The p sub-lattice of the   */
/*   q-lattice is the (e,f) lattice.                                    */
/************************************************************************/


/************************************************************************/
/*   Note:  for the smallest primes, we sieve only with respect to their*/
/*   powers. For the rational side, this is not a problem. However, even*/
/*   if f(x) has a root mod p,  it might not have one mod p^d.  We only */
/*  sieve those which have roots mod p^d  and skip the rest             */
/************************************************************************/


/************************************************************************/
/*  SCREWBALL NOTE:                                                     */
/*   The procedure that determines clockrate [get ticks, wait 3 sec,then*/
/*  compute difference DOES NOT WORK on my laptop.  Each time it returns*/
/*  a different !^&@#&@  number!!  Apparently, the laptop recognizes    */
/*   the Sleep() command and puts the processor into power save mode of */
/*  some kind.  However,  if a second program is running that keeps the */
/*  processor active while this siever is running, then the clockrate   */
/*  computation DOES work and consistently returns the correct clockrate*/
/*  This is too weird for words.  When nothing else is running, my      */
/*  laptop becomes a non-deterministic machine with respect to measuring*/
/*  clockrate.                                                          */
/************************************************************************/


/************************************************************************/
/*   Note:  I replaced ceil() with a macro. This macro is 9x faster, BUT*/
/*  if called with an exact integer as its argument it returns n+1, not */
/*  n.  However, all this means for the sieving is that if a feasible   */
/*   lattice point falls EXACTLY on the bottom side of a parallelogram, */
/*  it will not be hit. Such points are VERY rare. The speed improvement*/
/*   makes missing such points worthwhile. Only points in the strict    */
/*  interior get hit.                                                   */
/************************************************************************/

/************************************************************************/
/*   Note: I have added, as input, a bound on the size of cofactors that*/
/*  will be attempted to be split by QS/SQUFOF. While the LP bound can  */
/*  be as large as 2^30, it does not pay to try to split cofactors too  */
/*  close to 2^60 because the vast majority will not split into 2 primes*/
/*  less than the LP bound. Most will split into a small prime and one  */
/*  that is too big. This saves time by not attempting to split         */
/*  cofactors that are unlikely to yield two primes less than the bound */
/************************************************************************/

/************************************************************************/
/*                                                                      */
/*   Here is the overall code:  read data, looping over special q first */
/*   find main lattice. Then find start points and sub-lattices for     */
/*  primes. Init the sieve array, vector sieve big algebraic primes,    */
/*  line sieve small ones. Mark potential successes in sieve array and  */
/*   set up sieve array for rational sieve. Vector sieve big rational   */
/*  primes, then line sieve the small ones. Scan sieve array. Mark      */
/*  potential candidates and store them. Then factor-by-resieve the big */
/*   primes.If a prime hits a candidate point, save it and the pt coords*/
/*  in a linear list. The sieve is again by vectors. Then call trialdiv */
/*  on each candidate. Extract the factor-by-resieve primes by scanning */
/*   the save linear list. This should not be large. Divide out these   */
/*  primes, then finish trial factor by trial division on the small     */
/*  primes. Note:  indexing into the list is via hash                   */
/*                                                                      */
/************************************************************************/ 

/*  Remaining to do:
   (1) Possible Optimizizations

Divide out vector primes in trial_alg  first?  This will only
have a very minor effect, if it is noticeable at all....

Finish red-doing latred


======================================================
   To reduce this, I note the following:

      The biggest impact would come from reducing the vector sieve time.
      George Woltman added a bucket sieve approach that reduced this quite
      a bit already.  MANY Kudos and Thanks to George.

*/

#include <stdio.h>
#include <math.h>
#include <sys/types.h>
#include <stdlib.h>
#include <share.h>
#include <string.h>
#include <wtypes.h>
#include <winbase.h>
#include <windef.h>
#include <time.h>
static int gwcount = 0;    //GW
static int gwcount2 = 0;   //GW

/*      #include <sys/times.h>  (Unix only)  or calls to get_rusage           */

#define   FACTOR_BASE_ASCII    "ASCII_BASE"
#define   FACTOR_BASE_BIN      "NFS_BASE"
#define   NFS_DATA             "nfsl.dat"
#define   STOP_FILE            "stopfile"
#define   HOME_DIR              ""

/*         Debug flags and statistics keeping flags                           */

#define   DEBUG              0            /* general debugging                */
#define   DUMPNORM           0            /* display polynomial norms         */
#define   DEBUG_ALG_NORM     0
#define   NOTIFY_HITS        0            /* show when hits occur             */
#define   FALSESHOW          0            /* show when false hits occur       */
#define   SHOW_START         0            /* show sieve init values           */
#define   ALG_DEBUG          0            /* debug statements in trial div    */
#define   ALG_SMALL          0
#define   INT_DEBUG          0            /* debug statements in trial div    */
#define   INT_SMALL          0
#define   SCAN_DEBUG         0
#define   CHECK_INIT         0
#define   FALSECOUNT         0            /* count false hits?                */

#define   Q_OVERFLOW         0            /* count squfof queue failures      */
#define   RHO_TALLY          0            /* count calls to Pollard Rho       */
#define   SQUFOF             0            /* count calls to squfof            */
#define   SQUFOF_FAIL        0            /* count when squfof fails          */
#define   SPLITCHECK         0            /* if 1 small 1 big <- squfof       */

#define   HIT_COUNT          0            /* count sieving successes          */
#define   DEBUG_HITS         0            /* show sieve partitions            */
#define   DEBUG_FIND         0            /* show start point computations    */
#define   DEBUG_SCANNER      0            /* various debug flags.....         */
#define   DEBUG_VECTOR       0
#define   DEBUG_VECTOR_PT    0
#define   DEBUG_SUCCESS      0
#define   DEBUG_SCAN_VALS    0
#define   DEBUG_RESIEVE      0
#define   DEBUG_REF_ALG      0
#define   DEBUG_REF_INT      0
#define   DEBUG_INT_SUCCESS  0         
#define   DEBUG_INT_VECTOR   0
#define   DEBUG_INT_VECTOR_R 0
#define   DEBUG_ALG_VECTOR   0
#define   DEBUG_ALG_VECTOR_R 0
#define   DEBUG_OUTPUT       0
#define   DEBUG_SIEVE        0
#define   DEBUG_UPDATE       0
#define   DEBUG_NORM         0
#define   DEBUG_Q            0
#define   DEBUG_QS           0

#define   DISPLAY_REDUCTIONS 0
#define   CHECK_LATTICE      0
#define   SHOW_F             0
#define   SHOW_INT_REGION    0
#define   SHOW_INT_REGION_R  0
#define   SHOW_INT_LATTICE   0
#define   SHOW_INT_LATTICE_R 0
#define   SHOW_PREP          0
#define   SHOW_ALG_REGION    0
#define   SHOW_ALG_REGION_R  0
#define   SHOW_ALG_LATTICE   0
#define   SHOW_ALG_LATTICE_R 0
#define   SHOW_LARGE_PRIMES  0
#define   CHECK_BOUNDS       0
#define   CHECK_BOUNDS_R     0
#define   TIME_STATS         0
#define   GLOBAL_TIME_STATS  0
#define   SHOW_MAIN          1
#define   SHOW_SUCCESS_COUNT 0
#define   ALG_HASH_DEBUG     0
#define   COMPUTE_AREA       0
#define   CHECK_HASH         0
#define   SPECIAL_Q_TIME     1

#define   LAPTOP             0            /* if running on a LAPTOP     */
#define   LAPTOP_RATE        1600000.0    /* 1.6 GHz currently          */

/************************************************************************/
/*                                                                      */
/*      G L O B A L      V A R I A B L E S                              */
/*                                                                      */
/************************************************************************/

#include "alphad.h"         /* defined parameters for this program      */ 
#include "mpdim.h"
#include "alphag.h"
#include "functions.h"

int big_special_q[100000], big_roots[100000];

/************************************************************************/
/*                                                                      */
/*      F U N C T I O N   D E C L A R A T I O N S                       */
/* inline functions must be declared here; they can not be declared in  */
/* functions.h  for scoping reasons                                     */
/*                                                                      */
/************************************************************************/



/************************************************************************/
/*                                                                      */
/*                  M A I N                                             */ 
/*                                                                      */
/************************************************************************/
 
__cdecl main()
 
{   /* start of main */
register int i,q;
int line_num,totq,q1,r1,offset;
double seconds;
char namestring[128], outname[128], restname[128], rejname[128];
char hostname[256];
WORD wVersionRequested;             /* yech!  uSoft hack!               */
WSADATA wsaData;
FILE *linroots;

if ((MAX_SIEVE & 127) != 0) 
   {
   (void) printf("BAD VALUE FOR MAX_SIEVE\n");
   exit(0);
   }

(void) printf("Siever built on %s %s\n", __DATE__, __TIME__);

seconds = get_time();
Sleep(3000);                       /* this is 3 seconds                */
clockrate = (get_time() - seconds)/3000.0;
if (LAPTOP)  clockrate = LAPTOP_RATE;
(void) printf("Clock rate = %lf\n",clockrate);


int_bad =   0;
alg_bad =   0;
int_prime = 0;
alg_prime = 0;
int_cofactor =   0;
alg_cofactor =   0;
total_hits =     0;
squfof_fail =    0;
overflow_count = 0;
rho_calls =      0;
squfof_calls =   0;
totq =  0;

/*   WINDOZE requires this crazy call  WSAStartup before gethostname  */
/*   can be called;  hence the crazy %?*#%^!! code below              */

wVersionRequested = MAKEWORD(1, 1); 
if (WSAStartup(wVersionRequested, &wsaData) != 0) return(0);

gethostname(hostname,256);

input = fopen(NFS_DATA,"r");            /* open input data file       */
 
get_file(N,input);                      /* read the number            */
get_file(lhs_value,input);              /* Root of LHS                */
(void) fscanf(input,"%s",label);        /* identifying label          */
(void) fscanf(input,"%d",&dump_cutoff); /* prime output bound         */
(void) fscanf(input,"%lf",&LOGB);       /* logarithm base             */
(void) fscanf(input,"%d",&field_flag);  /* pure field? (x^d - k)      */
 
(void) printf("We try to factor: ");
write_num(N,numstr);
(void) fflush(stdout);

(void) strcat(label,hostname);      /* add hostname ID to label       */
 
if (!get_ascii_fbase())             /* read in ascii factor base      */
   get_fbase();                     /* if not ascii, then binary      */
 
if (field_flag != PUREFIELD)        /* get coeffs if field not 'pure' */
   {                                /* these are the coeffs of the    */
   fdegree = field_flag;            /* field extension polynomial     */
   for (i=0; i<=fdegree; i++)
       {
     (void) fscanf(input,"%d ",&poly_coef[i]);
     polyc_float[i] = (double)poly_coef[i];
     }
   pcoef0 = polyc_float[0];
   pcoefd = polyc_float[fdegree];
   }

(void) fscanf(input,"%d",&int_degree);
for (i=0; i<=int_degree; i++)      /* coeffs of norm of linear poly   */
   get_file(int_polyc[i],input);

(void) fscanf(input,"%d",&BPRIME_LIMIT);   /* limit on large primes   */
(void) fscanf(input,"%d",&SPLIT_LIMIT);    /* cofactor size limit     */
SIZE(bpsquared) = SINGLE;
FIRST(bpsquared) = SPLIT_LIMIT;
mul_single(bpsquared,SPLIT_LIMIT,bpsquared);


/*   Partition factor base according to #hits the primes take         */
/*   in a sub-interval; this is for loop unrolling during sieving     */

partition_loops();

if (DEBUG) 
   {
   (void) printf("int_pmax, alg_pmax = %d %d\n",int_pmax,alg_pmax);
   (void) fflush(stdout); 
   }

/*         compute squares of primes; useful for divide routine       */
 
i = 1;
while (i < int_linesieve)
   {
   int_squares[i] = int_fbase[i]*int_fbase[i];
   i++;
   }

i = 1;
while (i < alg_linesieve)
   {
   alg_squares[i] = alg_fbase[i]*alg_fbase[i]; 
   i++;
   }

/*      compute the scaled logarithms of the two factor bases         */

get_logs(int_fbase,size_int_base,int_logp,1);
get_logs(alg_fbase,size_alg_base,alg_logp,2);
 
for (i=0; i<NUMSMALLP; i++)            /* logs of small prime powers  */
   {
   intlogsmall[i] = (unsigned char)(log((double)smallp[i])/LOGB + .5);
   }
 
/*      compute  some constants that will help the computation of     */
/*      the sieve threshhold values.                                  */        

log_pmax_rhs = ALG_TOL * log((double)BPRIME_LIMIT)/LOGB;
log_pmax_lhs = INT_TOL * log((double)BPRIME_LIMIT)/LOGB;
log_PART1 = ((double)(SIZE(int_polyc[1])-2) * log((double)RADIX) + 
         log((double)FIRST(int_polyc[1])))/LOGB;
log_PART2 = ((double)(SIZE(int_polyc[0])-2) * log((double)RADIX) + 
         log((double)FIRST(int_polyc[0])))/LOGB;

if (DEBUG_REF_INT) (void) printf("Part1, Part2, log_lhs = %g %g %g\n",
                          log_PART1, log_PART2, log_pmax_lhs);

degree = (double)fdegree;
fconst = (double)poly_constant;

//sqfprep();                  /* prepare quad. res. table for squfof   */

/*   This new code gets its assignments for the range of q's to do     */
/*   from a file called 'range'. It finds an unassigned range,         */
/*   marks it as 'in use', does the sieving, then marks the range as   */
/*   'done' and gets another range.                                    */


/*   ASCII_BASE only holds roots for the small primes, not their       */
/*  powers.  Roots modulo the powers may/may not exist depending upon  */
/*  whether they divide the discriminant. It is faster, since there are*/
/*  so few to just try to find a root by trial division                */

prime_power_roots();

while (1)  {                 /* loop forever                           */

numtot = 0;                  /* set counts of successes to 0           */
numff =  0;
numfp =  0;
numpf =  0;
numppq = 0;
numqpp = 0;
numfpp = 0;
numpp =  0;
numppf = 0;
numpppp = 0;

get_range(&starting_q,&ending_q,&line_num);   /* get the range         */
if (starting_q > size_int_base)
   {
   linroots = fopen("LINEAR_ROOTS","r");
   offset = starting_q - size_int_base;  
   for (i=0; i<=ending_q - size_int_base; i++)
	  {
      fscanf(linroots,"%d %d",&q1, &r1);
      if (i >= offset) 
          {  
          big_special_q[i-offset] = q1;
          big_roots[i-offset] = r1;
          }
      }
   fclose(linroots);
   }

_itoa(starting_q,namestring,10);    /* construct name of output file   */
strcpy(outname,HOME_DIR);   
strcat(outname,"nfsl");
strcat(outname,namestring);
strcat(outname,".out");

strcpy(restname,HOME_DIR);
strcat(restname,"nfsl");            /* construct name of restart file  */
strcat(restname,namestring);
strcat(restname,".res");

strcpy(rejname,HOME_DIR);
strcat(rejname,"nfsl");             /* construct name of reject file   */
strcat(rejname,namestring);
strcat(rejname,".rej");

/*         Open data and restart output files                          */

rejfile = fopen(rejname,"w");
outfile = fopen(outname,"r");         /* open output file              */
if (outfile == NULL)                  /* if none open for writing      */
   {
   outfile = fopen(outname,"w");
   (void) fprintf(outfile,"%s\n",label); /* put label at top of file   */
   }
else                                  /* a partial output file already */
   {                                  /* exists; append to it          */
   fclose(outfile);
   outfile = fopen(outname,"a");
   }
 
restart = 1;

/*   note: we partition by giving different special q's to different   */
/*  processors. These are specified by INDEX into the factor base from */
/*  the 'range' file.                                                  */
/*   Thus,  q is an index into the factor base                         */

/*      M A I N      S I E V I N G     L O O P                         */

if (GLOBAL_TIME_STATS) global_time = get_time();

if (DISPLAY_REDUCTIONS) display_reductions();

for (q=starting_q; q<=ending_q; q++)   /* Let's do it; loop over q's   */
   {
   if (SPECIAL_Q_TIME) special_q_time = get_time();
   if (TIME_STATS) total_q_time = get_time();
   if (TIME_STATS) init_time_variables();
   if (FALSECOUNT) init_count_variables();

/*   find the reduced lattice                                          */
   if (q > size_int_base)
      { 
      special_q = big_special_q[q - starting_q];             
      special_q_root = big_special_q[q - starting_q] - big_roots[q - starting_q]; 
//      special_q_root = int_fbase[q] - int_roots[q]; 
      log_special_q = (char)(log((double)special_q)/LOGB);
      reduce_lattice(special_q, special_q_root, v1, v2); 
      }
   else
      {
      special_q = int_fbase[q];
      special_q_root = int_fbase[q] - int_roots[q]; 
      log_special_q = (char)(log((double)special_q)/LOGB);
      reduce_lattice(special_q, special_q_root, v1, v2);   
//    newest_latred_exp(special_q, special_q_root, v1, v2);   
//    latred(special_q, special_q_root, v1, v2);
      } 

   if (reject_this_q())
      {
      (void) fprintf(rejfile,"%d %d\n", q, special_q);
      (void) fflush(rejfile);
      continue;
      } 

//   if (SHOW_MAIN) (void) printf("qindex, special_q, root, log, = %d %d %d %d\n",
//                         q, special_q, special_q_root, log_special_q);
//   if (SHOW_MAIN) (void) printf("lattice = (%d,%d), (%d,%d)\n",v1[0],v1[1],v2[0],v2[1]);

/*  find the starting points for this lattice and special_q            */

   find_start_points();

   siever(q,1);                          /* sieve positive a's         */
 
   reset_startpts();                     /* update starting points     */

   siever(q,-1);                         /* now sieve negative a's     */
 
   dump_params(q,restname,outname);      /* write the current q etc.   */
   totq++;                               /*  into the restart file     */

   if (FALSECOUNT)
      {
      (void) printf("Int LP bad split: %d\n",int_bad);
      (void) printf("Alg LP bad split: %d\n",alg_bad);
      (void) printf("Int cofactor prime: %d\n",int_prime);
      (void) printf("Alg cofactor prime: %d\n",alg_prime);
      (void) printf("Int big cofactor: %d\n",int_cofactor);
      (void) printf("Alg big cofactor: %d\n",alg_cofactor);
      (void) printf("Total algebraic hits: %d\n",num_alg_success);
      (void) printf("Total combined hits: %d\n",total_q_successes);
      (void) printf("Not coprime: %d\n",not_coprime);
      (void) printf("Relations found: %d\n",relation_count);
      }
   if (RHO_TALLY)   (void) printf("calls to rho = %d\n",rho_calls);
   if (Q_OVERFLOW)  (void) printf("squfof overflows = %d\n",overflow_count);
   if (SQUFOF_FAIL) (void) printf("squfof failures = %d\n",squfof_fail);
   if (SQUFOF)        (void) printf("squfof calls = %d\n",squfof_calls);
   if (TIME_STATS)
      {
      (void) printf("======================================================\n");
      total_q_time = get_time() - total_q_time;
      seconds = total_q_time/(clockrate*1000.0);
      (void) printf("Total time to process %d : %lf\n",q,seconds);
      seconds = (int_linesieve_time + alg_linesieve_time + int_vector_time + alg_vector_time);
      seconds = seconds/(clockrate*1000.0);
      (void) printf("Total sieve time = %lf\n",seconds);
      (void) printf("   Int line: %lf\n",int_linesieve_time/(clockrate*1000.0));
      (void) printf("   Alg line: %lf\n",alg_linesieve_time/(clockrate*1000.0));
      (void) printf("   Int vector: %lf\n",int_vector_time/(clockrate*1000.0));
      (void) printf("   Alg vector: %lf\n",alg_vector_time/(clockrate*1000.0));
      seconds = (int_resieve_time + alg_resieve_time)/(clockrate*1000.0);
      (void) printf("Total resieve time = %lf\n",seconds);
      (void) printf("   Int resieve: %lf\n",int_resieve_time/(clockrate*1000.0));
      (void) printf("   Alg resieve: %lf\n",alg_resieve_time/(clockrate*1000.0));
      (void) printf("Trial Int time: %lf\n",trial_int_time/(clockrate*1000.0));
      (void) printf("Trial Alg time: %lf\n",trial_alg_time/(clockrate*1000.0));
      (void) printf("Find startpts time: %lf\n",find_startpt_time/(clockrate*1000.0));
      (void) printf("Alg scan time: %lf\n",alg_scan_time/(clockrate*1000.0));
      (void) printf("Lattice reduce time: %lf\n",latred_time/(clockrate*1000.0));
      (void) printf("QS/Squfof time: %lf\n",squfof_time/(clockrate*1000.0));
      (void) printf("Prepare regions time: %lf\n", vector_setup_time/(clockrate*1000.));
      (void) printf("Inverse time = %lf\n",invtime/(clockrate*1000.0));
      (void) printf("Prime Test time = %lf\n",total_ptime/(clockrate*1000.0));
      (void) printf("======================================================\n");
      }
   if (SPECIAL_Q_TIME)
	  {
	  (void) printf("==========\n");
      (void) printf("Time for Q = %lf\n", (get_time() - special_q_time)/(clockrate*1000.));
      (void) printf("==========\n");
	  }
   if (COMPUTE_AREA) show_areas();
   }      /* end of loop over q's */

/*  N.B.  the vector sieve times and the resieve times include      */
/*  the setup times. To get just sieving time, subtract 1/4 of      */
/*  prepare regions time for each                                   */

if (GLOBAL_TIME_STATS)
   {
   seconds = (get_time() - global_time)/(clockrate*1000.0);
   (void) printf("Total time for %d special_q = %lf\n",totq, seconds);
   (void) printf("This is %lf per special_q\n", seconds/(double)totq);
   }

mark_range_as_used(starting_q,line_num);
(void) fclose(rejfile);
(void) fclose(outfile);
}      /*   end of while loop for getting range assignments         */
 
return(0);

}   /* end of main */

      
/************************************************************************/
/*                                                                      */
/*   Routine to scan the sieve array after it has been sieved           */
/*   on both sides.                                                     */
/*   the a,b pairs of factored locations are then dumped to the         */
/*   output file aint with the large primes.                            */
/*                                                                      */
/*   We need only check those locations marked in the success           */
/*   array.                                                             */
/*                                                                      */
/************************************************************************/

void sieve_scanner(sign)
int sign;
 
{   /* start of sieve_scanner */
register int i;
int alg_number[MPDIM], int_number[MPDIM];
int intrem1,intrem2,alg_rem,loc,alg_rem2,c,d, column;
double stime;

if (DEBUG) { (void) printf("In main sieve scanner\n"); } 

/*  first, go through the array line by line and mark the good locations*/
/*  we do this line by line because we don't have to change the         */
/*  threshhold value too often if we do.                                */

find_success_locations();
if (num_total_success == 0) return;          /*  I hope not!            */

if (SHOW_SUCCESS_COUNT) (void) printf("num_total_success = %d\n",num_total_success);

if (DEBUG) (void) printf("Factor by resieve\n");

factor_by_resieve(sign);                   /* save primes that hit      */
                                           /* the success locations     */
if (DEBUG) (void) printf("Trial factoring; num_total = %d\n",num_total_success);
fflush(stdout);
for (i=0; i<num_total_success; i++)
   {
   loc = sign*success[i][LOC_STORE];
   column = success[i][COLUMN_STORE];

    if (DUMPNORM) (void) printf("Testing %d : (%d,%d)\n",i,loc,column);

    compute_coeffs(loc, column, &c, &d);   /* transform of (i,j)        */
                                           /* (c,d) = i V1 + j V2       */
   if (d < 0) { c = -c; d = -d; }          /* output convention         */

   if (DUMPNORM) (void) printf("(c,d) =  (%d,%d)\n",c,d);

   if (!relprime(abs(c),d)) { not_coprime++; continue;} /* skip if (c,d)!=1*/
   if (HIT_COUNT) total_hits++;

   if (DUMPNORM) (void) printf("Computing Alg Norm\n");

   alg_norm(d,c,alg_number);
   if (DUMPNORM)
      { 
      (void) printf("i,j= %d %d\n",loc, column);
      (void) printf("c,d= %d %d\n",c,d);
      (void) printf("initial alg norm = ");
      write_num(alg_number,numstr);
      (void) fflush(stdout);
      }
   if (TIME_STATS) stime = get_time();
   alg_rem = trial_alg_div(alg_number,sign,&alg_rem2,loc,column);
   if (TIME_STATS) trial_alg_time += (get_time() - stime);

   if (alg_rem < BPRIME_LIMIT && alg_rem2 < BPRIME_LIMIT)
                                      /* If algebraic side factored     */
      {                               /* try to factor integer side     */
      int_norm(d,c,int_number);       /* compute norm of a+b M          */
      if (DUMPNORM)
        {
        (void) printf("initial int norm = ");
        write_num(int_number,numstr);
        (void) fflush(stdout);
        }
      if (TIME_STATS) stime = get_time(); 
      trial_int_div(int_number,sign,&intrem1,&intrem2,loc,column);
      if (TIME_STATS) trial_int_time += (get_time() - stime); 
      }
 
/*      Check that the large primes (if any) are not too big            */

   if (SHOW_LARGE_PRIMES)
      {
      (void) printf("Large primes = %d %d %d %d\n",intrem1,intrem2,alg_rem,alg_rem2);
      fflush(stdout);
      }
   if (intrem1 > BPRIME_LIMIT || intrem2 > BPRIME_LIMIT || 
      alg_rem > BPRIME_LIMIT || alg_rem2 > BPRIME_LIMIT)
      {
      if (FALSESHOW)                    /* show false hits?            */
        {
        (void) printf("ALG big %d %d: %d %d\n",c,d, alg_rem,alg_rem2);
        (void) printf("INT big %d %d: %d %d\n",c,d,intrem1,intrem2);
        (void) fflush(stdout);
        }
      continue;                        /* too big, so ignore           */
      }
   relation_count++; 
   output_relation(intrem1,intrem2,alg_rem,alg_rem2,sign,c,d);
   }      /* end of i loop */

}   /* end of sieve_scanner */
 

/************************************************************************/
/*                                                                      */
/*   Routine to compute a norm on the integer side: a + bX              */
/*   b > 0                                                              */
/************************************************************************/
 
void int_norm(b,a, answer)
int a,b;
int answer[];
 
{   /* start of int_norm */
int sign, temp[MPDIM];

sign = SIGN(a);
a = abs(a);
 
mul_single(int_polyc[1],b,answer);
mul_single(int_polyc[0],a,temp);

if (sign > 0) add(answer,temp,answer);
else if (mpcmp(answer,temp) == 1) subt(answer,temp,answer);
else subt(temp,answer,answer);

}   /* end of int_norm */

 
/************************************************************************/
/*                                                                      */
/*   Routine to handle the trial division of a rational integer         */
/*                                                                      */
/************************************************************************/
 
void trial_int_div(number,sign,rem1,rem2,loc,column)
int number[],sign,loc,column;
int *rem1,*rem2;
 
{   /* start of trial_int_div */
register int i,rem,tval1;
int temp[MPDIM],result;
int primecount,count;
int index;
double stime,ptime;

primecount = 0;
int_factor_list[0] = 0;
*rem1 = 0;
*rem2 = 0;

//  pull out vector primes first;  these are known to divide a priori
//  it would make the cofactor slightly smaller and hence div_single
//  would run slightly faster???

/*   First do small primes; we can't determine if they divide just by   */
/*   looking at the sieve location, so we must do full trial division   */
/*   this is because we sieved with respect to their powers not the     */
/*   primes themselves.                                                 */

if (INT_DEBUG) { (void) printf("INT candidate: "); write_num(number,numstr); }

i = 0;
while (i < NUMSMALLP)
   {
   rem = div_single(number,int_fbase[i],temp);
   if (rem == 0)
      {
      primecount++;
      int_factor_list[primecount] = int_fbase[i];
      int_factor_list[0] = primecount;
//    mpcopy(temp,number);
      small_mpcopy(temp,number);
      if (INT_SMALL) (void) printf("div at: %d\n",int_fbase[i]);
      }
   else i++;
   }
 
/*   And now do the linesieve primes;                                   */
/*   If the loc = -column * root mod p,  then p divides.                */
/*   compute (loc + column*root) except proj primes, do those separately*/

if (INT_SMALL) { (void) printf("doing linesieve primes\n"); }
if (INT_SMALL) { (void) printf("column = %d\n",column); }

while (i < int_linesieve)
   {
   tval1 =  (loc + column * special_q_int_roots[i]) % int_fbase[i];
   if (tval1 == 0)                           /* successful divide       */
      {
      if (INT_SMALL) 
         { (void) printf("div at: %d\n",int_fbase[i]); fflush(stdout); }
      rem = div_single(number,int_fbase[i],temp); /* may not div twice  */
      if (rem == 0)                         /* so check remainder       */
         {
//       mpcopy(temp,number);               /* update cofactor          */
         small_mpcopy(temp,number);
         if (INT_SMALL)
            { 
            (void) printf("success; number = "); 
            write_num(number,numstr); 
            fflush(stdout);
            }
         primecount++;
         int_factor_list[primecount] = int_fbase[i];
         int_factor_list[0] = primecount;   
         }
      else i++;
      }
   else i++; 
   }

/*      Now do the projective primes, if any                           */

i = 1;
while (i<=proj_int_primes[0])
   {
   rem = div_single(number, int_fbase[proj_int_primes[i]], temp);
   if (rem == 0)                             /* so check remainder     */
      {
//    mpcopy(temp,number);                  /* update cofactor         */
      small_mpcopy(temp,number);
      if (INT_SMALL)
         { 
         (void) printf("int proj success; %d %d %d number = ",
                        i,proj_int_primes[0],int_fbase[proj_int_primes[i]]); 
         write_num(number,numstr); 
         fflush(stdout);
         }
      primecount++;
      int_factor_list[primecount] = int_fbase[proj_int_primes[i]];
      int_factor_list[0] = primecount;
      }
   else i++;
   }

/*      Now add the vector primes to the factor list; special_q first   */

rem = div_single(number,special_q,temp);

if (DEBUG_Q)
   {
   if (rem != 0)
      {
      (void) printf("special q div failure! %d %d %d\n",special_q, sign, loc);
      (void) printf("norm = "); write_num(number,numstr);
      }
   }
while (rem == 0)
   {
// mpcopy(temp,number);
   small_mpcopy(temp,number);
   primecount++;   
   int_factor_list[primecount] = special_q;
   int_factor_list[0] = primecount;
   rem = div_single(number,special_q,temp); /* maybe div more than once?*/
   }

/*  now we deal with any big primes stored via factor-by-resieve        */

index = int_hash(abs(loc),column);
count = int_faclist[index][0];   
if (INT_DEBUG)
   {
   (void) printf("loc,col, hash, count = %d %d %d %d\n",loc,column,index,count);
   for (i=1; i<=count; i++) (void) printf("stored: %d\n",int_faclist[index][2+i]);
   fflush(stdout);
   }

i = 0;
while (i<count)
   {
   if (INT_DEBUG) (void) printf("trying %d\n", int_faclist[index][3+i]);

   rem = div_single(number,int_faclist[index][3+i],temp);
   if (rem == 0)
      {
     if (INT_DEBUG) { (void) printf("success: number = "); write_num(temp,numstr); }
//   mpcopy(temp,number);
     small_mpcopy(temp,number);
     primecount++;
     int_factor_list[primecount] = int_faclist[index][3+i];
     int_factor_list[0] = primecount;
     continue;                    /* don't increase i, may div twice   */
     }                            /* else if rem != 0  BUG!!!!         */
   i++;
   }

/*      Now, all primes in the factor base have been tried. If the     */
/*      remaining cofactor is less than MAX_SPLIT^2 then check if      */
/*      it is prime. If it is prime and less than BPRIME_LIMIT it is   */
/*      a single large prime. If composite, try to factor it with      */
/*      SQUFOF.   If it is single precision we are done.               */

if (INT_DEBUG) { (void) printf("After divisions: "); write_num(number,numstr); }

if (SIZE(number) == SINGLE)                     /* check cofactor      */
   {                                            /* if single prec      */
   if (FIRST(number) == 1)
//   if (FIRST(number) <= int_pmax)             /* we are done         */
     {
     *rem1 = 1;                                 /* fully factored      */   
     return;
     }
   *rem1 = FIRST(number);                       /* one large prime     */
   *rem2 = 1;
   int_factor_list[0] = primecount;
   if (FALSECOUNT) 
      {
      if (*rem1 > BPRIME_LIMIT) int_prime++;  
      }
   return;
   }   

if (INT_DEBUG)
   {
   (void) printf("done w/ div: number = "); write_num(number,numstr);
   fflush(stdout);
   }
if (SIZE(number) > DOUBLE || mpcmp(number,bpsquared) > 0)
   {
   *rem1 = BPRIME_LIMIT+1; 
   if (FALSECOUNT) int_cofactor++;            /* cofactor too big      */
   return;                           
   }

if (mpcmp(number,int_pmax_squared) == -1)      /* cofactor < pmax^2    */
   {                                           /* thus it must be prime*/
   if (FALSECOUNT) int_prime++;                /* since it isn't single*/
   *rem1 = BPRIME_LIMIT+1;                     /* it is too big        */                        
   }

/*   OK,  the number is greater than pmax^2 and less than SPLIT_LIMIT^2*/

if (TIME_STATS) ptime = get_time();
result = fast_ptest(number);
if (TIME_STATS) total_ptime += (get_time() - ptime);
if (result == 1)                             /* if cofactor is prime   */               
   { 
   *rem1 = BPRIME_LIMIT+1;
   if (FALSECOUNT) int_prime++;              /* prime cofac too big    */
   return;
   }

/*      OK. The number is composite and less than SPLIT_LIMIT^2.       */
/*      Try   to do it by squfof/qs.                                   */

if (INT_DEBUG)
   {
   (void) printf("squfof is trying: ");
   write_num(number,numstr);
   fflush(stdout);
   }

if (TIME_STATS) stime = get_time();
result = squfof(number,rem1,rem2);
if (TIME_STATS) squfof_time += (get_time() - stime);

if (DEBUG) { (void) printf("squfof returns: %d %d %d\n",result,*rem1, *rem2); }

if (result == 0)                                        /* no luck   */
   {
   *rem1 = BPRIME_LIMIT+1;
   if (SQUFOF_FAIL) { (void) printf("squfof fail\n"); fflush(stdout); }
   return;
   }
if (FALSECOUNT)              /* factors returned by squfof too big   */
   {
   if (*rem1 > BPRIME_LIMIT || *rem2 > BPRIME_LIMIT) int_bad++;
   }
int_factor_list[0] = primecount;

}   /* end of trial_int_div */
 
 
/************************************************************************/
/*                                                                      */
/*   Routine to handle the trial division of an integer;                */
/*   algebraic factor base                                              */
/*                                                                      */
/************************************************************************/
 
trial_alg_div(number,sign,alg_rem2,loc,column)
int number[],sign,*alg_rem2,loc,column;
 
{   /* start of trial_alg_div */
register int i,j,rem,tval1;
int temp[MPDIM],count;
int primecount,alg_rem,result,index;
double stime,ptime;

primecount = 0;
alg_factor_list[0] = 0;

*alg_rem2 = 0;

/*   first do small primes; we can't determine if they divide just by   */
/*   looking at the sieve location, so we must do full trial division   */

if (ALG_DEBUG) 
   {
   (void) printf("ALG trial candidate = "); write_num(number,numstr); 
   fflush(stdout); 
   }
i = 0;
while (alg_fbase[i] <= SMALLP)
   {
   rem = div_single(number,alg_fbase[i],temp);
   if (ALG_SMALL) { (void) printf("tried %d, rem = %d\n",alg_fbase[i],rem); }
   if (rem == 0)
     {
     primecount++;
     alg_factor_list[primecount] = alg_fbase[i];
     alg_factor_list[0] = primecount;
//   mpcopy(temp,number);
     small_mpcopy(temp,number);
     if (ALG_SMALL) (void) printf("success at: %d\n",alg_fbase[i]);
     }
   else i++;
   }
 
/*            And now do the rest of the factor base                   */
/*  If the loc = -column * root mod p,  then p divides.                */
/*   compute (loc + column*root)except proj primes, do those separately*/


if (ALG_DEBUG) { (void) printf("doing proj primes\n"); fflush(stdout); }

/*   First do the projective primes, if any                            */

j = 1;
while (j<=proj_alg_primes[0])
   {
   rem = div_single(number, alg_fbase[proj_alg_primes[j]], temp);
   if (rem == 0)                             /* so check remainder     */
      {
//    mpcopy(temp,number);                  /* update cofactor         */
      small_mpcopy(temp,number);
      if (ALG_SMALL)
         { 
         (void) printf("alg proj success; %d %d %d number = ",j, 
                  proj_alg_primes[0],alg_fbase[proj_alg_primes[j]]);
 
         write_num(number,numstr); 
         fflush(stdout);
         }
      primecount++;
      alg_factor_list[primecount] = alg_fbase[proj_alg_primes[j]];
      alg_factor_list[0] = primecount;
      }
   else j++;
   }

if (ALG_DEBUG) { (void) printf("doing linesieve primes\n"); fflush(stdout); }

//   Note: can not get multiplication overflow just below because these primes
//   are small. (column*special_q_alg_roots[i])

while (i < alg_linesieve)
   {
   tval1 =  (loc + column * special_q_alg_roots[i]) % alg_fbase[i];
   if (tval1 == 0)                             /* successful divide    */
      {
      if (ALG_DEBUG) { printf("try div %d\n",alg_fbase[i]); fflush(stdout);}
      rem = div_single(number,alg_fbase[i],temp);  /* may not div twice*/
      if (rem == 0)                            /* so check remainder   */
         {
//       mpcopy(temp,number);                  /* update cofactor      */
         small_mpcopy(temp,number);
         if (ALG_DEBUG)
            { 
            (void) printf("success; number = "); 
            write_num(number,numstr); 
            fflush(stdout);
            }
         primecount++;
         alg_factor_list[primecount] = alg_fbase[i];
         alg_factor_list[0] = primecount;
         if (SIZE(number) == SINGLE)               /* check cofactor   */
            {                                      
            if (FIRST(number) < alg_squares[i])    /* if < p^2 we are  */
               {                                   /* done; cofactor   */
               if (FIRST(number) <= alg_pmax)      /* is in fbase      */      
                  {
                  primecount++;
                  alg_factor_list[primecount] = FIRST(number);
                  alg_factor_list[0] = primecount;
                  return(1);
                  }
               else
                  {                        /* cofactor is LP   */
                  alg_factor_list[0] = primecount;
                  if (FALSECOUNT) if (FIRST(number) > BPRIME_LIMIT) alg_prime++;
                  return(FIRST(number));
                  }
               }  /* check p^2 */
            }
         continue;         /* don't increase i, see if divides twice   */
         }                 /* rem = 0? */
      }
   i++;
   }

/*      Now check all the primes stored by the resieving process       */

index = alg_hash(abs(loc),column);
count = alg_faclist[index][0];   

if (ALG_DEBUG)
   {
   (void) printf("loc,col, hash, count = %d %d %d %d\n",loc,column,index,count);
   for (i=1; i<=count; i++) (void) printf("stored: %d\n",alg_faclist[index][2+i]);
   fflush(stdout);
   }

i = 0;
while (i<count)
   {
   if (ALG_DEBUG)
      {
      (void) printf("trying %d\n", alg_faclist[index][3+i]);
      }
   rem = div_single(number,alg_faclist[index][3+i],temp);
   if (rem == 0)
      {
     if (ALG_DEBUG) { (void) printf("success: number = "); write_num(number,numstr); }
//   mpcopy(temp,number);
     small_mpcopy(temp,number);
     primecount++;
     alg_factor_list[primecount] = alg_faclist[index][3+i];
     alg_factor_list[0] = primecount;
     continue;                /* don't increase i, may div twice      */
     }                        /* else if rem != 0  BUG!!!!            */
   i++;
   }
 
/*   Now, all primes in the factor base have been tried. If the       */
/*   remaining cofactor is less than BPRIME_LIMIT^2 then check if     */
/*   it is prime. If it is prime and less than BPRIME_LIMIT it is     */
/*   a single large prime. If composite, try to factor it with        */
/*   SQUFOF.   If it is single precision we are done.                 */

if (ALG_DEBUG) { (void) printf("After divisions: "); write_num(number,numstr); }

if (SIZE(number) == SINGLE)                       /* check cofactor   */
   {                                              /* if single prec   */
//   if (FIRST(number) <= alg_pmax) return(1);       /* we are done   */  
   if (FIRST(number) == 1) return(1);
   if (FALSECOUNT) if (FIRST(number) > BPRIME_LIMIT) alg_prime++;
   alg_factor_list[0] = primecount;
   return(FIRST(number));                        /* one large prime   */
   }   

if (SIZE(number) > DOUBLE || mpcmp(number,bpsquared) > 0)
   {
   if (FALSECOUNT) alg_cofactor++;               /* cofactor too big  */
   return(BPRIME_LIMIT+1);                  
   }

if (mpcmp(number,alg_pmax_squared) == -1)     /* cofactor < pmax^2    */
   {                                          /* thus it must be prime*/
   if (FALSECOUNT) alg_prime++;               /* since it isn't single*/
   return(BPRIME_LIMIT+1);                    /* it is too big        */                           
   }

if (ALG_DEBUG) (void) printf("calling ptest\n");
if (TIME_STATS) ptime = get_time();
result = fast_ptest(number);
if (TIME_STATS) total_ptime += (get_time() - ptime);
if (result == 1)                              /* if cofactor is prime */               
   {                                          /* then it is too big   */
   if (FALSECOUNT) alg_prime++;
   return(BPRIME_LIMIT+1);
   }

/*      OK. The number is composite and less than SPLIT_LIMIT^2.      */
/*      Try to do it by squfof.                                       */

if (ALG_DEBUG) 
   {
   (void) printf("squfof is trying: ");
   write_num(number,numstr);
   }

if (TIME_STATS) stime = get_time();
result = squfof(number,&alg_rem,alg_rem2);
if (TIME_STATS) squfof_time += (get_time() - stime);

if (ALG_DEBUG) { (void) printf("alg squfof: %d %d %d\n",result,alg_rem, *alg_rem2);}

if (result == 0)                               /* no luck             */
   {
   if (SQUFOF_FAIL) (void) printf("alg squfof fail\n");
   return(BPRIME_LIMIT+1);
   }
if (FALSECOUNT)               /* factors returned by squfof too big   */
   {
   if (alg_rem > BPRIME_LIMIT || *alg_rem2 > BPRIME_LIMIT) alg_bad++;
   }
alg_factor_list[0] = primecount;
return(alg_rem);

}   /* end of trial_alg_div */
 
 
/************************************************************************/
/*                                                                      */
/*   Routine to scan the sieve array after the algebraic side has       */
/*   been sieved. We count up the number of success and store           */
/*   their locations in the 'successes' array. We reset the sieve       */
/*   values for the start of the integer side sieveing. Every           */
/*   Non-success location is set to a large negative number. Every      */
/*   success location is set to zero.                                   */
/*   After the integer side is sieved, we need only scan the            */
/*   locations stored in the success array to check for relations       */
/*                                                                      */
/************************************************************************/
 
void do_delayed_alg_resieve_points ();

void algebraic_scan(sign) 
int sign;

{   /* start of algebraic_scan */
register int i, j, k, count, column;
int *ptr, add_val, test_val, check;
int refine_count,ival;
double stime;
int testcount;
double test_sum;

/*   if there are any projective roots that divide column, then update */
/*  the initial value; such primes divide every location in the column */
/*  otherwise the init value is just log(special-q)                    */

if (TIME_STATS) stime = get_time();
if (DEBUG)
   {
   testcount = 0;
   test_sum = 0;
   }

if (DEBUG) printf("Algebraic scan\n");

count = 0;

for (column=1; column <= NCOLUMNS; column++)
   {
   ival = log_special_q;
   for (i=1; i<=proj_int_primes[0]; i++)
      if (column % int_fbase[proj_int_primes[i]] == 0) ival += int_logp[proj_int_primes[i]];
   refine_count = 0;
   sieve = sieve_array[column];
   ptr = (int*)sieve_array[column];
   test_alg_val = refine_alg_log(column, 1, sign);
   if (DEBUG_SCANNER) 
      { (void) printf("column alg_test_val = %d %d\n",column, test_alg_val); fflush(stdout);}
   add_val = 128 - test_alg_val;
   test_val = (add_val << 24) + (add_val << 16) + (add_val << 8) + add_val;
   if (DEBUG) { test_sum+= test_alg_val; testcount++; }
   if (DEBUG_SCAN_VALS) 
      {
      (void) printf("ALG SCAN; column = %d\n",column);
      for (j=0; j<100; j++) (void) printf("%d ",sieve[j]);
      (void) printf("\n");
      (void) fflush(stdout);
      }

   for (i=0; i<(MAX_SIEVE>>2); i+=4)
      {
      refine_count++;
      if (refine_count == ALG_REFINE_FREQ)      /* refine threshold value      */
         {
         test_alg_val = refine_alg_log(column, 4*i, sign);
         add_val = 128 - test_alg_val;
         test_val = (add_val << 24) + (add_val << 16) + (add_val << 8) + add_val;
         if (DEBUG) { test_sum+= test_alg_val; testcount++; }
         }
      check = (ptr[i] + test_val) | (ptr[i+1] + test_val) |
               (ptr[i+2] + test_val) | (ptr[i+3] + test_val);
      if ((check & MASK) != 0)
         {
         for (k=i; k<i+4; k++)
            {
            if (((ptr[k]+test_val) & MASK) == 0)
               {
               ptr[k] = RESET; 
               continue;
               }
            else for (j=4*k; j<4*k+4; j++)
               {
               if (sieve[j] > test_alg_val)
                  { 
                  if (DEBUG_SCANNER) 
                     { 
                     if (count > SUCCESS_SIZE) 
                         {(void) printf("success array overflow at %d\n",count); exit(0); }
                     (void) printf("alg succ #%d at: %d %d\n",count, j, column); 
                     (void) printf("sieve val, test, new val = %d %d %d\n",sieve[j], test_alg_val,ival);
                     }
                  success[count][LOC_STORE] = j;
                  success[count][COLUMN_STORE] = column;
                  count++;
                  sieve[j] = ival;      
                  }
               else sieve[j] = (-BIGNUM);
               }
            }
         }
      else
         {  
         ptr[i]   = RESET;
         ptr[i+1] = RESET;
         ptr[i+2] = RESET;
         ptr[i+3] = RESET; 
         }
      }
   }   /* end of loop over columns */
 
num_alg_success = count;
 
if (SHOW_SUCCESS_COUNT) (void) printf("after scan, num_alg_success = %d\n",count);

do_delayed_alg_resieve_points ();

if (TIME_STATS) alg_scan_time += (get_time() - stime);
if (DEBUG)
   (void) printf("testcount, sum, average = %d %g %g\n",testcount, test_sum, test_sum/(double)testcount);
}   /* end of algebraic_scan */
 
/************************************************************************/
/*                                                                      */
/*      Sieve initializer; Highly unrolled.                             */
/*   N.B.   MAX_SIEVE   *MUST* be a multiple of 128                     */
/*                                                                      */
/************************************************************************/
 
void sieve_init()
 
{   /* start of sieve_init */

register int i, limit, ival, column;
int *sieve_val;

/*   if there are any projective roots that divide column, then update  */
/*  the initial value; such primes divide every location in the column  */
/*  N.B.  MAX_SIEVE better be exactly divisible by 128                  */

for (column = 0; column<=NCOLUMNS; column++)
   {
   ival = 0;
   for (i=1; i<=proj_alg_primes[0]; i++)
      { 
      if ((column % alg_fbase[proj_alg_primes[i]]) == 0) ival += alg_logp[proj_alg_primes[i]];
      }
   ival = (ival << 24 + ival << 16 + ival << 8 + ival);
   limit = (MAX_SIEVE >> 7);
   sieve_val = (int*) sieve_array[column];
   for (i=0; i<limit; i++)
      {
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      *(sieve_val++) = ival; *(sieve_val++) = ival;
      }
   }
}   /* end of sieve_init */
 
 
/************************************************************************/
/*                                                                      */
/*      Routine to perform the sieving                                  */
/*      line sieve the smallest primes                                  */
/*                                                                      */
/*   primes are partitioned into sets according to how many hits they   */
/*   take in each sub-interval. Those taking 4 or more have their loop  */
/*   partially unrolled.  Those taking (3 or 4),  (2 or 3), (1 or 2),   */
/*   hits have their loops completely unrolled.                         */
/*   speedup? check storage of temporary variables in loops             */
/*                                                                      */
/************************************************************************/

void siever(qindex,sign)
int qindex,sign;
 
{   /* start of siever */
register int j,incr,column;
unsigned char *sieve1, *sieve2, *sieve3,  val;
register int fourp;
register int stop,loc;
double itime,atime;

if (DEBUG)
   {
   (void) printf("In MAIN siever, sign = %d\n",sign);
   (void) printf("Now sieving special_q = %d %d\n",special_q, qindex);
   (void) printf("Vectors = (%d,%d)  (%d,%d)\n", v1[0], v1[1], v2[0], v2[1]);
   }

/*   We sieve over (i,column)  for (|i|, column) < (MAX_SIEVE, NCOLUMNS).  */
/*   Primes that   are >= MAX_SIEVE may not take a hit for a given column. */
/*   These, we sieve by vectors.  The smaller primes (those that take at   */
/*  least one hit, we line sieve, one column at a time.                    */

sieve_init();                         /* initialize sieve array            */
vector_alg_sieve(sign);               /* sieve big primes by vectors       */

/*   Now do line sieving for the small primes                              */

if (DEBUG_SIEVE) (void) printf("Line sieving algebraic primes\n");

for (column=1; column < NCOLUMNS; column++)  {  

if (DEBUG_SIEVE) (void) printf("sieving column %d\n",column);

if (!(column & 1)) 
   {
   even_alg_siever(column,sign);
   update_alg_start_points(sign);
   continue;
   }

if (TIME_STATS) atime = get_time();
sieve = sieve_array[column];              /* ptr to current column         */ 

/*      Sieve the algebraic side first; then scan for success.             */
/*      Start with the very small primes;                                  */

j = 0;
while (alg_fbase[j] <= SMALLP)
   {  
   loc = alg_startpts[j];
   val = alg_logp[j];
   incr = alg_smallp[alg_fbase[j]];
   fourp = incr << 2;
   stop = MAX_SIEVE - 2*fourp + incr;
   sieve1 = sieve + incr;
   sieve2 = sieve1 + incr;
   sieve3 = sieve2 + incr;
   while (loc < stop)
      {
      sieve[loc] += val;
      sieve1[loc] += val;
      sieve2[loc] += val;
      sieve3[loc] += val;
      loc += fourp;
      sieve[loc] += val;
      sieve1[loc] += val;
      sieve2[loc] += val;
      sieve3[loc] += val;
      loc += fourp;
      }
   while (loc < MAX_SIEVE)      
      {
      sieve[loc] += val;
      loc += incr;
      }
   j++;
   }

for (; j<alg_8hits; j++)      /* these primes take >= 8 hits in the   */
   {                          /* interval; unroll 8x                  */
   loc = alg_startpts[j];     /* starting sieve location              */
   val = alg_logp[j];         /* scaled log of the prime              */
   incr = alg_fbase[j];       /* stride                               */
   fourp = incr << 2;   
   stop = MAX_SIEVE - 2*fourp + incr;      /* stop location for sieve */
   sieve1 = sieve + incr;
   sieve2 = sieve1 + incr;
   sieve3 = sieve2 + incr; 
   while (loc < stop)
      {
      sieve[loc] += val;  
      sieve1[loc] += val; 
      sieve2[loc] += val; 
      sieve3[loc] += val; 
      loc += fourp;
      sieve[loc] += val; 
      sieve1[loc] += val; 
      sieve2[loc] += val; 
      sieve3[loc] += val; 
      loc += fourp;
      }
   while (loc < MAX_SIEVE)     /* note: maybe unroll this? check later.*/
      {                        /* should execute max of 3 times        */
      sieve[loc] += val;
      loc += incr;
      }
   }

for (; j<alg_4hits; j++)      /* these primes take >= 4 hits in the    */
   {                          /* interval; unroll 4x                   */
   loc = alg_startpts[j];
   val = alg_logp[j];
   incr = alg_fbase[j];
   fourp = incr << 2;
   stop = MAX_SIEVE - fourp + incr;
   sieve1 = sieve + incr;
   sieve2 = sieve1 + incr;
   sieve3 = sieve2 + incr; 
   while (loc < stop)
      {
      sieve[loc] += val; 
      sieve1[loc] += val;
      sieve2[loc] += val;
      sieve3[loc] += val;
      loc += fourp;
      }
   while (loc < MAX_SIEVE)      /* note: maybe unroll this? check later.*/
      {                         /* should execute max of 3 times        */
      sieve[loc] += val;
      loc += incr;
      }
   }

for (; j<alg_3hits; j++)      /* these primes take 3 or 4 hits in the   */
   {                          /* interval; unroll 3x                    */
   loc = alg_startpts[j]; 
   val = alg_logp[j];
   incr = alg_fbase[j];
   sieve[loc] += val; 
   loc += incr;
   sieve[loc] += val;
   loc += incr;
   sieve[loc] += val;
   loc += incr;
   if (loc < MAX_SIEVE)       /* perhaps 1 more?                      */
      {
      sieve[loc] += val;
      }
   }

for (; j<alg_2hits; j++)      /* these primes take 2 or 3 hits        */
   {                          /* unroll twice                         */
   loc = alg_startpts[j];
   sieve[loc] += alg_logp[j]; 
   loc += alg_fbase[j];
   sieve[loc] += alg_logp[j];
   loc += alg_fbase[j];
   if (loc < MAX_SIEVE)        /* do one more?                        */
      {
      sieve[loc] += alg_logp[j];
      }
   }

for (; j<alg_1hits; j++)       /* algebraic side  large primes        */
   {                           /* these take 1 or 2 hits              */
   loc = alg_startpts[j];
   sieve[loc] += alg_logp[j]; 
   loc += alg_fbase[j];
   if (loc < MAX_SIEVE)        /* do one more?                        */
      {
      sieve[loc] += alg_logp[j];
      }
   }

for (; j<alg_linesieve; j++)      /* these primes take AT MOST 1 hit  */
   {                              /* per line. Expect 1 hit per every */
   if (alg_startpts[j] < MAX_SIEVE)    /* FUZZY line                  */
      sieve[alg_startpts[j]] += alg_logp[j];
   }
   
/*   For all other primes,  we do not line sieve   at all.            */
/*   Instead, for such primes, we sieved by vectors   at the top      */

if (DEBUG) (void) printf("Finished alg sieve: %d\n",column);

update_alg_start_points(sign);
if (TIME_STATS) alg_linesieve_time += (get_time() - atime);

}   /* end of loop over columns */

/*   Now we scan the array for successes. Reinitialize the successful */
/*  points, set the unsuccessful to -oo, and record the success locs  */

algebraic_scan(sign);

if (DEBUG) (void) printf("in siever, after algebraic scan num succcess = %d\n",num_alg_success);

/*   Now we sieve the integer factor base                             */

vector_int_sieve(sign);            /* vector sieve the big primes     */

if (DEBUG) (void) printf("Doing LINE Int sieve\n");

for (column=1; column<=NCOLUMNS; column++) {

if (DEBUG) (void) printf("Column: %d\n",column);

if (!(column & 1))
   { 
   even_int_siever(column,sign); 
   update_int_start_points(sign);   
   continue; 
   }

if (TIME_STATS) itime = get_time();

sieve = sieve_array[column];         /* ptr to current column         */ 

for (j=0; j<NUMSMALLP; j++)         /* handle smallest primes         */
   {
    if (int_root_ok[j] == FALSE) {j++; continue;}
   loc = int_startpts[j];
   val = int_logp[j];
   incr = smallp[j];
   fourp = incr << 2;
   stop = MAX_SIEVE - 2*fourp + incr;
   sieve1 = sieve + incr;
   sieve2 = sieve1 + incr;
   sieve3 = sieve2 + incr;
   while (loc < stop)               /* unroll the loop                */
      {
      sieve[loc] += val;
      sieve1[loc] += val;
      sieve2[loc] += val;
      sieve3[loc] += val;
      loc += fourp;
      sieve[loc] += val;
      sieve1[loc] += val;
      sieve2[loc] += val;
      sieve3[loc] += val;
      loc += fourp;
      }
   while (loc < MAX_SIEVE)
      {
      sieve[loc] += val;
      loc += incr;
      }
   }

for (; j<int_8hits; j++)      /* these primes take >= 8 hits in the   */
   {                          /* interval; unroll 8x                  */
   loc = int_startpts[j];
   val = int_logp[j];
   incr = int_fbase[j];
   fourp = incr << 2;
   stop = MAX_SIEVE - 2*fourp + incr;
   sieve1 = sieve + incr;
   sieve2 = sieve1 + incr;
   sieve3 = sieve2 + incr;
   while (loc < stop)
      {
      sieve[loc] += val;
      sieve1[loc] += val;
      sieve2[loc] += val;
      sieve3[loc] += val;
      loc += fourp;
      sieve[loc] += val;
      sieve1[loc] += val;
      sieve2[loc] += val;
      sieve3[loc] += val;
      loc += fourp;
      }
   while (loc < MAX_SIEVE)    /* note: maybe unroll this? check later.*/
      {                       /* should execute max of 3 times        */
      sieve[loc] += val;
      loc += incr;
      }
   }

for (; j<int_4hits; j++)      /* these primes take >= 4 hits in the   */
   {                          /* interval; unroll 4x                  */
   loc = int_startpts[j];
   val = int_logp[j];
   incr = int_fbase[j];
   fourp = incr << 2;
   stop = MAX_SIEVE - fourp + incr;
   sieve1 = sieve + incr;
   sieve2 = sieve1 + incr;
   sieve3 = sieve2 + incr;
   while (loc < stop)
      {
      sieve[loc] += val;
      sieve1[loc] += val;
      sieve2[loc] += val;
      sieve3[loc] += val;
      loc += fourp;
      }
   while (loc < MAX_SIEVE)    /* note: maybe unroll this? check later.*/
      {                       /* should execute max of 3 times        */
      sieve[loc] += val;
      loc += incr;
      }
   }

for (; j<int_3hits; j++)      /* these primes take 3 or 4 hits in the */
   {                          /* interval; unroll 3x                  */
   loc = int_startpts[j]; if (loc > MAX_SIEVE) printf("BAD at %d %d %d\n",j,int_fbase[j], loc);
   val = int_logp[j];
   incr = int_fbase[j];
   sieve[loc] += val;
   loc += incr; if (loc > MAX_SIEVE) printf("BAD at %d %d %d\n",j,int_fbase[j], loc);
   sieve[loc] += val;
   loc += incr; if (loc > MAX_SIEVE) printf("BAD at %d %d %d\n",j,int_fbase[j], loc);
   sieve[loc] += val;
   loc += incr;
   if (loc < MAX_SIEVE)         /* perhaps 1 more?                   */
      {
      sieve[loc] += val;
      }
   }

for (; j<int_2hits; j++)      /* these primes take 2 or 3 hits       */
   {                          /* unroll twice                        */
   loc = int_startpts[j];
   sieve[loc] += int_logp[j];
   loc += int_fbase[j];
   sieve[loc] += int_logp[j];
   loc += int_fbase[j];
   if (loc < MAX_SIEVE)       /* do one more?                        */
      {
      sieve[loc] += int_logp[j];
      }
   }

for (; j<int_1hits; j++)         /* integer side  large primes       */
   {                             /* these take 1 or 2 hits           */
   loc = int_startpts[j];  if (loc > MAX_SIEVE) printf("BAD at %d %d %d\n",j,int_fbase[j], loc);
   sieve[loc] += int_logp[j];
   loc += int_fbase[j];
   if (loc < MAX_SIEVE)         /* do one more?                      */
      {
      sieve[loc] += int_logp[j];
      }
   }

for (; j<int_linesieve; j++)      /* these primes take AT MOST 1 hit */
   {                              /* per line. Expect 1 hit per every*/
   if (int_startpts[j] < MAX_SIEVE)   /* FUZZY line                  */
      sieve[int_startpts[j]] += int_logp[j];
   }
   
/*   For all other primes,  we do not line sieve   at all.           */
/*   Instead, for such primes, we sieved by vectors   at the top     */

update_int_start_points(sign);
if (TIME_STATS) int_linesieve_time += (get_time() - itime);

}   /* end of loop over columns */

/*   Now we scan the sieve array looking for successes. Sieve scanner */
/*   calls factor_by_resieve and the trial div routines               */

sieve_scanner(sign);

total_q_successes += num_total_success;

/*   If there were any projective primes we did not sieve them. instead */
/*  we sieved during initialization and set the logs to 0 and then      */
/*  treated them as ordinary primes.  We need to restore the log values */

restore_logs();   

}   /* end of siever */


/************************************************************************/
/*                                                                      */
/*      Routine to perform the sieving;   EVEN b; alg primes            */
/*                                                                      */
/*                                                                      */
/*   primes are partitioned into sets according to how many hits they   */
/*   take in each sub-interval. Those taking 4 or more have their loop  */
/*   partially unrolled.  Those taking (3 or 4),  (2 or 3), (1 or 2),   */
/*   (0 or 1) hits have their loops completely unrolled.                */
/*   speedup? check storage of temporary variables in loops             */
/*                                                                      */
/*   This handles the case when b = 0 mod 2.  We can skip sieving the   */
/*   even locations since (a,b) != 1. Thus, we double the strides.      */
/*                                                                      */
/************************************************************************/

void even_alg_siever(column,sign)
int column,sign;
 
{   /* start of even_alg_siever */
register int j;
unsigned char *sieve1, *sieve2, *sieve3, val;
int loc, incr, fourp, stop;
double stime;

if (DEBUG) (void) printf("In even alg siever\n");
 
sieve = sieve_array[column];
if (TIME_STATS) stime = get_time();

/*      Start with the very small primes; Since we only sieve with respect*/
/*      to their prime powers, there may not be a root.                   */

/*      Handle the prime 2 separately. If the starting loc is even        */
/*      then it will only hit even locations, so we can skip it           */
/*      entirely                                                          */

if (alg_smallroot[0] != -1)            /* if the root exists              */
   {
   loc = alg_startpts[0];
   if (!(loc & 1)) goto label1;        /* if startpoint isn't even        */
   val = alg_logp[0];
   incr = alg_smallp[alg_fbase[0]];
   fourp = incr << 2;
   stop = MAX_SIEVE - fourp + incr;
   sieve1 = sieve + incr;
   sieve2 = sieve1 + incr;
   sieve3 = sieve2 + incr;
   while (loc < stop)
      {
      sieve[loc] += val;
      sieve1[loc] += val;
      sieve2[loc] += val;
      sieve3[loc] += val;
      loc += fourp;
      }
   while (loc < MAX_SIEVE)            /* should execute max of 3 times   */
      {                               /* is it worth unrolling?          */
      sieve[loc] += val;
      loc += incr;
      }
   }

/*      Now handle the other small prime powers                          */
label1: ;
j = 1;
while (alg_fbase[j] <= SMALLP)
   {
   if (alg_smallroot[j] == -1) {j++;  continue;}         /* no root      */
   loc = alg_startpts[j];
   val = alg_logp[j];
   incr = alg_smallp[alg_fbase[j]];
   if (!(loc & 1)) loc += incr;   /* start point is even, increment it   */
   incr = incr << 1;              /* double the stride!!                 */
   fourp = incr << 2;
   stop = MAX_SIEVE - 2*fourp + incr;
   sieve1 = sieve + incr;
   sieve2 = sieve1 + incr;
   sieve3 = sieve2 + incr;
   while (loc < stop)
      {
      sieve[loc] += val;
      sieve1[loc] += val;
      sieve2[loc] += val;
      sieve3[loc] += val;
      loc += fourp;
      sieve[loc] += val;
      sieve1[loc] += val;
      sieve2[loc] += val;
      sieve3[loc] += val;
      loc += fourp;
      }
   while (loc < MAX_SIEVE)         /* should execute max of 3 times    */
      {                            /* is it worth unrolling?           */
      sieve[loc] += val;
      loc += incr;
      }
   j++;
   }

/*   This loop (and the corresponding one for the integer side) is     */
/*   the most expensive part of the entire code: sieving small primes  */


for (; j<alg_8hits; j++)       /* these primes take >= 4 hits in the   */
   {                           /* interval; unroll 4x. Note: normally  */
   loc = alg_startpts[j];      /* these take 8 hits; but here 2x stride*/
   val = alg_logp[j];
   incr = alg_fbase[j];
   if (!(loc & 1)) loc += incr;   /* start point is even, increment it */      
   incr = incr << 1;              /* double the stride                 */
   fourp = incr << 2;
   stop = MAX_SIEVE - fourp + incr;
   sieve1 = sieve + incr;
   sieve2 = sieve1 + incr;
   sieve3 = sieve2 + incr;
   while (loc < stop)
      {
      sieve[loc] += val;
      sieve1[loc] += val;
      sieve2[loc] += val;
      sieve3[loc] += val;
      loc += fourp;
      }
   while (loc < MAX_SIEVE)      /* note: maybe unroll this? check later. */
      {                         /* should execute max of 3 times         */
      sieve[loc] += val;
      loc += incr;
      }
   }


/*   We don't unroll larger primes in the "skip" even sieve because      */
/*  we don't know how many hits each one makes; instead we just dble     */
/*  all strides.                                                         */

for (; j<alg_1hits; j++)         /* algebraic side medium primes         */
   {
   loc = alg_startpts[j];  
   incr = alg_fbase[j];
   if (!(loc & 1)) loc += incr;
   incr = incr << 1;
   val = alg_logp[j];  
   while (loc < MAX_SIEVE)
      {
      sieve[loc] += val;
      loc += incr;
      }
   }

for (; j<alg_linesieve; j++)      /* these primes take AT MOST 1 hit   */
   {                              /* per line. Expect 1 hit per every  */
   if (alg_startpts[j] < MAX_SIEVE)       /* FUZZY line                */
      sieve[alg_startpts[j]] += alg_logp[j];
   }

if (TIME_STATS) even_alg_linesieve_time += (get_time() - stime);

}   /* end of even_alg_siever */
 


/************************************************************************/
/*                                                                      */
/*      Routine to perform the sieving;   EVEN b rational primes        */
/*      line sieve; small primes                                        */
/*                                                                      */
/*   primes are partitioned into sets according to how many hits they   */
/*   take in each sub-interval. Those taking 4 or more have their loop  */
/*   partially unrolled.  Those taking (3 or 4),  (2 or 3), (1 or 2),   */
/*   (0 or 1) hits have their loops completely unrolled.                */
/*   speedup? check storage of temporary variables in loops             */
/*                                                                      */
/*   This handles the case when b = 0 mod 2.  We can skip sieving the   */
/*   even locations since (a,b) != 1. Thus, we double the strides.      */
/*                                                                      */
/************************************************************************/

void even_int_siever(column,sign)
int column,sign;
 
{   /* start of even_int_siever */
register int j;
unsigned char *sieve1, *sieve2, *sieve3, val;
int loc, incr, fourp, stop;
double stime;

if (DEBUG) { (void) printf("In even int siever\n"); fflush(stdout); }
if (TIME_STATS) stime = get_time(); 
sieve = sieve_array[column];

/*   Handle the prime 2 separately                                */

loc = int_startpts[0];
if (loc & 1)                    /* only if startpt is odd         */
   {
   val = int_logp[0];
   incr = smallp[0] << 1;
   if (!(loc & 1)) loc += incr;
   incr = incr << 1;
   fourp = incr << 2;
   stop = MAX_SIEVE - fourp + incr;
   sieve1 = sieve + incr;
   sieve2 = sieve1 + incr;
   sieve3 = sieve2 + incr;
      while (loc < stop)            /* unroll the loop            */
        {
        sieve[loc] += val;
        sieve1[loc] += val;
        sieve2[loc] += val;
        sieve3[loc] += val;
        loc += fourp;
        }
   while (loc < MAX_SIEVE)
        {
        sieve[loc] += val;
        loc += incr;
        }
   }

for (j=1; j<NUMSMALLP; j++)         /* handle smallest primes      */
   {
   loc = int_startpts[j];
   val = int_logp[j];
   incr = smallp[j];
   if (!(loc & 1)) loc += incr;
   incr = incr << 1;
   fourp = incr << 2;
   stop = MAX_SIEVE - 2*fourp + incr;
   sieve1 = sieve + incr;
   sieve2 = sieve1 + incr;
   sieve3 = sieve2 + incr;
   while (loc < stop)               /* unroll the loop            */
      {
      sieve[loc] += val;
      sieve1[loc] += val;
      sieve2[loc] += val;
      sieve3[loc] += val;
      loc += fourp;
      sieve[loc] += val;
      sieve1[loc] += val;
      sieve2[loc] += val;
      sieve3[loc] += val;
      loc += fourp;
      }
   while (loc < MAX_SIEVE)
      {
      sieve[loc] += val;
      loc += incr;
      }
   }

for (; j<int_8hits; j++)      /* these primes take >= 4 hits in the   */
   {                          /* interval; unroll 4x. normally 8 hits */
   loc = int_startpts[j];     /* but we have 2x stride.               */
   val = int_logp[j];
   incr = int_fbase[j];
   if (!(loc & 1)) loc += incr;
   incr = incr << 1;
   fourp = incr << 2;
   stop = MAX_SIEVE - fourp + incr;
   sieve1 = sieve + incr;
   sieve2 = sieve1 + incr;
   sieve3 = sieve2 + incr;
   while (loc < stop)
      {
      sieve[loc] += val;
      sieve1[loc] += val;
      sieve2[loc] += val;
      sieve3[loc] += val;
      loc += fourp;
      }
   while (loc < MAX_SIEVE)      /* note: maybe unroll this? check later.*/
      {                         /* should execute max of 3 times        */
      sieve[loc] += val;
      loc += incr;
      }
   }

/*   We don't unroll larger primes in the "skip" even sieve because     */
/*  we don't know how many hits each one makes; instead we just dble    */
/*  all strides. We do, however, unroll the outer loop.                 */

for (; j<int_1hits; j++)         /* integer side medium primes          */
   {                             /* these take 0 or 1 hits              */
   loc = int_startpts[j];
   incr = int_fbase[j];
   if (!(loc & 1)) loc += incr;
   incr = incr << 1;
   val = int_logp[j];
   while (loc < MAX_SIEVE)
      {
      sieve[loc] += val;
      loc += incr;
      }
   }

for (; j<int_linesieve; j++)      /* these primes take AT MOST 1 hit   */
   {                              /* per line. Expect 1 hit per every  */
   if (int_startpts[j] < MAX_SIEVE)       /* FUZZY line                */
      sieve[int_startpts[j]] += int_logp[j];
   }

if (TIME_STATS) even_int_linesieve_time += (get_time() - stime);

}   /* end of even_int_siever */
 
/************************************************************************/
/*                                                                      */
/*   The following GROUP of routines manage the bucket sieve. Basically */
/*  as we stride through a parallelogram, rather than update sieve[d][c]*/
/*  at every stride, we store an entire vector of points (if we can) and*/
/*  individual points if we can not in a linked-list bucket. Then before*/
/*  emptying the bucket we use an in-cache quicksort to sort the entries*/
/*  by column.  Then, when we empty, it allows us to stride along a row */
/*  (C is row-major) to maintain good locality and cache access.        */
/*                                                                      */
/*   Kudos!!! to George Woltman for providing these routines.           */
/*  Except for subroutine headers, and indentation, I have left his     */
/*  style intact; including inter-mixing code, data declarations,       */
/*  defines, etc.                                                       */
/*                                                                      */
/************************************************************************/

//
// These routines accumulate vectors to use to update the sieve at a later
// time in a cache-friendly manner.  We use a bucket sort to update the
// find the individual points on all the vectors in sorted order.
//

// Defines to tune the delayed vector updating process.
//
// It is best if the vector update info fits in the L2 cache.  It may even be
// better if it is small enough to eliminate TLB cache misses.
//
// NOTE:  DO NOT select 396KB to 536KB as the jump from 2 byte indexes
// to 4 byte indexes makes this a no-man's zone.  The macros below will
// compensate for this, but you should be aware of the phenomenon.

#define MEM_FOR_VECTOR_UPDATES   (256*1024)      /* 256 kilobytes */

// This tells us when to switch from saving dense vectors to saving individual
// points.  This is an important crossover point!  Set it carefully.

#define NUM_VECTOR_UPDATE_PRIMES   9500

// We go to great lengths to cram as much data in the vector buffer.  The
// more vectors we can save, the more cache hits we will get when updating
// the sieve array.
//
// These are the structures and variables used to save the vector info.
//
// Note we try to keep structures a power of two size.  In calculating
// how many vectors we can support, we assume we'll get lots of vectors
// per small prime making the size of the saved_a1_b1 data unimportant.

struct vector_info 
   {
   short   c;
   short   i_diff;      // Add to base_i to get the "true" i
   };
struct a1_b1 
   {         // Save a1,b1 the dist between vector points
   short   a1;
   short   b1;
   };

// MAX_VECTORS calculations are a little complicated.
// First we must see if using a short results in more than 65535 vectors.
// Second we must see if using a long would result in less than 65535 vectors.

#define NUM_BUCKETS			(NCOLUMNS+1)
#define SHORT_MAX_VECTORS	((MEM_FOR_VECTOR_UPDATES-NUM_BUCKETS*2)/6)
#define LONG_MAX_VECTORS	((MEM_FOR_VECTOR_UPDATES-NUM_BUCKETS*4)/8)

// Define the data type for the linked list and calculate the true MAX_VECTORS.
// To conserve memory, use short if possible.

#if SHORT_MAX_VECTORS <= 65535
typedef unsigned short bucket_ptr;
#define MAX_VECTORS   SHORT_MAX_VECTORS
#elif LONG_MAX_VECTORS <= 65535
typedef unsigned short bucket_ptr;
#define MAX_VECTORS   65535
#else
typedef   unsigned long bucket_ptr;
#define MAX_VECTORS   LONG_MAX_VECTORS
#endif

// Guess at a maximum number of A1/B1 values we'll need to save

#define MAX_A1_B1   (MAX_VECTORS >> 2)

// Global variables for delayed vector updating

struct a1_b1 saved_a1_b1[MAX_A1_B1];          // Saved a1,b1 values
struct vector_info vectors[MAX_VECTORS+1];    // c,i for each vector
bucket_ptr next_bucket_ptrs[MAX_VECTORS+1];   // Next bucket for each vector
bucket_ptr buckets[NUM_BUCKETS] = {0};        // Buckets for sorting vectors
int   num_vectors = 0;
int   base_i = 0;

/************************************************************************/
/*                                                                      */
/* Routine to process the delayed vectors in sorted order. This         */
/* maximizes our cache hits in the sieve array.                         */
/*                                                                      */
/************************************************************************/
                                                                     
void do_delayed_vectors (char *logs)

{   /* start of do_delayed_vectors */
int   d;
double   itime;
static   int   num_sets = 0;

// If we haven't accumulated any vectors to process, then we are done

   if (num_vectors == 0) return;

   if (TIME_STATS) itime = get_time();

// Loop through all the buckets
   
   for (d = 0; d <= NCOLUMNS; d++) 
     {
      int   idx;
      idx = buckets[d];
      buckets[d] = 0;
      while (idx) 
      {
         int   i, c, next_d, next_idx;
         i = vectors[idx].i_diff;
         c = vectors[idx].c;
         sieve_array[d][c] += logs[base_i+i];
         //gwcount2++;

         c += saved_a1_b1[i].a1;
         next_d = d + saved_a1_b1[i].b1;
         next_idx = next_bucket_ptrs[idx];
         if (c >= 0 && c < MAX_SIEVE && next_d <= NCOLUMNS) 
         {
            vectors[idx].c = c;
            next_bucket_ptrs[idx] = buckets[next_d];
            buckets[next_d] = idx;
            }
         idx = next_idx;
       }
     }

// Reset count of saved vectors

   num_vectors = 0;

// Output statistics
/*
   if (TIME_STATS) 
      {
      itime = (get_time() - itime);
      num_sets++;
      (void) printf("Vector array #%d (%d,%d): %lf\n",
                       num_sets,gwcount,gwcount2,itime/(clockrate*1000.0));
      } 
*/
}   /* end of do_delayed_vectors */

/************************************************************************/
/*                                                                      */
/* Routine to process the delayed vectors in sorted order.  This        */
/* maximizes our cache hits in the sieve array.                         */
/*                                                                      */
/************************************************************************/

void do_delayed_int_resieve_vectors ()
{   /* start of do_delayed_int_resieve_vectors */
   int   d;

// If we haven't accumulated any vectors to process, then we are done

   if (num_vectors == 0) return;

// Loop through all the buckets
   
   for (d = 0; d <= NCOLUMNS; d++) {
      int   idx;
      idx = buckets[d];
      buckets[d] = 0;
      while (idx) {
         int   i, c, next_d, next_idx;
         i = vectors[idx].i_diff;
         c = vectors[idx].c;
         if (sieve_array[d][c] == GOOD && c != 0 && d != 0) {
            int   index;
            index = int_hash (c,d);
            int_count = int_faclist[index][0];
            int_faclist[index][1] = c;
            int_faclist[index][2] = d;
            int_faclist[index][3+int_count] = int_fbase[base_i+i];
            int_faclist[index][0]++;
         }
         c += saved_a1_b1[i].a1;
         next_d = d + saved_a1_b1[i].b1;
         next_idx = next_bucket_ptrs[idx];
         if (c >= 0 && c <= MAX_SIEVE && next_d <= NCOLUMNS) {
            vectors[idx].c = c;
            next_bucket_ptrs[idx] = buckets[next_d];
            buckets[next_d] = idx;
         }
         idx = next_idx;
      }
   }

// Reset count of saved vectors

   num_vectors = 0;
}    /* end of do_delayed_int_resieve_vectors */

/************************************************************************/
/*                                                                      */
/* Routine to process the delayed vectors in sorted order. This         */ 
/* maximizes our cache hits in the sieve array.                         */
/*                                                                      */
/************************************************************************/
void do_delayed_alg_resieve_vectors ()

{   /* start of do_delayed_alg_resieve_vectors */
   int   d;

// If we haven't accumulated any vectors to process, then we are done

   if (num_vectors == 0) return;

// Loop through all the buckets
   
   for (d = 0; d <= NCOLUMNS; d++) {
      int   idx;
      idx = buckets[d];
      buckets[d] = 0;
      while (idx) {
         int   i, c, next_d, next_idx;
         i = vectors[idx].i_diff;
         c = vectors[idx].c;
         if (sieve_array[d][c] == GOOD && c != 0 && d != 0) {
            int   index;
            index = alg_hash (c,d);
            alg_count = alg_faclist[index][0];
            alg_faclist[index][1] = c;
            alg_faclist[index][2] = d;
            alg_faclist[index][3+alg_count] = alg_fbase[base_i+i];
            alg_faclist[index][0]++;
         }
         c += saved_a1_b1[i].a1;
         next_d = d + saved_a1_b1[i].b1;
         next_idx = next_bucket_ptrs[idx];
         if (c >= 0 && c <= MAX_SIEVE && next_d <= NCOLUMNS) {
            vectors[idx].c = c;
            next_bucket_ptrs[idx] = buckets[next_d];
            buckets[next_d] = idx;
         }
         idx = next_idx;
      }
   }

// Reset count of saved vectors

   num_vectors = 0;
}   /* end of do_delayed_alg_resieve_vectors */

/************************************************************************/
/*                                                                      */
/* Routine to add one delayed vector update entry.  Returns TRUE if     */
/* successful, FALSE if caller should switch to storing individual pts  */
/*                                                                      */
/************************************************************************/

int __inline add_vector_entry (int d, int c, int b1, int a1, int i, char *logs, int alg)

{   /* start of add_vector_entry */

// If this vector is a complete miss, then return TRUE

   if (c < 0 || c > MAX_SIEVE || d > NCOLUMNS) return (TRUE);
//gwcount++;
//if (i>size_int_base-1000)
//printf ("D=%d,C=%d,i=%d,b1=%d,a1=%d\n",d,c,i,b1,a1);

// When primes get too large, switch to saving individual points

   if (i > NUM_VECTOR_UPDATE_PRIMES && num_vectors == 0) return (FALSE);

// optimization to try - quickly see if this vector generates a lot of hits
// send lots of hits to the bucket sort, send sparse ones to the array sort.
// May be a very bad idea as it may very well lead to lots of contention
// for the L2 cache.

// If we are starting a new set of vectors, then remember the base i value.
// From now on we can use the distance from this base.

   if (num_vectors == 0) base_i = i;
   i -= base_i;

// Add this vector to our arrays.  Note we do not use the first array entry
// so that zero can be used as an end-of-linked-list value.

// BUG - check that i does not exceed MAX_A1_B1
   saved_a1_b1[i].a1 = (short) a1;
   saved_a1_b1[i].b1 = (short) b1;
   num_vectors++;
   vectors[num_vectors].c = (short) c;
   vectors[num_vectors].i_diff = (short) i;
   next_bucket_ptrs[num_vectors] = buckets[d];
   buckets[d] = num_vectors;

// If we've filled up this array, go process it.

   if (num_vectors == MAX_VECTORS)
      if (logs != NULL) do_delayed_vectors (logs);
      else if (alg) do_delayed_alg_resieve_vectors ();
      else do_delayed_int_resieve_vectors ();
   return (TRUE);
}   /* end of add_vector_entry */


/****************************************************************************/
/*                                                                          */
/*                                                                          */
/* These routines accumulate points c,d to update sieve entries at a        */
/* later time.  When done accumulating points, they are sorted and the      */
/* sieve array is updated in a cache-friendly manner.                       */
/*                                                                          */
/* Defines to tune the accumulating / sorting process.                      */
/*                                                                          */
/* It is best if the sort_chunk_size fits in L2 cache.  It is probably even */
/* better if it is small enough to eliminated TLB cache misses.  It is also */
/* probably better if it is not a multiple of 4KB as the updating process   */
/* will read one cache line from each sort_chunk and we want to avoid L2    */
/* cache collisions during that process.                                    */
/*                                                                          */
/****************************************************************************/

// Structure describing one update entry

struct point_data 
   {
   short   d;
   short   c;
   unsigned long payload;
   };

// Point arrays / sorting global variables

struct point_data *point_array = NULL;
int    num_points = 0;
int    unsorted_count = 0;
unsigned short radix_sort_counts[NUM_SORT_KEYS] = {0};

   // Can we safely use an unsigned char array instead?  Can we easily
   // detect overflows and do a quicksort instead?

double   sort_time = 0.0;
double   add_point_time = 0.0;

/************************************************************************/
/*                                                                      */
/* Routine to init the delayed vector and point sieve updating code     */
/*                                                                      */
/************************************************************************/

void delayed_init ()
{   /* start of delayed_int */
   if (point_array == NULL) {
      // optimization to do - alloc on 8 byte boundary....
      point_array = (struct point_data *) malloc (MEM_FOR_POINT_UPDATES);
      if (point_array == NULL) {
         printf ("Unable to alloc sort area\n");
         exit(0);
      }
   }
   num_points = 0;
   unsorted_count = 0;
   num_vectors = 0;
gwcount = 0;
gwcount2 = 0;
}   /* end of delayed_init */

/************************************************************************/
/*                                                                      */
/*   QUICKSORT                                                          */
/*                                                                      */
/************************************************************************/

#define qswap(x,y) { struct point_data qtmp; qtmp = x; x = y; y = qtmp; }

void quicksort (struct point_data arr[], int beg, int end) 

{   /* start of quicksort */
   if (end > beg + 1) {
      int pivot = arr[beg].d, l = beg + 1, r = end;
      while (l < r) {
         if (arr[l].d <= pivot) l++;
         else {
            --r;
            qswap (arr[l], arr[r]);
         }
      }
      --l;
      qswap (arr[l], arr[beg]);
      quicksort (arr, beg, l);
      quicksort (arr, r, end);
   }
}   /* end of quicksort */

/************************************************************************/
/*                                                                      */
/* Radix-like sort.  During construction of the point_array we count #  */
/* times each sort key occurs.This lets us do an one-pass, in-place sort*/
/*                                                                      */
/************************************************************************/


void radix_sort (struct point_data *sort_array)

{   /* start of radix_sort */
   int   i;
   struct point_data *write_ptr[NUM_SORT_KEYS];

// The first step is to form an array that tells us where to store each
// sort key.

   write_ptr[0] = sort_array;
   for (i = 1; i < NUM_SORT_KEYS; i++)
      write_ptr[i] = write_ptr[i-1] + radix_sort_counts[i-1];

// Now loop until the array is sorted

   for (i = 0; i < NUM_SORT_KEYS; ) {
      struct point_data tmp, tmp2;

// Skip this bucket if it is already sorted

      if (radix_sort_counts[i] == 0) {
         i++;
         continue;
      }

// We've found the first entry that may not be in its correct position.
// Write it in its proper position, repeat using the ousted entry.
// When we finally write an ousted entry in the bucket we are working on
// then we can stop looping (and resume outer loop looking for next
// entry that may be out of position).

      tmp = *write_ptr[i];
      for ( ; ; ) {
         int   j;
         j = MAKE_SORT_KEY (tmp.d);
         if (j == i) {
            *(write_ptr[i])++ = tmp;
            break;
         }
         tmp2 = *write_ptr[j];
         *(write_ptr[j])++ = tmp;
         radix_sort_counts[j]--;

         j = MAKE_SORT_KEY (tmp2.d);
         if (j == i) {
            *(write_ptr[i])++ = tmp2;
            break;
         }
         tmp = *write_ptr[j];
         *(write_ptr[j])++ = tmp2;
         radix_sort_counts[j]--;
      }
      radix_sort_counts[i]--;
   }
}   /* end of radix_sort */


/************************************************************************/
/*                                                                      */
/* Sort one chunk of update data. Done every time we enough data        */
/* that we are worried it might start falling out of the L2 cache.      */
/*                                                                      */
/************************************************************************/

void do_sort ()
{   /* start of do_sort */
   struct point_data *sort_array;

   if (TIME_STATS) sort_time -= get_time ();

   sort_array = &point_array[num_points - unsorted_count];
//   quick_sort (sort_array, 0, unsorted_count-1);
   radix_sort (sort_array);
   unsorted_count = 0;

   if (TIME_STATS) sort_time += get_time();
}


// Routine to add one delayed sieve update entry.

void __inline add_point_entry (int d, int c, int i, char log)
{
   if (TIME_STATS && unsorted_count == 0) add_point_time -= get_time ();

   point_array[num_points].d = (short) d;
   point_array[num_points].c = (short) c;
   point_array[num_points].payload = (i << 8) | log;
   num_points++;

// Accumulate data to make sorting easier.

   radix_sort_counts[MAKE_SORT_KEY(d)]++;

// This attempt at prefetching a little ahead did not help
//   if ((num_points & 7) == 0)
//   _mm_prefetch (((char *) &point_array[num_points])+128, _MM_HINT_T1);

// Sort a chunk once it fills up

   unsorted_count++;
   if (unsorted_count == POINTS_PER_SORT_CHUNK) {
      if (TIME_STATS) add_point_time += get_time ();
      do_sort ();
   }

// BUG - detect overflow of point_array and either exit or do the updates
// and start another batch of delayed updates

}   /* end of do_sort */

/************************************************************************/
/*                                                                      */
/* This routine applies the delayed sieve updates                       */
/*                                                                      */
/************************************************************************/

void do_delayed_points ()

{   /* start of do_delayed_points */
   int   i, j, num_chunks;
   double   itime;
   struct point_data *read_ptrs[NUM_SORT_CHUNKS];

// Return if no delayed updates have been accumulated

   if (num_points == 0) return;

// Print some useful data
/*
   printf ("Doing %d delayed points: %dMB of %dMB used.\n", num_points,
      ((num_points * sizeof (struct point_data)) >> 20) + 1,
      MEM_FOR_POINT_UPDATES >> 20); */

// Sort the last chunk -- we never finished filling it up.

   if (unsorted_count) {
      if (TIME_STATS) add_point_time += get_time ();
      do_sort ();
   }

// Determine the number of sorted chunks and create pointers into each

   if (TIME_STATS) itime = get_time ();
   num_chunks = (num_points + POINTS_PER_SORT_CHUNK - 1) /
                  POINTS_PER_SORT_CHUNK;
   for (i = 0; i < num_chunks; i++)
      read_ptrs[i] = point_array + i * POINTS_PER_SORT_CHUNK;

// Create a dummy last entry so that when the last chunk is completely
// processed it runs into an entry that doesn't match any of the sort keys.
// NOTE: we assume that the first entry in each sort chunk is less than
// the last entry in the previous sort chunk.  We could eliminate this
// restriction by creating dummy entries between chunks or by expanding the
// if statement below to make sure we process no more than
// POINTS_PER_SORT_CHUNK for each read_ptrs.

   point_array[num_points].d = -1;

// Apply the sieve updates, one sort key at a time.
   
   for (i = 0; i <= NUM_SORT_KEYS; i++) {

// Work on each sort chunk processing all the matching sort keys

      for (j = 0; j < num_chunks; j++) {
         struct point_data *u;
         for (u = read_ptrs[j]; MAKE_SORT_KEY (u->d) == i; u++)
            sieve_array[u->d][u->c] += (char) u->payload;
         read_ptrs[j] = u;
      // optimization to try - we could throw in a prefetch here to
      // grab one or more cache lines in this chunk for use
      // on the next sort key.
      }
   }

// All done
/*
   if (TIME_STATS) {
      itime = (get_time() - itime);
      (void) printf ("Point array build time: %lf\n",add_point_time/(clockrate*1000.0));
      (void) printf ("Point sort time: %lf\n",sort_time/(clockrate*1000.0));
      (void) printf ("Point update time: %lf\n",itime/(clockrate*1000.0));
      add_point_time = 0.0;
      sort_time = 0.0;
   } */
}   /* end of do_delayed_points */

/************************************************************************/
/*                                                                      */
/* For int resieving, this routine applies the delayed sieve updates    */ 
/* saving the i values for sieve entries marked GOOD.                   */
/*                                                                      */
/************************************************************************/

void do_delayed_int_resieve_points ()
{
   int   i, j, num_chunks;
   double   itime;
   struct point_data *read_ptrs[NUM_SORT_CHUNKS];

// Return if no delayed updates have been accumulated

   if (num_points == 0) return;

// Determine the number of sorted chunks and create pointers into each

   if (TIME_STATS) itime = get_time ();
   num_chunks = (num_points + POINTS_PER_SORT_CHUNK - 1) /
                  POINTS_PER_SORT_CHUNK;
   for (i = 0; i < num_chunks; i++)
      read_ptrs[i] = point_array + i * POINTS_PER_SORT_CHUNK;

// Apply the sieve updates, one sort key at a time.
   
   for (i = 0; i <= NUM_SORT_KEYS; i++) {

// Work on each sort chunk processing all the matching sort keys

      for (j = 0; j < num_chunks; j++) {
         struct point_data *u;
         for (u = read_ptrs[j]; MAKE_SORT_KEY (u->d) == i; u++) {
            int   index;
            if (sieve_array[u->d][u->c] != GOOD) continue;
            index = int_hash (u->c,u->d);
            int_count = int_faclist[index][0];
            int_faclist[index][1] = u->c;
            int_faclist[index][2] = u->d;
            int_faclist[index][3+int_count] = int_fbase[u->payload >> 8];
            int_faclist[index][0]++;
         }
         read_ptrs[j] = u;
      // optimization to try - we could throw in a prefetch here to
      // grab one or more cache lines in this chunk for use
      // on the next sort key.
      }
   }

// All done
/*
   if (TIME_STATS) {
      itime = (get_time() - itime);
      (void) printf ("Point int resieve time: %lf\n",itime/(clockrate*1000.0));
   } */
}

/************************************************************************/
/*                                                                      */
/* For alg resieving (step 1 of 2).  This routine examines the delayed  */
/* sieve  updates saving the entries that passed the algebraic sieve    */
/* (not marked -BIGNUM in the sieve) .                                  */
/*                                                                      */
/************************************************************************/

void do_delayed_alg_resieve_points ()
{
   int   i, j, num_chunks;
   double   itime;
   struct point_data *read_ptrs[NUM_SORT_CHUNKS];
   struct point_data *write_ptr;

// Return if no delayed updates have been accumulated

   if (num_points == 0) return;

// Determine the number of sorted chunks and create pointers into each

   if (TIME_STATS) itime = get_time ();
   num_chunks = (num_points + POINTS_PER_SORT_CHUNK - 1) /
                  POINTS_PER_SORT_CHUNK;
   for (i = 0; i < num_chunks; i++)
      read_ptrs[i] = point_array + i * POINTS_PER_SORT_CHUNK;

// Test the sieve updates, one sort key at a time.

   write_ptr = LAST_POINT_PTR;
   for (i = 0; i <= NUM_SORT_KEYS; i++) {

// Work on each sort chunk processing all the matching sort keys

      for (j = 0; j < num_chunks; j++) {
         struct point_data *u;
         for (u = read_ptrs[j]; MAKE_SORT_KEY (u->d) == i; u++)
            if (sieve_array[u->d][u->c] !=
                  (unsigned char) -BIGNUM)
               *write_ptr-- = *u;
         read_ptrs[j] = u;
      // optimization to try - we could throw in a prefetch here to
      // grab one or more cache lines in this chunk for use
      // on the next sort key.
      }
   }

// Output a dummy ending record

   write_ptr->d = 0;
   write_ptr->c = 0;
   write_ptr--;
      
// All done
/*
   if (TIME_STATS) {
      int   count;
      count = LAST_POINT_PTR - write_ptr;
      itime = (get_time() - itime);
      printf ("Point alg resieve time (additional %dMB used): %lf\n",
         ((count * sizeof (struct point_data)) >> 20) + 1,
         itime/(clockrate*1000.0));
   } */
}

/************************************************************************/
/*  For alg resieving (step 2 of 2).  This routine looks at the saved   */
/*  updates  from step 1 that passed the algebraic sieve and saves the  */
/*  values that also  passed the integer sieve (marked GOOD in the      */
/*  sieve array).                                                       */
/************************************************************************/

void do_twice_delayed_alg_resieve_points ()
{
   struct point_data *read_ptr;
   int   index;
   double   itime;

   if (TIME_STATS) itime = get_time ();

// Test the saved sieve updates
   
   for (read_ptr = LAST_POINT_PTR; ; read_ptr--) {
      if (read_ptr->d == 0 && read_ptr->c == 0) break;
      if (sieve_array[read_ptr->d][read_ptr->c] != GOOD) continue;
      index = alg_hash (read_ptr->c,read_ptr->d);
      alg_count = alg_faclist[index][0];
      alg_faclist[index][1] = read_ptr->c;
      alg_faclist[index][2] = read_ptr->d;
      alg_faclist[index][3+alg_count] =
         alg_fbase[read_ptr->payload >> 8];
      alg_faclist[index][0]++;
   }

/*
   if (TIME_STATS) {
      itime = (get_time() - itime);
      (void) printf ("Point alg resieve2 time: %lf\n",itime/(clockrate*1000.0));
   }
*/
}

/************************************************************************/
/*                                                                      */
/*   This routine sieves the rational large primes by vector sieving    */
/*   N.B. try unrolling inner loop                                      */
/*                                                                      */
/************************************************************************/

void vector_int_sieve(sign)
int sign;

{   /* start of vector_int_sieve */
int i,c0,d0,c1,d1,c,d,e_min, e_max;
int e,f;                                   /*  Pollard's notation  */
int a1, b1, a0, b0;
int f_lower,f_upper;
int upper_hits,lower_hits;
double itime;
int maxf;

delayed_init ();                        // Init delayed sieve updating code

/*   For each prime >= MAX_SIEVE we need to sieve by vectors. int_v1[i]   */
/*  and int_v2[i] are reduced bases for the ith prime.  We need to find   */
/*  all  (e,f) pairs such that e*V1[i] + f*V2[i]  < (MAX_SIEVE,NCOLUMNS)  */
/*  i.e. inside the reduce q-lattice. sign indicates if we are sieving    */
/*  the q-lattice (c,d) with c > 0  or c < 0.                             */
/*  Let v1 = (a0, b0),  v2 = (a1,b1).  We require                         */
/*      -MAX_SIEVE < e*a0 + f*a1 < MAX_SIEVE                              */
/*      0 <  e*b0 + f*b1 < NCOLUMNS            to get bounds on e,f       */

if (DEBUG_INT_VECTOR) 
		{ (void) printf("Vector INT sieve %d\n",sign); fflush(stdout); }
if (TIME_STATS) itime = get_time();

if (sign == 1) {

upper_hits = 0;
for (i=int_linesieve;  i<size_int_base; i++)
   {
   if (DEBUG_INT_VECTOR) 
		{ 
		(void) printf("int vector[%d] = ",int_fbase[i]); 
		show_int_vector(i); 
		}

   if (int_v1[i][0] == 0) continue;      /* no hits in lattice      */

   if (COMPUTE_AREA) int_upper_area[i] = 0;
   prepare_bounds(int_v1[i], int_v2[i], sign);
   lc1 += .50;							/* to compute ceiling		*/
   lc2 += .50;

   if (SHOW_INT_REGION)
      {
      (void) printf("upper int region %d\n",i);
      show_region(int_fbase[i]);
      }

   a0 = int_v1[i][0];
   b0 = int_v1[i][1];
   a1 = int_v2[i][0];
   b1 = int_v2[i][1];
   e_min = (int)emin;
   e_max = (int)emax;

   if (CHECK_LATTICE) if (i < 50000) check_lattice("INT", i, int_v1[i], int_v2[i]);
   if (DEBUG_INT_VECTOR)
      {
      (void) printf("emin, emax = %d %d\n",e_min,e_max);
      if (emin > emax) {(void) printf("FATAL ERROR\n"); exit(0);}
      }                           
   c0 = e_min * a0;
   d0 = e_min * b0;

   for (e=e_min; e<=e_max; e++)
      { 
 //     f_lower = (int)lower_bound((double)e);       older code

	  if (e < lower_midpt) f_lower = round_double(lm1 * (double)e + lc1);
	  else                 f_lower = round_double(lm2 * (double)e + lc2);

      c1 = f_lower * a1;
      d1 = f_lower * b1; 
      c = c0 + c1;
      d = d0 + d1;
      if (c == 0 && d == 0) c = a1, d = b1;

// Try to add this entire vector to the delayed update list

      if (! add_vector_entry (d, c, b1, a1, i, int_logp, 0)) {

// We couldn't add the entire vector, add each point instead

      if (SHOW_INT_LATTICE) 
		 (void) printf("e, lower,upper,c,d = %d %d %d %d %d\n",e,f_lower,f_upper,c,d);

      while (c >= 0 && c <= MAX_SIEVE && d <= NCOLUMNS)
         {
         if (SHOW_INT_LATTICE) 
			(void) printf("i, e, f, c, d = %d %d %d %d %d\n",i, e, f, c, d);

         if (CHECK_BOUNDS)
            {
            if (c < 0 || c > MAX_SIEVE) 
					{ (void) printf("BOUND c at (%d,%d), (%d,%d)\n",c,d,e,f); exit(0); }
            if (d < 0 || d > NCOLUMNS)  
					{ (void) printf("BOUND d at (%d,%d), (%d,%d)\n",c,d,e,f); exit(0); }
            }                                          
//       sieve_array[d][c] += int_logp[i];

         add_point_entry (d, c, i, int_logp[i]);

         if (COMPUTE_AREA) int_upper_area[i]++;
//       upper_hits++;
         c += a1;      
         d += b1;      
         }
      }
         
      c0 += a0; 
      d0 += b0;
      } 
   if (SHOW_F) if (i < 50000) (void) printf("max f for %d = %d\n",int_fbase[i], maxf);
   }

}  /* end of sign = 1 */

else {

lower_hits = 0;
for (i=int_linesieve;  i<size_int_base; i++)
   {
   if (DEBUG_INT_VECTOR) 
		{ (void) printf("int vector[%d] = ", int_fbase[i]); show_int_vector(i); }

   if (int_v1[i][0] == 0) continue;      /* no hits in lattice      */
   
   if (COMPUTE_AREA) int_lower_area[i] = 0;

   prepare_bounds(int_v1[i], int_v2[i], sign);
   lc1 += .50;
   lc2 += .50;

   if (SHOW_INT_REGION)
      {
      (void) printf("lower int region %d\n",i);
      show_region(int_fbase[i]);
      }

   a0 = int_v1[i][0];
   b0 = int_v1[i][1];
   a1 = int_v2[i][0];
   b1 = int_v2[i][1];
   e_min = (int)emin;
   e_max = (int)emax;

   if (DEBUG_INT_VECTOR) (void) printf("emin, emax = %d %d\n",e_min,e_max);
                           
   c0 = e_min * a0;
   d0 = e_min * b0;
   for (e=e_min; e<=e_max; e++)
      {    
//      f_lower = (int)lower_bound((double)e);  older code

	  if (e < lower_midpt) f_lower = round_double(lm1 * (double)e + lc1);
	  else                 f_lower = round_double(lm2 * (double)e + lc2);

      c1 = f_lower * a1;
      d1 = f_lower * b1; 
      c = c0 + c1;
      d = d0 + d1;
      if (c == 0 && d == 0) c = a1, d = b1;

// Try to add this entire vector to the delayed update list

      if (! add_vector_entry (d, -c, b1, -a1, i, int_logp, 0)) {

// We couldn't add the entire vector, add each point instead

      if (SHOW_INT_LATTICE) 
		(void) printf("e, lower,upper,c,d = %d %d %d %d %d\n",e,f_lower,f_upper,c,d);

      while (c <= 0 && c >= -MAX_SIEVE && d <= NCOLUMNS)
         {
         if (SHOW_INT_LATTICE) (void) printf("f, c, d = %d %d %d\n",f,c,d);
         if (CHECK_BOUNDS)
            {
            if (c < -MAX_SIEVE || c > 0) 
				{ (void) printf("BOUND c at (%d,%d), (%d,%d)\n",c,d,e,f); exit(0); }
            if (d < 0 || d > NCOLUMNS)   
				{ (void) printf("BOUND d at (%d,%d), (%d,%d)\n",c,d,e,f); exit(0); }
            }                                          
//         sieve_array[d][-c] += int_logp[i]; 
//         lower_hits++;
         add_point_entry (d, -c, i, int_logp[i]);
         if (COMPUTE_AREA) int_lower_area[i]++;
         c += a1;      
         d += b1;      
         }

      }

      c0 += a0; 
      d0 += b0;
      }
   }
}   /* end of sign = -1 */

do_delayed_vectors (int_logp);
do_delayed_points ();

if (NOTIFY_HITS)
   {
   if (sign == 1) (void) printf("int upper hits = %d\n",upper_hits);
   else (void) printf("int lower hits = %d\n",lower_hits);
   }
if (TIME_STATS) int_vector_time += (get_time() - itime);

}   /* end of vector_int_sieve */


/************************************************************************/
/*                                                                      */
/*   This routine sieves the algebraic large primes by vector sieving   */
/*                                                                      */
/************************************************************************/

void vector_alg_sieve(sign)
int sign;

{   /* start of vector_alg_sieves */
int i,c,d,c0,d0,c1,d1;
int e,f,e_min,e_max;                     /*  Pollard's notation  */
int a0, b0, a1, b1;
int f_lower,f_upper;
int upper_hits,lower_hits;
double atime;

delayed_init ();      // Init delayed sieve updating code

/*   For each prime >= MAX_SIEVE we need to sieve by vectors. alg_v1[i]   */
/*  and alg_v2[i] are reduced bases for the ith prime.  We need to find   */
/*  all  (e,f) pairs such that e*V1[i] + f*V2[i]  < (MAX_SIEVE,NCOLUMNS)  */
/*  i.e. inside the reduce q-lattice. sign indicates if we are sieving    */
/*  the q-lattice (c,d) with c > 0  or c < 0.                             */
/*  Let v1 = (a0, b0),  v2 = (a1,b1).  We require                         */
/*      -MAX_SIEVE < e*a0 + f*a1 < MAX_SIEVE                              */
/*      0 <  e*b0 + f*b1 < NCOLUMNS         to get bounds on e,f          */

if (TIME_STATS) atime = get_time();
if (DEBUG_ALG_VECTOR) 
		{ (void) printf("Vector ALG sieve %d\n",sign); fflush(stdout); }

if (sign == 1) {

upper_hits = 0;
for (i=alg_linesieve;  i<size_alg_base; i++)
   {
   if (DEBUG_ALG_VECTOR) 
      { 
      (void) printf("alg vector[%d] = ", alg_fbase[i]);
      show_alg_vector(i);
      fflush(stdout);
      }
   if (alg_v1[i][0] == 0) continue;      /* no hits in lattice      */

   if (COMPUTE_AREA) alg_upper_area[i] = 0;
   prepare_bounds(alg_v1[i], alg_v2[i], sign);
   lc1 += .50;                           /* for ceiling function    */
   lc2 += .50;

   if (SHOW_ALG_REGION)
      {
      (void) printf("alg region %d\n",i);
      show_region(alg_fbase[i]);
      } 
   a0 = alg_v1[i][0];
   b0 = alg_v1[i][1];
   a1 = alg_v2[i][0];
   b1 = alg_v2[i][1];
   e_min = (int)emin;
   e_max = (int)emax;

   if (DEBUG_ALG_VECTOR)

      {
      (void) printf("emin, emax, midpt = %d %d %lf\n",e_min,e_max,lower_midpt);
      if (emin > emax) {(void) printf("FATAL ERROR\n"); exit(0);}
      }
                           
   c0 = e_min * a0;
   d0 = e_min * b0;
   for (e=e_min; e<=e_max; e++)
      {       
//      f_lower = (int)lower_bound((double)e);   older code

	  if (e < lower_midpt) f_lower = round_double(lm1 * (double)e + lc1);
	  else                 f_lower = round_double(lm2 * (double)e + lc2);

      c1 = f_lower * a1;
      d1 = f_lower * b1; 
      c = c0 + c1;
      d = d0 + d1;
      if (c == 0 && d == 0) c = a1, d = b1;

// Try to add this entire vector to the delayed update list

      if (! add_vector_entry (d, c, b1, a1, i, alg_logp, 1)) {

// We couldn't add the entire vector, add each point instead

      if (SHOW_ALG_LATTICE) 
         (void) printf("e, lower,upper,c,d = %d %d %d %d %d\n",e,f_lower,f_upper,c,d);

      while (c >= 0 && c <= MAX_SIEVE && d <= NCOLUMNS)
         {   
         if (SHOW_ALG_LATTICE) (void) printf("f, c, d = %d %d %d\n",f,c,d);
         if (CHECK_BOUNDS)
            {
            if (c < 0 || c > MAX_SIEVE) 
				{ (void) printf("BOUND c at (%d,%d), (%d,%d)\n",c,d,e,f); exit(0); }
            if (d < 0 || d > NCOLUMNS)  
				{ (void) printf("BOUND d at (%d,%d), (%d,%d)\n",c,d,e,f); exit(0); }
            }                                          
//       sieve_array[d][c] += alg_logp[i]; 
         add_point_entry (d, c, i, alg_logp[i]);
//       upper_hits++;
         if (COMPUTE_AREA) alg_upper_area[i]++;
         c += a1;      
         d += b1;      
         }
      }
      c0 += a0; 
      d0 += b0;
      } 
   }

}   /* end of sign = 1 */

else {

lower_hits = 0;
for (i=alg_linesieve;  i<size_alg_base; i++)
   {
   if (DEBUG_ALG_VECTOR) 
		{ (void) printf("vector[%d] = ", alg_fbase[i]); show_alg_vector(i); }

   if (alg_v1[i][0] == 0) continue;      /* no hits in lattice      */

   if (COMPUTE_AREA) alg_lower_area[i] = 0;
   prepare_bounds(alg_v1[i], alg_v2[i], sign);
   lc1 += .50;
   lc2 += .50;

   if (SHOW_ALG_REGION) 
      {
      (void) printf("alg region %d\n",i);
      show_region(alg_fbase[i]);
      }

   a0 = alg_v1[i][0];
   b0 = alg_v1[i][1];
   a1 = alg_v2[i][0];
   b1 = alg_v2[i][1];
   e_min = (int)emin;
   e_max = (int)emax;

   if (DEBUG_ALG_VECTOR) (void) printf("emin, emax = %d %d\n",e_min,e_max);
                           
   c0 = e_min * a0;
   d0 = e_min * b0;
   for (e=e_min; e<=e_max; e++)
      {       
	  if (e < lower_midpt) f_lower = round_double(lm1 * (double)e + lc1);
	  else                 f_lower = round_double(lm2 * (double)e + lc2);

      c1 = f_lower * a1;
      d1 = f_lower * b1; 
      c = c0 + c1;
      d = d0 + d1;
      if (c == 0 && d == 0) c = a1, d = b1;

// Try to add this entire vector to the delayed update list

      if (! add_vector_entry (d, -c, b1, -a1, i, alg_logp, 1)) {

// We couldn't add the entire vector, add each point instead

      if (SHOW_ALG_LATTICE) 
			(void) printf("e, lower,upper,c,d = %d %d %d %d %d\n",e,f_lower,f_upper,c,d);

      while (c <= 0 && c >= -MAX_SIEVE && d <= NCOLUMNS)
         {   
         if (SHOW_ALG_LATTICE) (void) printf("f, c, d = %d %d %d\n",f,c,d);
         if (CHECK_BOUNDS)
            {
            if (c < -MAX_SIEVE || c > 0) 
				{ (void) printf("BOUND c at (%d,%d), (%d,%d)\n",c,d,e,f); exit(0); }
            if (d < 0 || d > NCOLUMNS)   
				{ (void) printf("BOUND d at (%d,%d), (%d,%d)\n",c,d,e,f); exit(0); }
            }                                          
//       sieve_array[d][-c] += alg_logp[i]; 
         add_point_entry (d, -c, i, alg_logp[i]);
//       lower_hits++;
         if (COMPUTE_AREA) alg_lower_area[i]++;
         c += a1;
         d += b1;
         }
      }
      c0 += a0; 
      d0 += b0;
      } 
   }

}   /* end of sign = -1 */

do_delayed_vectors (alg_logp);
do_delayed_points ();

if (NOTIFY_HITS)
   {
   if (sign == 1) (void) printf("alg upper hits = %d\n",upper_hits);
   else (void) printf("alg lower hits = %d\n",lower_hits);
   }

if (TIME_STATS) alg_vector_time += (get_time() - atime); 

}   /* end of vector_alg_sieves */
 
/************************************************************************/
/*                                                                      */
/*   This routine computes the initial starting points for sieving      */
/*  based upon the current special q  and V1, V2                        */
/*                                                                      */
/************************************************************************/
 
void find_start_points()

 
{   /* start of find_start_points */
int i,numer, denom, inverse;
int prime,root;
double dtime;

/*   startpt = -(a1 + b1 root) * (a0 + b0 root)^-1  mod p              */
/*   note:  if denom = 0,  then there is no root;  do not sieve.       */
/*   instead, this prime takes hits only when p | f. Deal with that    */
/*   during sieve inititalization. Add log(p) to EVERY location where  */
/*    p | the column.                                                  */
/*   if numer = 0, then the sub-lattice is (p,0), (0,1) with root = 0  */
/*   if this happens to a vector sieve prime, then it takes no hits in */
/*   the lattice because (p > MAX_SIEVE)                               */               

if (TIME_STATS) dtime = get_time();

proj_int_primes[0] = 0;
proj_alg_primes[0] = 0;       /* count of projective primes            */

/*   Do the roots for the algebraic and rational small primes. These   */
/*  are special cases, since we deal with prime powers                 */

if (DEBUG_FIND) (void) printf("DOING FIND_START_POINTS\n"); fflush(stdout);
if (DEBUG_FIND) (void) printf("lat: (%d,%d), (%d,%d)\n",v1[0],v1[1],v2[0],v2[1]);

/*   This loop takes negligible time; NUMSMALLP is tiny                */

for (i=0; i<NUMSMALLP; i++)
   {
   root = v2[1] * int_roots[i] + v2[0];
   numer = (int)(root % smallp[i]);
   if (numer < 0) numer += smallp[i];

   root = v1[1] * int_roots[i] + v1[0];
   denom = (int)(root % smallp[i]);
   if (denom < 0) denom += smallp[i]; 
   if (denom == 0)              /* projective root                     */
      {
      if (DEBUG_FIND) 
			(void) printf("Int smallp proj root at %d %d\n",i, int_fbase[i]);

      proj_int_primes[0]++;                  /* make a list of these   */
      proj_int_primes[proj_int_primes[0]] = i;
      special_q_int_roots[i] = 0;
      int_startpts[i] = 0;
      proj_int_logs[proj_int_primes[0]] = int_logp[i];   /* save log   */
      int_logp[i] = 0;                       /* set to 0 temporarily   */
      }                                      /* set it back later      */
   else
      { 
      inverse = single_modinv(denom, smallp[i]); 
      special_q_int_roots[i] = mod_mult_asm(numer, inverse, smallp[i]);
      int_startpts[i] = smallp[i] -  special_q_int_roots[i];
      }   
   }

/* Now do the primes that we line sieve                                */

if (DEBUG_FIND) (void) printf("Finding startpts for Int linesieve primes\n");

/*   This loop does not take much time; int_linesieve is small         */

for (; i<int_linesieve; i++)
   {
   root = v2[1] * int_roots[i] + v2[0];         /* can't overflow      */
   numer = (int)(root % int_fbase[i]);
   if (numer < 0) numer += int_fbase[i];

   root = v1[1] * int_roots[i] + v1[0];
   denom = (int)(root % int_fbase[i]);  
   if (denom < 0) denom += int_fbase[i];
   if (denom == 0)              /* projective root                     */
      {
      if (DEBUG_FIND) 
		(void) printf("int linesieve proj root at %d %d\n",i, int_fbase[i]);

      proj_int_primes[0]++;
      proj_int_primes[proj_int_primes[0]] = i;
      special_q_int_roots[i] = 0;
      int_startpts[i] = 0;
      proj_int_logs[proj_int_primes[0]] = int_logp[i];   /* save log   */
      int_logp[i] = 0;                       /* set to 0 temporarily   */
      }
   else
      { 
      inverse = single_modinv(denom, int_fbase[i]); 
      special_q_int_roots[i] = mod_mult_asm(numer, inverse, int_fbase[i]);
      int_startpts[i] = int_fbase[i] -  special_q_int_roots[i];
      }
   if (DEBUG_FIND) 
		(void) printf("Int find: p, num, den, inv, start: %d %d %d %d %d\n",
						int_fbase[i], numer,denom,inverse, int_startpts[i]);
   }

/*   Finally, the primes to be sieved by vectors                       */
/*   This loop takes lots of time.                                     */
   
if (DEBUG_FIND) (void) printf("Doing Int vector primes\n");

for (; i<size_int_base; i++)
   { 
   numer = modmultadd_asm(v2[1], int_roots[i], v2[0], int_fbase[i]);  
   denom = modmultadd_asm(v1[1], int_roots[i], v1[0], int_fbase[i]);  
   if (denom == 0 || numer == 0)           /* projective root          */
      {                                    /* takes no hits            */
      int_v1[i][0] = 0;
      }
   else
      { 
      inverse = single_modinv(denom, int_fbase[i]); 
      root = int_fbase[i] - mod_mult_asm(numer, inverse, int_fbase[i]);
      reduce_lattice(int_fbase[i], root, int_v1[i], int_v2[i]);
//      newest_latred_exp(int_fbase[i], root, int_v1[i], int_v2[i]);
//      latred(int_fbase[i], root, int_v1[i], int_v2[i]);
      }
   if (DEBUG_FIND) 
      {
      (void) printf("Int vector: p, num, den, inv start: %d  %d %d %d %d  Vector: ", 
                  int_fbase[i],numer,denom,inverse, root);
      show_int_vector(i);
      }
   }

/*   Now do the algebraic small prime powers                            */
/*   If the prime has  root mod a power > 1, then use the power in place*/
/*  of the prime. Else, just use the prime                              */

if (DEBUG_FIND) (void) printf("Doing ALG small prime startpts\n");

i = 0;
while (alg_fbase[i]  <= SMALLP)
   {
   prime = alg_fbase[i];
   if (alg_root_ok[i] == TRUE)      /* root mod prime power exists      */
      {
      root = v2[1] * alg_prime_power_root[i] + v2[0];
      numer = (int) (root % alg_smallp[prime]);
      if (numer < 0) numer += alg_smallp[prime];

      root = v1[1] * alg_prime_power_root[i] + v1[0];
      denom = (int)(root % alg_smallp[prime]);
      if (denom < 0) denom += alg_smallp[prime];
      if (denom == 0)                 /* projective root                */
         {
         if (DEBUG_FIND) (void) printf("alg proj root at %d\n",alg_fbase[i]);
         proj_alg_primes[0]++;
         proj_alg_primes[proj_alg_primes[0]] = i;
         special_q_alg_roots[i] = 0;
         alg_startpts[i] = 0;
         proj_alg_logs[proj_alg_primes[0]] = alg_logp[i];   /* save log */
         alg_logp[i] = 0;                       /* set to 0 temporarily */
         }
      else
         {
         inverse = single_modinv(denom, alg_smallp[prime]); 
         special_q_alg_roots[i] = mod_mult_asm(numer, inverse, alg_smallp[prime]);
         alg_startpts[i] = alg_smallp[prime] - special_q_alg_roots[i];
         }
      }
   else
      {
      root = v2[1] * alg_roots[i] + v2[0];
      numer = (int) (root % alg_fbase[i]);
      if (numer < 0) numer += alg_fbase[i];

      root = v1[1] * alg_roots[i] + v1[0];
      denom = (int)(root % alg_fbase[i]);
      if (denom < 0) denom += alg_fbase[i];
      if (denom == 0)                 /* projective root                */
         {
         if (DEBUG_FIND) (void) printf("alg proj root at %d %d\n",i, alg_fbase[i]);
         proj_alg_primes[0]++;
         proj_alg_primes[proj_alg_primes[0]] = i;
         special_q_alg_roots[i] = 0;
         alg_startpts[i] = 0;
         proj_alg_logs[proj_alg_primes[0]] = alg_logp[i];   /* save log */
         alg_logp[i] = 0;                       /* set to 0 temporarily */
         }
      else
         {
         inverse = single_modinv(denom, alg_fbase[i]);
         special_q_alg_roots[i] = mod_mult_asm(numer, inverse, alg_fbase[i]);
         alg_startpts[i] = alg_fbase[i] - special_q_alg_roots[i];
         }
      }
   if (DEBUG_FIND) 
      {
       (void) printf("ALG smallp start: p, num, den, inv, start: %d %d %d %d %d\n\n", 
					prime, numer,denom,inverse, alg_startpts[i]);
      fflush(stdout);
      }
   i++;
   }
 
/*   Now do the line sieve primes                                       */

for (; i<alg_linesieve; i++)
   {
   root = v2[1] * alg_roots[i] + v2[0];              /* can't overflow  */
   numer = (int) (root % alg_fbase[i]);
   if (numer < 0) numer += alg_fbase[i];

   root = v1[1] * alg_roots[i] + v1[0];
   denom = (int)(root % alg_fbase[i]);
   if (denom < 0) denom += alg_fbase[i];
   if (denom == 0)               /* projective root                     */
      {
      if (DEBUG_FIND) (void) printf("alg proj root at %d\n",alg_fbase[i]);
      proj_alg_primes[0]++;
      proj_alg_primes[proj_alg_primes[0]] = i;
      special_q_alg_roots[i] = 0;
      alg_startpts[i] = 0;
      proj_alg_logs[proj_alg_primes[0]] = alg_logp[i];    /* save log   */
      alg_logp[i] = 0;                        /* set to 0 temporarily   */
      }
   else
      {
      inverse = single_modinv(denom, alg_fbase[i]);
      special_q_alg_roots[i] = mod_mult_asm(numer, inverse, alg_fbase[i]);
      alg_startpts[i] = alg_fbase[i] - special_q_alg_roots[i];
      }
   if (DEBUG_FIND) printf("Alg line p, num, den, inv, start: %d %d %d %d %d\n", alg_fbase[i], 
                     numer,denom,inverse, alg_startpts[i]);
   }

/*   And finally, the vector sieve primes                               */

if (DEBUG_FIND) (void) printf("Doing Alg vector primes\n");

for (; i<size_alg_base; i++)
   {
   numer = modmultadd_asm(v2[1], alg_roots[i], v2[0], alg_fbase[i]);  
   denom = modmultadd_asm(v1[1], alg_roots[i], v1[0], alg_fbase[i]);
   if (denom == 0 || numer == 0)            /* projective root          */
      {
      alg_v1[i][0] = 0;
      }
   else
      {
      inverse = single_modinv(denom, alg_fbase[i]);
      root = alg_fbase[i] - mod_mult_asm(numer, inverse, alg_fbase[i]);
      reduce_lattice(alg_fbase[i], root, alg_v1[i], alg_v2[i]);
//      newest_latred_exp(alg_fbase[i], root, alg_v1[i], alg_v2[i]);
//      latred(alg_fbase[i], root, alg_v1[i], alg_v2[i]);
      }
   if (DEBUG_FIND) printf("Alg vector: p, num, den, inv, start: %d %d %d %d %d\n",alg_fbase[i], 
                     numer,denom,inverse, root);
   }

if (TIME_STATS) find_startpt_time += (get_time() - dtime);

if (DEBUG_FIND) check_start_points();

}   /* end of find_start_points */

 

/************************************************************************/
/*                                                                      */
/*   This routine updates the starting points as column increases by    */
/*   1. i.e. b --> b+1 requires updating the initial sieve start        */
/*   points.   This requires just a modular add. Done only for linesieve*/
/*   primes                                                             */
/************************************************************************/
 
void update_int_start_points(sign)
int sign;

{   /* start of update_int_start_points */

register int i,j;

if (sign == 1) {

for (i=0; i<NUMSMALLP; i++)      /* small prime powers first   */
   {
   int_startpts[i] -= special_q_smalliroot[i];
   if (int_startpts[i] < 0) int_startpts[i] += smallp[i];
   }
 
for (; i<int_linesieve-4; i+=4)
   {
   int_startpts[i] -= special_q_int_roots[i];
   if (int_startpts[i] < 0) int_startpts[i] += int_fbase[i];

   j = i+1;
   int_startpts[j] -= special_q_int_roots[j];
   if (int_startpts[j] < 0) int_startpts[j] += int_fbase[j];

   j = i+2;
   int_startpts[j] -= special_q_int_roots[j];
   if (int_startpts[j] < 0) int_startpts[j] += int_fbase[j];

   j = i+3;
   int_startpts[j] -= special_q_int_roots[j];
   if (int_startpts[j] < 0) int_startpts[j] += int_fbase[j];
   }
 
for (; i<int_linesieve; i++)
   {
   int_startpts[i] -= special_q_int_roots[i];
   if (int_startpts[i] < 0) int_startpts[i] += int_fbase[i];
   }

}   /* end of sign = 1 (upper half plane*/

else {

for (i=0; i<NUMSMALLP; i++)      /* small prime powers first   */
   {
   int_startpts[i] += special_q_smalliroot[i];
   if (int_startpts[i] > smallp[i]) int_startpts[i] -= smallp[i];
   }
 
for (; i<int_linesieve-4; i+=4)
   {
   int_startpts[i] += special_q_int_roots[i];
   if (int_startpts[i] > int_fbase[i]) int_startpts[i] -= int_fbase[i];

   j = i+1;
   int_startpts[j] += special_q_int_roots[j];
   if (int_startpts[j] > int_fbase[i]) int_startpts[j] -= int_fbase[j];

   j = i+2;
   int_startpts[j] += special_q_int_roots[j];
   if (int_startpts[j] > int_fbase[i]) int_startpts[j] -= int_fbase[j];

   j = i+3;
   int_startpts[j] += special_q_int_roots[j];
   if (int_startpts[j] > int_fbase[i]) int_startpts[j] -= int_fbase[j];
   }
 
for (; i<int_linesieve; i++)
   {
   int_startpts[i] += special_q_int_roots[i];
   if (int_startpts[i] > int_fbase[i]) int_startpts[i] -= int_fbase[i];
   }

}   /* end of sign = -1 */

if (DEBUG_UPDATE)
   {
   (void) printf("int update\n");
   for (i=0; i<30; i++) (void) printf("%d ",int_startpts[i]);
   (void) printf("\n");
   }

}   /* end of update_int_start_points */


/************************************************************************/
/*                                                                      */
/*   This routine updates the starting points as column increases by    */
/*   1. i.e. b --> b+1 requires updating the initial sieve start        */
/*   points.   This requires just a modular add                         */
/*                                                                      */
/************************************************************************/
 
void update_alg_start_points(sign)
int sign;

{   /* start of update_alg_start_points */

register int i,j;

if (sign == 1) {

i = 0;
while (alg_fbase[i] <= SMALLP)
   {
   alg_startpts[i] -= special_q_alg_roots[i];
   if (alg_startpts[i] < 0) alg_startpts[i] += alg_smallp[alg_fbase[i]];
   i++;
   }

for (; i<alg_linesieve-4; i+=4)
   {
   alg_startpts[i] -= special_q_alg_roots[i];
   if (alg_startpts[i] < 0) alg_startpts[i] += alg_fbase[i];

   j = i+1;
   alg_startpts[j] -= special_q_alg_roots[j];
   if (alg_startpts[j] < 0) alg_startpts[j] += alg_fbase[j];

   j = i+2;
   alg_startpts[j] -= special_q_alg_roots[j];
   if (alg_startpts[j] < 0) alg_startpts[j] += alg_fbase[j];

   j = i+3;
   alg_startpts[j] -= special_q_alg_roots[j];
   if (alg_startpts[j] < 0) alg_startpts[j] += alg_fbase[j];
   }
 
for (; i<alg_linesieve; i++)
   {
   alg_startpts[i] -= special_q_alg_roots[i];
   if (alg_startpts[i] < 0) alg_startpts[i] += alg_fbase[i];
   }

}   /* end of sign = 1 */

else {

i = 0;
while (alg_fbase[i] <= SMALLP)
   {
   alg_startpts[i] += special_q_alg_roots[i];
   if (alg_startpts[i] > alg_smallp[alg_fbase[i]]) alg_startpts[i] -= alg_smallp[alg_fbase[i]];
   i++;
   }

for (; i<alg_linesieve-4; i+=4)
   {
   alg_startpts[i] += special_q_alg_roots[i];
   if (alg_startpts[i] > alg_fbase[i]) alg_startpts[i] -= alg_fbase[i];

   j = i+1;
   alg_startpts[j] += special_q_alg_roots[j];
   if (alg_startpts[j] > alg_fbase[i]) alg_startpts[j] -= alg_fbase[j];

   j = i+2;
   alg_startpts[j] += special_q_alg_roots[j];
   if (alg_startpts[j] > alg_fbase[i]) alg_startpts[j] -= alg_fbase[j];

   j = i+3;
   alg_startpts[j] += special_q_alg_roots[j];
   if (alg_startpts[j] > alg_fbase[i]) alg_startpts[j] -= alg_fbase[j];
   }
 
for (; i<alg_linesieve; i++)
   {
   alg_startpts[i] += special_q_alg_roots[i];
   if (alg_startpts[i] > alg_fbase[i]) alg_startpts[i] -= alg_fbase[i];
   }

}

if (DEBUG_UPDATE)
   {
   (void) printf("alg update\n");
   for (i=0; i<25; i++) (void) printf("%d ",alg_startpts[i]);
   (void) printf("\n");
   }

}   /* end of update_alg_start_points */



/************************************************************************/
/*                                                                      */
/*   Routine to update the sieve starting points when the sieve         */
/*   changes direction. i.e. lower half plane                           */
/*                                                                      */
/************************************************************************/
 
void reset_startpts()

{   /* start of reset_startpts */
register int i;

for (i=0; i<NUMSMALLP; i++)
   {
   if (int_startpts[i] != 0)
     int_startpts[i] = special_q_int_roots[i];  
   }
for (; i<int_linesieve-7; i+=4)
   {
   if (int_startpts[i] != 0)
      int_startpts[i] = special_q_int_roots[i];

   if (int_startpts[i+1] != 0)
      int_startpts[i+1] = special_q_int_roots[i+1];

   if (int_startpts[i+2] != 0)
      int_startpts[i+2] = special_q_int_roots[i+2];

   if (int_startpts[i+3] != 0)
      int_startpts[i+3] = special_q_int_roots[i+3];

   }
for (; i<int_linesieve; i++)
   {
   if (int_startpts[i] != 0)
      int_startpts[i] = special_q_int_roots[i];
   }

i = 0;
while (alg_fbase[i] <= SMALLP)
   {
   if (alg_startpts[i] != 0)
     alg_startpts[i] = special_q_alg_roots[i];
   i++;
   }
for (; i<alg_linesieve-7; i+=4)
   {
   if (alg_startpts[i] != 0)
      alg_startpts[i] = special_q_alg_roots[i];

   if (alg_startpts[i+1] != 0)
      alg_startpts[i+1] = special_q_alg_roots[i+1];

   if (alg_startpts[i+2] != 0)
      alg_startpts[i+2] = special_q_alg_roots[i+2];

   if (alg_startpts[i+3] != 0)
      alg_startpts[i+3] = special_q_alg_roots[i+3];

   }
for (; i<alg_linesieve; i++)
   {
   if (alg_startpts[i] != 0)
      alg_startpts[i] = special_q_alg_roots[i];
   }

if (DEBUG_FIND)
   {
   (void) printf("int reset\n");
   for (i=0; i<20; i++) (void) printf("%d ",int_startpts[i]);
   (void) printf("\n");
   (void) printf("alg reset\n");
   for (i=0; i<20; i++) (void) printf("%d ",alg_startpts[i]);
   (void) printf("\n");
   }

}   /* end of reset_startpts */
 
 
/************************************************************************/
/*                                                                      */
/*   routine to get scaled logarithms of the primes                     */
/*                                                                      */
/************************************************************************/
 
void get_logs(fbase,size,logs,which)
int size;
int fbase[];
unsigned char logs[];
int which;
 
{   /* start of get_logs */
 
register int i;
 
if (which == 1)  
   {
   for (i=0; i<size; i++)
     {
     if (i < NUMSMALLP)
        logs[i] = (unsigned char)(log((double)smallp[i])/LOGB + .5);
     else
        logs[i] = (unsigned char)(log((double)fbase[i])/LOGB + .5);
     }
   }
else
   {
   for (i=0; i<size; i++)
     { 
     logs[i] = (unsigned char)(log((double)fbase[i])/LOGB + .5);
     }
   }

}   /* end of get_logs */
 

/************************************************************************/
/*                                                                      */
/*   routine to get scaled log of test value for factorization          */
/*   computes algebraic side tolerance.                                 */
/*                                                                      */
/************************************************************************/
 
unsigned char get_alg_log(b,a)
int b,a;
 
{   /* start of get_alg_log */
double temp1,temp2,temp3;
 
/*   we require that the sum of the logs exceeds                */
/*   log(Norm(a + b alpha)).                                    */

a += MAX_SIEVE/2;                  /* midpoint of the interval  */

temp1 = log((double)a);
temp2 = log((double)b);
 
if (field_flag == GENFIELD)
   {
   temp3 = fabs(pcoef0*exp(temp1*degree) + pcoefd*exp(temp2*degree));
   temp3 = log(temp3)/LOGB - log_pmax_rhs;
   return((unsigned char)temp3);
   }

/*   We have a pure field. The threshold value is               */
/*   log(a^d - b^d field_base) - tol*log(pmax)                  */
/*   (except for very small b). for very small b we use a fudge */
/*   factor.                                                    */

temp3 = fabs(exp(temp1*degree) + fconst*exp(temp2*degree));
temp3 = log(temp3)/LOGB - log_pmax_rhs;
 
 
return((unsigned char)temp3);
 
}   /* end of get_alg_log */
  
/************************************************************************/
/*                                                                      */
/*         probable prime testing routine                               */
/*         works on arbitrary sized numbers                             */
/*                                                                      */
/************************************************************************/
 
int prp(number)
int number[];

{   /* start of prp */
 
//register int j;
int temp[MPDIM],temp1[MPDIM],base[MPDIM];

//COPY1(number,temp,j);
small_mpcopy(number,temp);
SIZE(base) = SINGLE;
LAST(base) = 101;               /* 101 as base for Fermat test      */

temp[1]--;                      /* use N-1 as the exponent          */

mpower(base,temp,number,temp1);      /* base^(N-1) MOD N            */

/*   check the residue: if it isn't 1 the number is definitely      */
/*   composite.                                                     */

if (SIZE(temp1) > 2 || LAST(temp1) != 1) return(1);
else return(0);
 
}   /* end of prp */


/************************************************************************/
/*                                                                      */
/*      routine to read an MP number from a file                        */
/*      max MPCDIM digits                                               */
/*                                                                      */
/************************************************************************/
 
void get_file(number,ptr)
int number[];
FILE *ptr;
 
{   /* start of get_file */

char string[MPCDIM];

(void) fscanf(ptr,"%s",string);
mpctoi(string,number);

}   /* end of get_file */
 
/************************************************************************/
/*                                                                      */
/*   routine to write parameters into restart file                      */
/*                                                                      */
/************************************************************************/

void dump_params(val,restname,outname)
int val;
char *restname,*outname;

{   /* start of dump_params */

rest = fopen(restname,"w");
(void) fprintf(rest,"%d\n",val);
(void) fprintf(rest,"%d %d %d %d %d %d %d %d %d %d\n",numtot,numff,numfp,numpf,numppf,
            numfpp, numpp, numppq, numqpp, numpppp);
(void) fflush(rest);
(void) fclose(rest);
(void) fclose(outfile);          /* close & reopen output; windows hack  */
outfile = fopen(outname,"a");

}   /* end of dump_params */
 
 
/************************************************************************/
/*                                                                      */
/*   routine to read in the factor base                                 */
/*                                                                      */
/************************************************************************/
 
get_ascii_fbase()
 
{   /* start of get_ascii_fbase */
int i;
 
FILE *binpc;
 
binpc  =  fopen(FACTOR_BASE_ASCII,"r");
if (binpc == NULL) return(0);

(void) fscanf(binpc,"%d\n",&field_base);
(void) fscanf(binpc,"%d\n",&poly_constant);
(void) fscanf(binpc,"%d\n",&lhs_degree);
(void) fscanf(binpc,"%d\n",&fdegree);
(void) fscanf(binpc,"%d\n",&size_int_base);
(void) fscanf(binpc,"%d\n",&size_alg_base);
 
for (i=0; i<size_int_base; i++)
   (void) fscanf(binpc,"%d %d\n",&int_fbase[i],&int_roots[i]);
 
for (i=0; i<size_alg_base; i++)
   (void) fscanf(binpc,"%d %d\n",&alg_fbase[i],&alg_roots[i]);
 
(void) fclose(binpc); 
 
int_pmax = int_fbase[size_int_base - 1];   /* largest int factor   */
alg_pmax = alg_fbase[size_alg_base - 1];   /* largest alg factor   */
 
return(1);
}   /* end of get_ascii_fbase */
 
/************************************************************************/
/*                                                                      */
/*      routine to read in the factor base                              */
/*                                                                      */
/************************************************************************/
 
void get_fbase()
 
{   /* start of get_fbase */

fbase = fopen(FACTOR_BASE_BIN,"r+b");
 
(void) fread((char*)&field_base,sizeof(int),1,fbase);
(void) fread((char*)&poly_constant,sizeof(int),1,fbase);
(void) fread((char*)&lhs_degree,sizeof(int),1,fbase);
(void) fread((char*)&fdegree,sizeof(int),1,fbase);
(void) fread((char*)&size_int_base,sizeof(int),1,fbase);
(void) fread((char*)&size_alg_base,sizeof(int),1,fbase);
(void) fread((char*)int_fbase,sizeof(int),(int)size_int_base,fbase);
(void) fread((char*)int_roots,sizeof(int),(int)size_int_base,fbase);
(void) fread((char*)alg_fbase,sizeof(int),(int)size_alg_base,fbase);
(void) fread((char*)alg_roots,sizeof(int),(int)size_alg_base,fbase);
(void) fclose(fbase);
 
int_pmax = int_fbase[size_int_base - 1];   /* largest int factor   */
alg_pmax = alg_fbase[size_alg_base - 1];   /* largest alg factor   */
 
}   /* end of get_fbase */
 
/************************************************************************/
/*                                                                      */
/*   Routine to test if two integers are relatively prime               */
/*   Uses binary method.                                                */
/*                                                                      */
/************************************************************************/
 
relprime(a,b)
int a,b;

{   /* start of relprime */
register unsigned int temp;

if (!(a & 1))                    /* if a is even       */
   {
   if (!(b & 1)) return(0);      /* and b is even      */
   else 
     {                           /* one is odd; swap   */
     temp = a;
     a = b;
     b = temp;
     }
   }

while (b != 0)
   {
   while (!(b & 1)) b >>= 1;     /* remove powers of two*/
   if (a > b)
     {
     temp = b;
     b = a - b;
     a = temp;
     }
   else b -= a;
   }

return (a == 1);

}   /* end of relprime */
 

/************************************************************************/
/*                                                                      */
/*   Routine to access the clock. On Pentium prectime returns ticks     */
/*                                                                      */
/************************************************************************/
 
double get_time()
 
{   /* start of get_time */

return((double)prectime());
   
}   /* end of get_time */
 
 
/************************************************************************/
/*                                                                      */
/*   Routine to compute the norm of a + b alpha.                        */
/*                                                                      */
/*   inputs:                                                            */
/*      b; obvious,  a; obvious                                         */
/*   returns:                                                           */
/*      answer; a multi-precision variable holding abs(norm)            */
/*      function returns the sign of the norm                           */
/*  can probably do this faster, but this takes very little time        */
/*                                                                      */
/************************************************************************/
int alg_norm(b,a, answer)
int a,b,answer[];
 
{   /* start of alg_norm */
int terms[MAX_DEGREE+1][MPDIM],termsign[MAX_DEGREE+1],
     a_powers[MAX_DEGREE+1][MPDIM],
     b_powers[MAX_DEGREE+1][MPDIM],
     negterms[MPDIM],plusterms[MPDIM];
int i,asign,bsign,sign;

if (DEBUG_ALG_NORM) (void) printf("norm(%d,%d): \n",a,b);

sign = SIGN(a);
a = abs(a);
b = abs(b);

a_powers[0][0] = SINGLE;
a_powers[0][1] = 1;
 
b_powers[0][0] = SINGLE;
b_powers[0][1] = 1;
 
/*      compute a^i and b^i  for i = 1...fdegree         */

for (i=1; i<=fdegree; i++)
   {
   mul_single(a_powers[i-1],a,a_powers[i]);
   mul_single(b_powers[i-1],b,b_powers[i]);
   if (DEBUG_ALG_NORM)
      {
      (void) printf("i, a^i = %d\n",i); write_num(a_powers[i],numstr);
      (void) printf("i, b^i = %d\n",i ); write_num(b_powers[i],numstr);
      }
   }

/*      coef * a^i * b^(d-i)  for i=1....fdegree         */

for (i=0; i<=fdegree; i++)
   {
   mult(a_powers[fdegree-i],b_powers[i],terms[i]);
   mul_single(terms[i],abs(poly_coef[i]),terms[i]);
   if (DEBUG_ALG_NORM) { (void) printf("%d : ",i);  write_num(terms[i],numstr); }
   }

/*   determine the sign of each term of the polynomial   */

if (fdegree & 1) bsign = -1;
else bsign = 1;
for (i=0; i<=fdegree; i++)
   {
   if (sign < 0) asign = 1;
   else 
      {
      if ((fdegree-i) & 1) asign = -1;
      else asign = 1;
      }
   termsign[i] = SIGN(poly_coef[i])*bsign*asign;
   if (DEBUG_ALG_NORM) 
     (void) printf("asign, bsign,i,termsign[i]: %d %d %d %d\n",asign,bsign,i,termsign[i]);
   }


/* add up all the terms; add the positive and negative terms separately   */
/* then subtact the negative term sum from the positive term sum          */

negterms[0] = SINGLE;
negterms[1] = 0;
plusterms[0] = SINGLE;
plusterms[1] = 0;
for (i=0; i<=fdegree; i++)
   {
   if (termsign[i] == 1) add(plusterms,terms[i],plusterms);
   else add(negterms,terms[i],negterms);
   }
if (DEBUG_ALG_NORM)
   {
   (void) printf("norm/ plusterms = "); write_num(plusterms, numstr);
   (void) printf("norm/ -terms = "); write_num(negterms, numstr);
   }
 
/*      subtract the odd terms from the even ones or vice versa           */

if (mpcmp(plusterms,negterms) == 1)
   {
   subt(plusterms,negterms,answer); 
   return(1);
   }
else
   {
   subt(negterms,plusterms,answer);
   return(-1);
   }

}   /* end of alg_norm */
 
/************************************************************************/
/*                                                                      */
/*   Subroutine to refine the threshhold value for an algebraic         */
/*   integer. The purpose is (hopefully) that a refined value will      */
/*   help reject some false hits. This will cut down on trial           */
/*   division                                                           */
/*                                                                      */
/*   we compute (-d)^deg f(-c/d) in floating point and take the log     */
/*   to base 2.                                                         */
/************************************************************************/
 
unsigned char refine_alg_log(j, index, sign)
int j, index, sign;

{   /* start of refine_alg_log */
double eval, sum;
int c,d,i;

compute_coeffs(sign*index, j, &c, &d);

if (DEBUG_REF_ALG) (void) printf("alg ref: (%d,%d), (%d,%d)\n",j, sign*index, c,d); 

eval = -(double)c/(double)d;

if (DEBUG_REF_ALG) (void) printf("eval,degree = %g %g\n",eval,degree);

sum = pcoef0;         /* use Horner's method for f(-a/b)   */
for (i=1; i<=fdegree; i++)
   {
   sum = sum*eval + polyc_float[i];
   }

sum = fabs(sum);

log_b = log(fabs((double)d))/LOGB;
eval = degree*log_b + log(sum)/LOGB - log_pmax_rhs;
if (eval < 10.0) eval = 10.0;

if (DEBUG_REF_ALG) (void) printf("sum, log_b, eval = %g %g %g\n",sum,log_b, eval);

return((unsigned char) eval);

}   /* end of refine_alg_log */


/************************************************************************/
/*                                                                      */
/*   Subroutine to refine the threshhold value for a rational           */
/*   integer. The purpose is (hopefully) that a refined value will      */
/*   help reject some false hits. This will cut down on trial           */
/*   division                                                           */
/*                                                                      */
/*  Int norm =  a * int_poly[0] + b * int_poly[1]                       */
/*   normally, one term will be larger. we computer both and compare    */
/*                                                                      */
/************************************************************************/
 
unsigned char refine_int_log(loc, column)
int column, loc;

{   /* start of refine_int_log */
double part1, part2, eval;
int c,d;

compute_coeffs(loc, column, &c, &d);

part1 = log(fabs((double)d))/LOGB + log_PART1 - log_pmax_lhs;
part2 = log(fabs((double)c))/LOGB + log_PART2 - log_pmax_lhs;

eval = part2;
if (part1 > part2) eval = part1;

if (DEBUG_REF_INT)
   {
   (void) printf("refine int value: (%d,%d), (%d,%d) = %g\n",loc,column,c,d,eval);
   }

return((unsigned char) eval);

}   /* end of refine_int_log */

/************************************************************************/
/*                                                                      */
/*   Routine to factor a number via Shanks' SQUFOF                      */
/*                                                                      */
/*   Returns squfof = 1 and factors in fac1 <= fac2 if successful       */
/*   Returns squfof = 0 if unsuccessful                                 */
/*                                                                      */
/*   Changed in 11/04. First call QS. (it is faster). Then call SQUFOF  */
/*  if QS fails                                                         */
/*                                                                      */
/************************************************************************/

squfof(n,fac1,fac2)
int n[];
int *fac1, *fac2;

{   /* start of squfof: factor n as fac1*fac2  faster in FP?????*/

int temp[MPDIM], temp1[MPDIM];
register int iq,ll,l2,p,pnext,q,qlast,r,s,t,i;
int jter,iter;

if (SQUFOF) squfof_calls++;
#define   DEBUG_QSC   0
if (DEBUG_QSC) { (void) printf("trying qs on: "); write_num(n,numstr); fflush(stdout); }
jter = small_qs(n,fac1,fac2);
if (DEBUG_QSC)
   {
   (void) printf("qs returns %d\n",jter);
   (void) printf("tried qs on "); write_num(n, numstr);
   (void) printf("fac1, fac2 = %d %d, squfof returning %d\n", *fac1, *fac2, jter);
   fflush(stdout);
   }         /* -1 should be very rare; QS found no non-trivial GCD's   */
if (jter != -1) return(jter); 

/*   Notes: could replace calls to the MP routines [on Pentium] by using*/
/*  __int64  for temp,temp1, convert n to __int64 etc. But the MP stuff */
/*  is not inside a loop. It is not worth doing                         */

qlast = 1;
mpsqrt(n,temp);

s = LAST(temp);
p = s;
mult(temp,temp,temp1);
subt(n,temp1,temp);                    /* temp = n - floor(sqrt(n))^2   */
if (SIZE(temp) <= SINGLE && LAST(temp) == 0)
   {                                   /* Here n is a square            */
   *fac1 = s;
   *fac2 = s;
   return(1);
   }

q = LAST(temp);              /* q = excess of n over next smaller square */
ll = 1 + 2*(int)sqrt((double)(p+p));
l2 = ll/2;
qpoint = 0;

/*   In the loop below, we need to check if q is a square right before   */
/*  the end of the loop.  Is there a faster way? The current way is      */
/*   EXPENSIVE! (many branches and double prec sqrt)                     */

for (jter=0; jter < 800000; jter++)      /* I see no way to speed this   */
   {                                     /*  main loop                   */
   iq = (s + p)/q;   
   pnext = iq*q - p;
   if (q <= ll)
      {
      if ((q & 1) == 0) enqu(q/2,&jter);
      else if (q <= l2) enqu(q,&jter);
      if (jter < 0)
          {                        /* queue overflow: try pollard rho   */
         if (Q_OVERFLOW) overflow_count++;    /* should be rare         */
          return(pollardrho2(n,fac1,fac2));   /* also try squfof(3N)?   */
          }
      }
   t = qlast + iq*(p - pnext);
   qlast = q;
   q = t;
   p = pnext;                          /* check for square; even iter   */
   if (jter & 1) continue;             /* jter is odd:omit square test  */
//   if (t = (q & 3) > 1) continue;    /* skip t = 2,3 mod 4            */
//   if (q & 7 == 5) continue;         /* skip q = 5 mod 8              */
//   if (t == 0 && (q & 15) > 5) continue;   /* skip 8, 12 mod 16       */
//   if (s3465[q % 3465] == 'n') continue;   /* skip QNR's mod 3465     */
   r = (int)sqrt((double)q);                 /* r = floor(sqrt(q))      */
   if (q != r*r) continue;
   if (qpoint == 0) goto gotit;
   for (i=0; i<qpoint-1; i+=2)      /* treat queue as list for simplicity*/
      {
      if (r == qqueue[i]) goto contin;
      if (r == qqueue[i+1]) goto contin;
      }
   if (r == qqueue[qpoint-1]) continue;
   goto gotit;
contin:;   
   }   /* end of main loop */

if (DEBUG) { (void) printf("calling rho on: "); write_num(n,numstr); }
if (RHO_TALLY) rho_calls++;

return(pollardrho2(n,fac1,fac2));         /* very rare                 */

gotit:   ;
qlast = r;
p = p + r*((s - p)/r);
SIZE(temp1) = SINGLE;
LAST(temp1) = p;
mult(temp1,temp1,temp);
subt(n,temp,temp);
div_single(temp,qlast,temp1);
q = LAST(temp1);                  /* q = (n - p*p)/qlast (div is exact)*/
for (iter=0; iter<40000; iter++)
   {                              /* begin second main loop            */
   iq = (s + p)/q;                /* unroll it, of course              */
   pnext = iq*q - p;
   if (p == pnext) goto gotfac;
   t = qlast + iq*(p - pnext);
   qlast = q;
   q = t;
   p = pnext;
   iq = (s + p)/q;
   pnext = iq*q - p;
   if (p == pnext) goto gotfac;
   t = qlast + iq*(p - pnext);
   qlast = q;
   q = t;
   p = pnext;
   iq = (s + p)/q;
   pnext = iq*q - p;
   if (p == pnext) goto gotfac;
   t = qlast + iq*(p - pnext);
   qlast = q;
   q = t;
   p = pnext;
   iq = (s + p)/q;
   pnext = iq*q - p;
   if (p == pnext) goto gotfac;
   t = qlast + iq*(p - pnext);
   qlast = q;
   q = t;
   p = pnext;
   }

if (SQUFOF_FAIL) squfof_fail++;
return(0);                               /* this shouldn't happen      */

gotfac:   ; if ((q & 1) == 0) q/=2;      /* q was factor or 2*factor   */
*fac1 = q;
t = div_single(n,q,temp);
if (SIZE(temp) > SINGLE || t != 0) 
   {
   if (SPLITCHECK) badsplit++;
   return (0);
   }
*fac2 = LAST(temp);
return(1);
}

/************************************************************************/
/*                                                                      */
/*      Routines to prepare tables for squfof algorithm                 */
/*                                                                      */
/************************************************************************/
 

void sqfprep()
{
int i;
for (i=0; i<3465; i++) s3465[i] = 'n';
for (i=0; i<1733; i++) s3465[(i*i) % 3465] = 'r';
}

void enqu(q,iter)
int q;
int *iter;
{
qqueue[qpoint] = q;
if (++qpoint > 50) *iter = -1;
}


/************************************************************************/
/*                                                                      */
/*   Routine to output relations to the nfs.out file                    */
/*                                                                      */
/************************************************************************/

void output_relation(ilp1, ilp2, alp, alp2, sign, a, b)
int ilp1,ilp2;
int alp,alp2,sign,a,b;

{   /* start of output_relation */
int whichcase,tmp;
int icount,acount, ifirst, afirst, i;
int nint, nalg;

if (ilp1 < ilp2)
   {
   tmp = ilp2;
   ilp2 = ilp1;
   ilp1 = tmp;
   }

if (alp < alp2)
   {
   tmp = alp2;
   alp2 = alp;
   alp = tmp;
   }

if (DEBUG_OUTPUT) 
   (void) printf("ilp1,ilp2,alp = %d %d %d %d\n",ilp1,ilp2,alp,alp2);

whichcase = 0;
if (ilp1 > 1) whichcase += IL1;
if (ilp2 > 1) whichcase += IL2;
if (alp  > 1) whichcase += AL1;
if (alp2 > 1) whichcase += AL2;

if (DEBUG_OUTPUT) (void) printf("whichcase = %d\n",whichcase);
if (DEBUG_OUTPUT) (void) printf("a,b = %d %d\n",a,b);

numtot++;               /* update total number of factorizations   */

icount = int_factor_list[0];
acount = alg_factor_list[0];

if (DEBUG_OUTPUT)
   {
   (void) printf("icount,acount = %d %d\n",icount,acount);
   for (i=1; i<=icount; i++) (void) printf("%d ",int_factor_list[i]); printf("\n");
   for (i=1; i<=acount; i++) (void) printf("%d ",alg_factor_list[i]); printf("\n");
   (void) fflush(stdout);
   }

sort_int_factor_list(int_factor_list);

if (DEBUG_OUTPUT)
   {
   (void) printf("sorted list: ");
   for (i=0; i<=int_factor_list[0]; i++) (void) printf("%d ",int_factor_list[i]);
   (void) printf("\n");
   }
 
ifirst = icount+1;
for (i=1; i<=icount; i++)
   {
   if (int_factor_list[i] > dump_cutoff)
   {
   ifirst = i;
   break;
   }
   }
 
afirst = acount+1;
for (i=1; i<=acount; i++)
   {
   if (alg_factor_list[i] > dump_cutoff)
   {
   afirst = i;
   break;
   }
   }
 
if (DEBUG_OUTPUT) { (void) printf("ifirst,afirst = %d %d\n",ifirst,afirst); (void) fflush(stdout); }

nint = icount - ifirst + 1;
nalg = acount - afirst + 1;

if (DEBUG_OUTPUT) { (void) printf("#int, #alg = %d %d\n",nint,nalg); (void) fflush(stdout); }

switch(whichcase)
   {
   case FF:
     numff++;
      (void) fprintf(outfile,"01%c%c %d %d",nint+'0',nalg+'0',-a,b);
     dump_int_factors(ifirst,icount);
     dump_alg_factors(afirst,acount);
     (void) fprintf(outfile,";\n");
     break;
   case FP:
     numfp++;
      (void) fprintf(outfile,"01%c%c %d %d",nint+'0',nalg+'0'+1,-a,b);
     dump_int_factors(ifirst,icount);
      (void) fprintf(outfile," %d",alp);
     dump_alg_factors(afirst,acount);
     (void) fprintf(outfile,";\n");
     break;
   case PF:
     numpf++;
      (void) fprintf(outfile,"01%c%c %d %d %d",nint+'0'+1,nalg+'0',-a,b,ilp1);
     dump_int_factors(ifirst,icount);
     dump_alg_factors(afirst,acount);
     (void) fprintf(outfile,";\n");
     break;
   case PPF:
     numppf++;
      (void) fprintf(outfile,"01%c%c %d %d %d %d",nint+'0'+2,nalg+'0',-a,b,ilp1,ilp2);
     dump_int_factors(ifirst,icount);
     dump_alg_factors(afirst,acount);
     (void) fprintf(outfile,";\n");
     break;
   case PP:
     numpp++;
      (void) fprintf(outfile,"01%c%c %d %d %d",nint+'0'+1,nalg+'0'+1,-a,b,ilp1);
     dump_int_factors(ifirst,icount);
      (void) fprintf(outfile," %d",alp);
     dump_alg_factors(afirst,acount);
     (void) fprintf(outfile,";\n");
     break;
   case PPQ:
     numppq++;
      (void) fprintf(outfile,"01%c%c %d %d %d %d",nint+'0'+2,nalg+'0'+1,-a,b,ilp1,ilp2);
     dump_int_factors(ifirst,icount);
      (void) fprintf(outfile," %d",alp);
     dump_alg_factors(afirst,acount);
     (void) fprintf(outfile,";\n");
     break;
   case QPP:
     numqpp++;
      (void) fprintf(outfile,"01%c%c %d %d %d",nint+'0'+1,nalg+'0'+2,-a,b,ilp1);
     dump_int_factors(ifirst,icount);
      (void) fprintf(outfile," %d %d",alp,alp2);
     dump_alg_factors(afirst,acount);
     (void) fprintf(outfile,";\n");
     break;
   case FPP:
     numfpp++;
      (void) fprintf(outfile,"01%c%c %d %d",nint+'0',nalg+'0'+2,-a,b);
     dump_int_factors(ifirst,icount);
      (void) fprintf(outfile," %d %d",alp,alp2);
     dump_alg_factors(afirst,acount);
     (void) fprintf(outfile,";\n");
     break;
   case PPPP:
     numpppp++;
      (void) fprintf(outfile,"01%c%c %d %d %d %d",nint+'0'+2,nalg+'0'+2,-a,b,ilp1,ilp2);
     dump_int_factors(ifirst,icount);
      (void) fprintf(outfile," %d %d",alp,alp2);
     dump_alg_factors(afirst,acount);
     (void) fprintf(outfile,";\n");
     break;
   }

(void) fflush(outfile);
 
}   /* end of output_relation */
 
/************************************************************************/
/*                                                                      */
/*   Routine to dump integer side large primes to output file           */
/*                                                                      */
/************************************************************************/
 
void dump_int_factors(start,end)
int start,end;

{   /* start of dump_int_factors */
int i;

if (DEBUG_OUTPUT) { (void) printf("dumping int: start,end = %d %d\n",start,end); }

for (i=end; i>=start; i--) (void) fprintf(outfile," %d",int_factor_list[i]);

}   /* end of dump_int_factors */
 
/************************************************************************/
/*                                                                      */
/*   Routine to dump algebraic side large primes to output file         */
/*                                                                      */
/************************************************************************/
 
void dump_alg_factors(start,end)
int start,end;

{   /* start of dump_alg_factors */
int i;

if (DEBUG_OUTPUT) 
   {
   (void) printf("dumping alg: start,end = %d %d\n",start,end); 
   (void) fflush(stdout); 
   }
 
for (i=end; i>=start; i--) (void) fprintf(outfile," %d",alg_factor_list[i]);

}   /* end of dump_alg_factors */

 
/************************************************************************/
/*                                                                      */
/*   Debug routine; now obsolete/not needed                             */
/*                                                                      */
/************************************************************************/
void dump_list(a)
int a[];

{
int i;

(void) printf("list %d:\n",a[0]);
for (i=1; i<=a[0]; i++) (void) printf("%d ",a[i]);
(void) printf("\n");
}


/************************************************************************/
/*                                                                      */
/*   Routine to compute poly roots for powers of algebraic              */
/*   small primes. We do this by direct search since the primes are tiny*/
/*                                                                      */
/************************************************************************/

void prime_power_roots()

{   /* start of prime_power_roots */
int i,j,smallpow;
int prime,rem;
int temp[MPDIM],num[MPDIM];

/*   We try to find a root of  f(x) mod alg_smallp[i]                   */
/*  If a root mod a power exists, use the power instead of the prime    */
/*   itself. If the root exists reset the log value to log(p^d)         */

power_root_count = 0;
i = 0;
while(alg_fbase[i] <= SMALLP)
   {
   prime = alg_fbase[i];
   smallpow = alg_smallp[prime];
   j = prime - alg_roots[i];
    alg_root_ok[i] = FALSE;       /* assume no root until we find one   */
   while (j < smallpow)
      {
      alg_norm(1,j,num);
      rem = div_single(num, smallpow, temp);
      if (rem == 0)
         {
         if (DEBUG_FIND) (void) printf("Alg prime power root: %d %d\n",alg_fbase[i],j);
         alg_prime_power_root[i] = j;
         alg_root_ok[i] = TRUE;
         power_root_count++;
         alg_logp[i] = (unsigned char)(log((double)alg_smallp[alg_fbase[i]])/LOGB + .5);
         break;
         }
      j += prime;
      }
   i++;
   }         

if (DEBUG) (void) printf("power root count = %d\n",power_root_count);

}   /* end of prime_power_roots */


/************************************************************************/
/*                                                                      */
/*   routine to initialize list to hold resieve factors                 */
/*                                                                      */
/************************************************************************/

void init_faclist()

{   /* start of init_faclist */ 
int i;

for (i=0; i<TSIZE; i++)      
   {
   int_faclist[i][0] = 0; 
   alg_faclist[i][0] = 0; 
   }
}   /* end of init_faclist */


 
/************************************************************************/
/*                                                                      */
/*   routine to factor successes by re-sieving                          */
/*                                                                      */
/************************************************************************/

void factor_by_resieve(sign) 
int sign;                

{   /* start of factor_by_resieve */
int i,c,d,c1,d1,c0,d0;
int e,e_min,e_max;                     /*  Pollard's notation  */
int a0, b0, a1, b1;
int f_lower;
double itime,atime;


/*   we resieve only with respect to primes >= MAXSIEVE ~ 8K      */
/*   we therefore resieve by vectors                              */

if (DEBUG_RESIEVE) { (void) printf("Vector INT re-sieve\n"); fflush(stdout); }

init_faclist();
 
if (sign == 1)   {         /* upper half plane; e > 0               */

if (TIME_STATS) itime = get_time();

do_delayed_int_resieve_points ();

for (i=int_linesieve;  i<size_int_base; i++)
   {
   if (DEBUG_INT_VECTOR_R) { (void) printf("resieve vector[%d] = ", int_fbase[i]); show_int_vector(i); }
   if (int_v1[i][0] == 0) continue;   /* no hits      */

   prepare_bounds(int_v1[i], int_v2[i], sign);
   lc1 += .50;
   lc2 += .50;

   if (SHOW_INT_REGION_R) show_region(int_fbase[i]);

   a0 = int_v1[i][0];
   b0 = int_v1[i][1];
   a1 = int_v2[i][0];
   b1 = int_v2[i][1];
   e_min = (int)emin;
   e_max = (int)emax;

   if (DEBUG_INT_VECTOR_R) (void) printf("INT resieve emin, emax = %d %d\n",e_min,e_max);
   if (DEBUG_INT_VECTOR_R) (void) printf("i, a0,b0, a1, b1 = %d %d %d %d %d\n",i,a0,b0,a1,b1);
                           
   c0 = e_min * a0;
   d0 = e_min * b0;
   for (e=e_min; e<=e_max; e++)
      {       
//    f_lower = (int)lower_bound((double)e);

	  if (e < lower_midpt) f_lower = round_double(lm1 * (double)e + lc1);
	  else                 f_lower = round_double(lm2 * (double)e + lc2);

      c1 = f_lower * a1;
      d1 = f_lower * b1; 
      c = c0 + c1;
      d = d0 + d1;
      if (c == 0 && d == 0) c = a1, d = b1;

// Try to add this entire vector to the delayed update list

      if (! add_vector_entry (d, c, b1, a1, i, NULL, 0))
         goto int_done;
      c0 += a0; 
      d0 += b0; 
      } 
   }
   do_delayed_int_resieve_vectors ();
int_done:

if (TIME_STATS) int_resieve_time += (get_time() - itime);

if (DEBUG_RESIEVE) { (void) printf("Vector ALG re-sieve\n"); fflush(stdout); }

if (TIME_STATS) atime = get_time();

do_twice_delayed_alg_resieve_points ();

for (i=alg_linesieve;  i<size_alg_base; i++)
   {
   if (DEBUG_ALG_VECTOR_R) { (void) printf("vector[%d] = ", alg_fbase[i]); show_alg_vector(i); }
   if (alg_v1[i][0] == 0) continue;   /* no hits   */

   prepare_bounds(alg_v1[i], alg_v2[i], sign);
   lc1 += .50;
   lc2 += .50;

   if (SHOW_ALG_REGION) show_region(alg_fbase[i]);

   a0 = alg_v1[i][0];
   b0 = alg_v1[i][1];
   a1 = alg_v2[i][0];
   b1 = alg_v2[i][1];
   e_min = (int)emin;
   e_max = (int)emax;

   if (DEBUG_ALG_VECTOR_R) (void) printf("ALG emin, emax = %d %d\n",e_min,e_max);
   if (DEBUG_ALG_VECTOR_R) (void) printf("i, a0,b0, a1, b1 = %d %d %d %d %d\n",i,a0,b0,a1,b1);

   c0 = e_min * a0;
   d0 = e_min * b0;
   for (e=e_min; e<=e_max; e++)
      {       
//    f_lower = (int)lower_bound((double)e);

	  if (e < lower_midpt) f_lower = round_double(lm1 * (double)e + lc1);
	  else                 f_lower = round_double(lm2 * (double)e + lc2);

      c1 = f_lower * a1;
      d1 = f_lower * b1; 
      c = c0 + c1;
      d = d0 + d1;
      if (c == 0 && d == 0) c = a1, d = b1;

// Try to add this entire vector to the delayed update list

      if (! add_vector_entry (d, c, b1, a1, i, NULL, 1))
         goto alg_done;
      c0 += a0; 
      d0 += b0;
      }                            }
   do_delayed_alg_resieve_vectors ();
alg_done:

if (TIME_STATS) alg_resieve_time += (get_time() - atime);

}   /* end of sign = 1 */

else {                              /* lower half plane; sign = -1   */

if (TIME_STATS) itime = get_time();

do_delayed_int_resieve_points ();

for (i=int_linesieve;  i<size_int_base; i++)
   {
   if (DEBUG_INT_VECTOR_R) { (void) printf("resieve vector[%d] = ", int_fbase[i]); show_int_vector(i); }
   if (int_v1[i][0] == 0) continue;

   prepare_bounds(int_v1[i], int_v2[i], sign);
   lc1 += .50;
   lc2 += .50;

   if (SHOW_INT_REGION_R) show_region(int_fbase[i]);

   a0 = int_v1[i][0];
   b0 = int_v1[i][1];
   a1 = int_v2[i][0];
   b1 = int_v2[i][1];
   e_min = (int)emin;
   e_max = (int)emax;

   if (DEBUG_INT_VECTOR_R) (void) printf("INT resieve emin, emax = %d %d\n",e_min,e_max);
   if (DEBUG_INT_VECTOR_R) (void) printf("i, a0,b0, a1, b1 = %d %d %d %d %d\n",i,a0,b0,a1,b1);                        
   c0 = e_min * a0;
   d0 = e_min * b0;
   for (e=e_min; e<=e_max; e++)
      {       
//    f_lower = (int)lower_bound(e);

	  if (e < lower_midpt) f_lower = round_double(lm1 * (double)e + lc1);
	  else                 f_lower = round_double(lm2 * (double)e + lc2);

      c1 = f_lower * a1;
      d1 = f_lower * b1; 
      c = c0 + c1;
      d = d0 + d1;
      if (c == 0 && d == 0) c = a1, d = b1;

// Try to add this entire vector to the delayed update list

      if (! add_vector_entry (d, -c, b1, -a1, i, NULL, 0))
         goto int_done2;
      c0 += a0; 
      d0 += b0;
      } 
   }
   do_delayed_int_resieve_vectors ();
int_done2:

if (TIME_STATS) int_resieve_time += (get_time() - itime);

if (TIME_STATS) atime = get_time();

do_twice_delayed_alg_resieve_points ();

for (i=alg_linesieve;  i<size_alg_base; i++)
   {
   if (DEBUG_ALG_VECTOR_R) { (void) printf("vector[%d] = ", alg_fbase[i]); show_alg_vector(i); }
   if (alg_v1[i][0] == 0) continue;

   prepare_bounds(alg_v1[i], alg_v2[i], sign);
   lc1 += .50;
   lc2 += .50;

   if (SHOW_ALG_REGION) show_region(alg_fbase[i]);

   a0 = alg_v1[i][0];
   b0 = alg_v1[i][1];
   a1 = alg_v2[i][0];
   b1 = alg_v2[i][1];
   e_min = (int)emin;
   e_max = (int)emax;

   if (DEBUG_ALG_VECTOR_R) (void) printf("ALG emin, emax = %d %d\n",e_min,e_max);
   if (DEBUG_INT_VECTOR_R) (void) printf("i, a0,b0, a1, b1 = %d %d %d %d %d\n",i,a0,b0,a1,b1);
                           
   c0 = e_min * a0;
   d0 = e_min * b0;
   for (e=e_min; e<=e_max; e++)
      {       
//    f_lower = (int)lower_bound(e);

	  if (e < lower_midpt) f_lower = round_double(lm1 * (double)e + lc1);
	  else                 f_lower = round_double(lm2 * (double)e + lc2);

      c1 = f_lower * a1;
      d1 = f_lower * b1; 
      c = c0 + c1;
      d = d0 + d1;
      if (c == 0 && d == 0) c = a1, d = b1;

// Try to add this entire vector to the delayed update list

      if (! add_vector_entry (d, -c, b1, -a1, i, NULL, 1))
         goto alg_done2;
      c0 += a0;
      d0 += b0;
      }          
   }
   do_delayed_alg_resieve_vectors ();
alg_done2:

if (TIME_STATS) alg_resieve_time += (get_time() - atime);

}   /* end of lower half plane */

if (DEBUG_RESIEVE) { (void) printf("Finished ReSieve\n"); fflush(stdout); }

}   /* end of factor_by_resieve */

/************************************************************************/
/*                                                                      */
/*   routine to read the RANGE file and get a range                     */
/*   records are of the form start_q end_q  flag                        */
/*   where flag = 0,1,2 for 'available', 'in use', 'done'               */
/*                                                                      */
/************************************************************************/

void get_range(start_q, end_q, line_num)
int *start_q, *end_q,  *line_num;

{   /* start of get_range */
int fail_count,i,start,end,flag;
FILE *rangefile;

fail_count = 0;

while (1)
   {                            /* if accessing over a LAN, may not be able */
   if (fail_count == 100)       /* to open if another job has this open     */
      {                         /* so try again.                            */
      (void) printf("Open failure on RANGE file\n");
      exit(0);
      }
   rangefile = fopen("lrange","r+");
   if (rangefile == NULL)      /* can't open; sleep and try again         */
      {
      Sleep(SLEEPTIME);
      fail_count++;
      continue;
      }
   else break;
   }

for (i=0; i<MAX_LINES; i++)
   {
   (void) fscanf(rangefile,"%d %d %d",&start,&end,&flag);
   if (flag == 0)                  /* check for availability         */
      {
      *start_q = start;
      *end_q = end;
      *line_num = i+1;
      (void) fseek(rangefile, -(long)sizeof(char), SEEK_CUR);   /* mark as 'in use' */
      (void) fprintf(rangefile,"1");
      (void) fflush(rangefile);
      (void) fclose(rangefile);
      return;
      }
   }

(void) printf("Out of ranges\n");
exit(0);

}   /* end of get_range */


/************************************************************************/
/*                                                                      */
/*   routine to mark the RANGE file that a range has been finished      */
/*                                                                      */
/************************************************************************/

void mark_range_as_used(q_value,line_num)
int q_value,line_num;

{   /* start of mark_range_as_used */
int fail_count,i,start,end,flag;
FILE *rangefile;

fail_count = 0;

while (1)
   {
   if (fail_count == 200)
      {
      (void) printf("Open failure on marking RANGE file\n");
      exit(0);
      }
   rangefile = fopen("lrange","r+");
   if (rangefile == NULL)         /* sleep for .3 sec and try again   */
      {
      Sleep(SLEEPTIME);
      fail_count++;
      continue;
      }
   else break;
   }

for (i=0; i<MAX_LINES; i++)            /* mark this range as done      */
   {
   (void) fscanf(rangefile,"%d %d %d",&start,&end,&flag);
   if (start == q_value) 
      {
      (void) fseek(rangefile, -(long)sizeof(char), SEEK_CUR);
      (void) fprintf(rangefile,"2");
      (void) fflush(rangefile);
      (void) fclose(rangefile);
      return;
      }
   }

(void) printf("Could not mark %d as done.\n",q_value);

}   /* end of mark_range_as_used */


/************************************************************************/
/*                                                                      */
/*   routine to update the global count file (of successes)             */
/*                                                                      */
/************************************************************************/

void update_global_count()

{   /* start of update_global_count */
FILE *global;
int  gnumtot,gnumff,gnumfp,gnumpf,gnumppf, gnumfpp, gnumpp, gnumppq, gnumqpp, gnumpppp;
int fail_count;

fail_count = 0;

while (1)
   {
   if (fail_count == 200)
     {
     (void) printf("Open failure on global count file\n");
     exit(0);
     }
   global = fopen("global","r");
   if (global == NULL)      /* can't open; sleep for 3 sec and try again */
     {
     Sleep(SLEEPTIME);
     fail_count++;
     continue;
     }
   else break;
   }

(void) fscanf(global,"%d %d %d %d %d %d %d %d %d %d\n",&gnumtot,&gnumff,&gnumfp,&gnumpf,&gnumppf,
            &gnumfpp, &gnumpp, &gnumppq, &gnumqpp, &gnumpppp);

(void) fclose(global);

gnumtot += numtot;
gnumff  += numff;
gnumfp  += numfp;
gnumpf  += numpf;
gnumppf += numppf;
gnumfpp += numfpp; 
gnumpp  += numpp;
gnumppq += numppq;
gnumqpp += numqpp;
gnumpppp+= numpppp;

fail_count = 0;

while (1)
   {
   if (fail_count == 200)
      {
      (void) printf("Open failure on global count file\n");
      exit(0);
      }
   global = fopen("global","w");
   if (global == NULL)         /* can't open; sleep for 3 sec and try again */
      {
      Sleep(SLEEPTIME);
      fail_count++;
      continue;
      }
   else break;
   }

(void) fprintf(global,"%d %d %d %d %d %d %d %d %d %d\n",gnumtot,gnumff,gnumfp,gnumpf,gnumppf,
            gnumfpp, gnumpp, gnumppq, gnumqpp, gnumpppp);

(void) fclose(global);


}  /* end of update_global_count */


/************************************************************************/
/*                                                                      */
/*   routine to compute (c,d) = i V1 + j V2                             */
/*                                                                      */
/************************************************************************/

void compute_coeffs(i, j, c, d)
int i,j, *c, *d;

{   /* start of compute_coeffs */

*c = i*v1[0] + j*v2[0];
*d = i*v1[1] + j*v2[1];


}   /* end of compute_coeffs */

#define   FEWERDIVLATRED   1
#define   NODIVLATRED      0

/************************************************************************/
/*                                                                      */
/*   routine to compute two reduced lattice vectors                     */
/*  borrowed from AKL; which explains stylistic differences in code     */
/*  this routine is old;  I re-wrote it;  the new one is much faster    */
/*                                                                      */
/************************************************************************/

reduce_lattice(special_q, root, v1, v2)
int special_q, root;
int *v1, *v2;

{   /* start; of reduce_lattice */
double stime;
long x1; 
long yy1 = 0;
long x2;
long y2 = 1;
long x3;
long y3; 
double s2;
double s3;
long q; 
long sx2;
long parity = 0;
//int v1a[2], v2a[2];

#if (FEWERDIVLATRED | NODIVLATRED)
double s2_2;
double s2_4;
double s2_8;
double s2_16;
#endif

if (TIME_STATS) stime = get_time();

x1 = special_q;
x2 = root;

s2 = x2 * (double) x2 + y2 * (double) y2;

for (;;)
   {
   x3 = x1;
   y3 = yy1;
   sx2 = (x2 << 4);
   while (x3 > sx2 )
      {
      x3 -= sx2;
      y3 -= (y2 << 4);
      }
   sx2 >>= 1;
   if (x3 > sx2 )
      {
      x3 -= sx2;
      y3 -= (y2 << 3);
      }
   sx2 >>= 1;
   if (x3 > sx2 )
      {
      x3 -= sx2;
      y3 -= (y2 << 2);
      }
   sx2 >>= 1;
   if (x3 > sx2 )
      {
      x3 -= sx2;
      y3 -= (y2 << 1);
      }
   if (x3 > x2 )
      {
      x3 -= x2;
      y3 -= y2;
      }
   if (x3 > (x2 >> 1))
      {
      x3 -= x2;
      y3 -= y2;
      }
   if (x3 < 0)
      {
      x3 = -x3;
      y3 = -y3;
      parity ^= 1;
      }

   s3 = x3 * (double) x3 + y3 * (double) y3;
   if (s3 >= s2) break;
   s2 = s3;
   x1 = x2;
   yy1 = y2;
   x2 = x3;
   y2 = y3;
   parity ^= 1;
   }

#if (FEWERDIVLATRED | NODIVLATRED)
s3 = x2 * (double) x3 + y2 * (double) y3;
if (s3 < 0)
   {
   yy1 = 1;
   s3 = -s3;
   }
else
   {
   yy1 = 0;
   }
s2_2 = s2 + s2;
s2_4 = s2_2 + s2_2;
s2_8 = s2_4 + s2_4;
s2_16 = s2_8 + s2_8;
q = 0;
#ifdef NODIVLATRED
while (s3 > s2_16)
#else
again:
if (s3 > s2_16)
#endif
   {
   s3 -= s2_16;
   q += 16;
   }
if (s3 > s2_8)
   {
   s3 -= s2_8;
   q += 8;
   }
if (s3 > s2_4)
   {
   s3 -= s2_4;
   q += 4;
   }
if (s3 > s2_2)
   {
   s3 -= s2_2;
   q += 2;
   }
if (s3 > s2)
   {
   s3 -= s2;
   q ++;
   }
#ifndef NODIVLATRED
if (s3 > s2)
   {
   if (s3 < (s2_16+s2_16)) goto again;
   q += (long)mrint(s3/s2);
   }
else
#endif
if ((s3+s3) >= s2) q ++;
if (yy1) q = -q;
#else
q = (long) mrint( (x2 * (double) x3 + y2 * (double) y3) / s2 );
#endif

v2[0] = x2;
v2[1] = y2;
v1[0] = x3 - x2 * q;
v1[1] = y3 - y2 * q;

if (v2[1] < 0)      /* make 'b'  positive */
   {
   v2[0] = -v2[0];
   v2[1] = -v2[1];
   }

if (v1[1] < 0)      /* make 'b'  positive */
   {
   v1[0] = -v1[0];
   v1[1] = -v1[1];
   }

if (TIME_STATS) latred_time += (get_time() - stime);
//qtime = (get_time() - stime);

//time = get_time();
//newest_latred_exp(special_q, root, v1a, v2a);
//ttime = (get_time() - stime);

//printf("old, new = %lf %lf\n", qtime, ttime);
//printf("RATIO: %lf\n", qtime/ttime);
//tot1 += qtime;
//tot2 += ttime;

//cmp_latred(special_q, root, v1,v2, v1a, v2a);

return(parity);
}   /* end of reduce_lattice */


/************************************************************************/
/*                                                                      */
/*   Find success locations. Search among the ALGEBRAIC successes only. */
/*  When the result of the rational sieve gives a success for that      */
/*   location, mark it as GOOD.   Note that algebraic successes were    */
/*  stored in order of increasing column, so the test value can be      */
/*   updated at regular intervals.                                      */
/*                                                                      */
/************************************************************************/

void find_success_locations()

{   /* start of find_success_locations */
int loc, column;
int i, refine_count;
unsigned char testval;

if (DEBUG_INT_SUCCESS) printf("Finding success locations\n");

refine_count = 0;
num_total_success = 0;
for (i=0; i<num_alg_success; i++)
   {
   loc = success[i][LOC_STORE];
   column = success[i][COLUMN_STORE];
   testval = refine_int_log(loc,column);   /* is this right? do every time? */
    if (DEBUG_SUCCESS) { (void) printf("scanning: %d %d %d %d\n",i, column,loc,testval); fflush(stdout); }

   refine_count++;
   if (refine_count == INT_REFINE_FREQ) 
      {
      testval = refine_int_log(loc,column);
      refine_count = 0;
      if (DEBUG_REF_INT) { (void) printf("refining int log: i, testval = %d %d\n",i,testval); fflush(stdout); } 
      }
   if (sieve_array[column][loc] >= testval) 
      {
      if (DEBUG_SUCCESS) 
         { 
         (void) printf("found one at #%d (%d,%d) sum=%d, vs:%d %d\n", i, loc,column, sieve_array[column][loc], 
                         testval,num_total_success); 
         fflush(stdout); 
         }
      sieve_array[column][loc] = GOOD;
      success[num_total_success][LOC_STORE] = loc;
      success[num_total_success][COLUMN_STORE] = column;
      num_total_success++;
      }
   }

if (DEBUG) (void) printf("FINAL SUCCESS COUNT = %d\n",num_total_success);

}   /* end of find_success_locations */


/************************************************************************/
/*                                                                      */
/*   hash function; computes index into list for saving primes when     */ 
/*   resieving                                                          */
/*                                                                      */
/************************************************************************/

int_hash(c,d)
int c,d;

{   /* start of int_hash */
int index, i;

index = (c << SHIFT) + (d);

if (CHECK_HASH)
   {
   if (index > TSIZE || index < 0) (void) printf("Hash fail: %d %d %d\n",index,c,d);
   exit(0);
   }

/*  check to see if occupied */

if (int_faclist[index][0] == 0)  return(index);

/*  if occupied, check for match with c & d */

if (int_faclist[index][1] == c  && int_faclist[index][2] == d) return(index);

/*  collision; check next few locations */

for (i=index+1; i<index+200; i++)       // I hope 100 is enough!
   { 
   if (int_faclist[i][0] == 0) return(i);
   if (int_faclist[i][1] == c  && int_faclist[i][2] == d) return(i);
   }

(void) printf("Integer hash failure at %d %d %d\n",index, c, d);
for (i=index; i<index+20; i++)
   (void) printf("%d %d %d\n", int_faclist[i][0],int_faclist[i][1], int_faclist[i][2]);
exit(0); 

}   /* end of int_hash */


/************************************************************************/
/*                                                                      */
/*   hash function; computes index into list for saving primes when     */ 
/*   resieving                                                          */
/*                                                                      */
/************************************************************************/

alg_hash(c,d)
int c,d;

{   /* start of alg_hash */
int index, i;

if (ALG_HASH_DEBUG) (void) printf("hashing (%d,%d)\n",c,d);

index = (c << SHIFT) + (d);

if (CHECK_HASH)
   {
   if (index > TSIZE || index < 0) (void) printf("Hash fail: %d %d %d\n",index,c,d);
   exit(0);
   }

/*  check to see if occupied */

if (alg_faclist[index][0] == 0)  return(index);

/*  if occupied, check for match with c & d */

if (alg_faclist[index][1] == c  && alg_faclist[index][2] == d) return(index);

/*  collision; check next few locations */

if (ALG_HASH_DEBUG) 
   {
   (void) printf("collision at %d\n",index);
   for (i=index+1; i<index+10; i++) 
      (void) printf("stored: %d %d %d\n",alg_faclist[i][0],alg_faclist[i][1], alg_faclist[i][2]);
   }

for (i=index+1; i<index+200; i++)       // I hope 100 is enough!
   {
   if (alg_faclist[i][0] == 0) return(i);
   if (alg_faclist[i][1] == c  && alg_faclist[i][2] == d) return(i);
   }

(void) printf("Algebraic hash failure at %d %d %d\n",index, c,d);
for (i=index; i<index+20; i++)
   (void) printf("%d %d %d\n", alg_faclist[i][0], alg_faclist[i][1], alg_faclist[i][2]);
exit(0); 

}   /* end of alg_hash */


/************************************************************************/
/*                                                                      */
/*   Routine to reset log values after sieving for projective primes    */ 
/*                                                                      */
/************************************************************************/

void restore_logs()

{   /* start of restore_logs */
int i,index;

for (i=1; i<=proj_int_primes[0]; i++)
   {
   index = proj_int_primes[i];
   int_logp[index] = proj_int_logs[i];
   }

for (i=1; i<=proj_alg_primes[0]; i++)
   {
   index = proj_alg_primes[i];
   alg_logp[index] = proj_alg_logs[i];
   }


}   /* end of restore_logs */


/************************************************************************/
/*                                                                      */
/*   Routine to compute intersection points and slopes that define the  */ 
/*   vector sieve region                                                */
/*  V1 = (a0, b0), V2 = (a1, b1). We need  (-M,0) < e V1 + f V2 < (M,N) */
/*  Comment: YECH!   there has GOT to be a better way!!!!!!!!!!!        */
/*  There are 4 lines, designated as L1, L2, L3, L4.   L2 & L4 go thru  */
/*  the origin. L1 & L2  and  L3 & L4 are parallel.                     */
/*   Note: this is called both by vector sieve and the resiever. We     */
/*   could call it just once and store the results, but it would be a   */
/*   LOT of storage; 5 doubles/ideal                                    */
/*                                                                      */
/************************************************************************/

void prepare_bounds(v1, v2, sign)
int v1[2], v2[2], sign;

{   /* start of prepare_bounds */

double m_over_a1, n_over_b1, b0_over_b1, a0_over_a1;
double a0, a1, b0, b1, one_over_a1, one_over_b1;
double l1and4, l2and3, l1and3, determ, inv_determ;
double stime;
double temp;

if (TIME_STATS) stime = get_time();

a0 = (double)v1[0];
a1 = (double)v2[0];
b0 = (double)v1[1];
b1 = (double)v2[1];

one_over_a1 = 1.0/a1;      
a0_over_a1 = a0 * one_over_a1;
one_over_b1 = 1.0/b1;      
b0_over_b1 = b0 * one_over_b1;

/*  We have a parallelogram.  One point is always (0,0).  f is the   */
/*  vertical axis,  e the horizontal.  Compute emin & emax           */
/*   Also, compute intersections and boundary slopes                 */
/*   Note that determ = p (up to sign) so could precompute: but it   */
/*   would require xtra storage to hold the sign bit                 */

determ = (double)(v1[0] * v2[1] - v2[0] * v1[1]);   
inv_determ = 1.0/determ;                        

if (SHOW_PREP)
   {
   (void) printf("Prep: a0,a1,b0,b1 = %g %g %g %g\n",a0,a1,b0,b1);
   (void) printf("m_over_a1, n_over_b1 = %g %g\n",m_over_a1, n_over_b1);
   (void) printf("a0_over_a1 , b0_over_b1 = %g %g\n",a0_over_a1, b0_over_b1);
   (void) printf("determ, inv = %g %g\n",determ,inv_determ);
   (void) printf("L1and4, L2and3, L1and3 = %g %g %g\n",l1and4, l2and3, l1and3);
   }

if (sign == 1) {

if (a0 > 0  && a1 > 0)
   {
   l1and4 = (double)(MAX_SIEVE * v2[1]) * inv_determ;
   l2and3 = (double)(-NCOLUMNS * v2[0]) * inv_determ;

   emin = iceil(min(l2and3, l1and4));  
   emax = max(l2and3, l1and4);
    lower_midpt = 0.0;
   lm1 = -max(a0_over_a1, b0_over_b1);
   lm2 = -min(a0_over_a1, b0_over_b1);
   lc1 = 0.0;
   lc2 = 0.0;
   }
else if (a0 < 0 && a1 < 0)
   {
   temp = (double)(MAX_SIEVE * v2[1]);
   l1and4 = temp*inv_determ;
   l1and3 = (temp - (double)(NCOLUMNS* v2[0])) * inv_determ;
   m_over_a1 = DMAX_SIEVE * one_over_a1;
   emin = iceil(min(0.0,l1and3));
   emax = max(0.0, l1and3);
   lower_midpt = l1and4;
   lm1 = -max(a0_over_a1, b0_over_b1);
   lm2 = -min(a0_over_a1, b0_over_b1);
   if (a0_over_a1 < b0_over_b1)
      {
      lc1 = 0.0;
      lc2 = m_over_a1;
      }
   else
      {
      lc1 = m_over_a1;
      lc2 = 0.0;
      }
   }
else if (a0 < 0 && a1 > 0)
   {
   l1and4 = (double)(MAX_SIEVE * v2[1]) * inv_determ;
   l2and3 = (double)(-NCOLUMNS * v2[0]) * inv_determ;
   emin = iceil(l1and4);
   emax = l2and3;
   lower_midpt = 0.0;
   lm1 = -b0_over_b1;
   lm2 = -a0_over_a1;
   lc1 = 0.0;
   lc2 = 0.0;
   }
else  // a0 > 0,  a1 < 0
   {
   l1and4 = (double)(MAX_SIEVE * v2[1]) * inv_determ;
   l1and3 = (double)(MAX_SIEVE * v2[1] - NCOLUMNS * v2[0]) * inv_determ;
   m_over_a1 = DMAX_SIEVE * one_over_a1;

   emin = 0.0;
   emax = l1and3;
   lower_midpt = l1and4;
   lm1 = -b0_over_b1; 
   lm2 = -a0_over_a1;
   lc1 = 0.0;
   lc2 = m_over_a1;
   }

}  /* end of sign == 1 */

else {   /* sign == -1 */

if (a0 > 0  && a1 > 0)    
   {
   l1and4 = (double)(-MAX_SIEVE * v2[1]) * inv_determ;
   l1and3 = (double)(-MAX_SIEVE * v2[1] - NCOLUMNS * v2[0]) * inv_determ;
   m_over_a1 = DMAX_SIEVE * one_over_a1;

   emin = iceil(min(l1and3, 0.0));
   emax = max(l1and3, 0.0);
    lower_midpt = l1and4;
   lm1 = -max(a0_over_a1, b0_over_b1);
   lm2 = -min(a0_over_a1, b0_over_b1);
   if (a0_over_a1 < b0_over_b1)
      {
      lc1 = 0.0;
      lc2 = -m_over_a1;
      }
   else
      {
      lc1 = -m_over_a1;
      lc2 = 0.0;
      }
   }
else if (a0 < 0 && a1 < 0)
   {
   l1and4 = (double)(-MAX_SIEVE * v2[1]) * inv_determ;
   l2and3 = (double)(-NCOLUMNS* v2[0]) * inv_determ;
   emin = iceil(min(l1and4, l2and3));
   emax = max(l1and4, l2and3);
   lower_midpt = 0.0;
   lm1 = -max(a0_over_a1, b0_over_b1);
   lm2 = -min(a0_over_a1, b0_over_b1);
   lc1 = 0.0;
   lc2 = 0.0;
   }
else if (a0 < 0 && a1 > 0)
   {
   l1and4 = (double)(-MAX_SIEVE * v2[1]) * inv_determ;
   l1and3 = (double)(-MAX_SIEVE * v2[1] - NCOLUMNS * v2[0]) * inv_determ;
   m_over_a1 = DMAX_SIEVE * one_over_a1;

   emin = 0.0;
   emax = l1and3;
   lower_midpt = l1and4;
   lm1 = -b0_over_b1;
   lm2 = -a0_over_a1;
   lc1 = 0.0;
   lc2 = -m_over_a1;
   }
else // a0 > 0,  a1 < 0
   {
   l1and4 = (double)(-MAX_SIEVE * v2[1]) * inv_determ;
   l2and3 = (double)(-NCOLUMNS * v2[0]) * inv_determ;
   emin = iceil(l1and4);
   emax = l2and3;
   lower_midpt = 0.0;
   lm1 = -b0_over_b1; 
   lm2 = -a0_over_a1;
   lc1 = 0.0;
   lc2 = 0.0;
   }

}  /* end of sign == -1 */

if (TIME_STATS) vector_setup_time += (get_time() - stime);

}   /* end of prepare_bounds */


/************************************************************************/
/*                                                                      */
/*   Routine to sort list of factors; special_q gets inserted out of    */
/*  order, so fix it here.                                              */
/*                                                                      */
/************************************************************************/

void sort_int_factor_list(list)
int list[];

{   /* start of sort_int_fator_list */
int i,j,t1;
int len,min,place;

len = list[0];
for (i=1; i<=len; i++)
   {
   min = list[i];
   place = i;
   for (j=i+1; j<=len; j++)      /* find min */
      {
      if (list[j] < min) {min = list[j]; place = j;}
      }
   t1 = list[i];                /* swap      */
   list[i] = min;
   list[place] = t1;
   }

}   /* end of sort_int_factor_list */



/************************************************************************/
/*                                                                      */
/*   Routine to initialize time statistic variables                     */
/*                                                                      */
/************************************************************************/

void init_time_variables()

{   /* start of init_time_variables */

vector_setup_time  = 0.0;
int_linesieve_time = 0.0;
alg_linesieve_time = 0.0;
even_int_linesieve_time = 0.0;
even_alg_linesieve_time = 0.0;
int_vector_time = 0.0;
alg_vector_time = 0.0;
int_resieve_time = 0.0;
alg_resieve_time = 0.0;
trial_int_time = 0.0;
trial_alg_time = 0.0;
find_startpt_time = 0.0;
alg_scan_time = 0.0;
squfof_time = 0.0;
latred_time = 0.0;
find_success_time = 0.0;  
total_ptime = 0.0; 
invtime = 0.0;

}   /* end of init_time_variables */


/************************************************************************/
/*                                                                      */
/*   Routine to initialize count statistic variables                    */
/*                                                                      */
/************************************************************************/

void init_count_variables()

{   /* start of init_count_variables */

int_bad =  0;
alg_bad =  0;
int_prime =  0;
alg_prime =  0;
int_cofactor =  0;
alg_cofactor =  0;
total_hits =  0;
squfof_fail =  0;
overflow_count =  0;
squfof_calls =  0;
badsplit = 0;
rho_calls =  0;
total_q_successes = 0;
num_alg_success = 0;
relation_count = 0;
not_coprime = 0;

}   /* end of init_count_variables */

/************************************************************************/
/*                                                                      */
/*   Routine to partition sieve loops; for unrolling                    */
/*                                                                      */
/************************************************************************/

void partition_loops()

{   /* start of partition_loops */
int i;

for (i=0; i<size_int_base; i++)        
   {
   if (int_fbase[i] > (MAX_SIEVE >> 3)) {int_8hits = i;  break;}
   }

for (i=int_8hits; i<size_int_base; i++)
   {
   if (int_fbase[i] > (MAX_SIEVE >> 2)) {int_4hits = i; break;}
   }

for (i=int_4hits; i<size_int_base; i++)
   {
   if (int_fbase[i] > (MAX_SIEVE/3)) {int_3hits = i; break;}
   }

for (i=int_3hits; i<size_int_base; i++)
   {
   if (int_fbase[i] > (MAX_SIEVE >> 1)) {int_2hits = i; break;}
   }

int_1hits = size_int_base; 
for (i=int_2hits; i<size_int_base; i++)
   {
   if (int_fbase[i] > (MAX_SIEVE)) {int_1hits = i; break;}
   }

int_linesieve = size_int_base;
for (i=int_1hits; i<size_int_base; i++)
   {
   if (int_fbase[i] > (FUZZY*MAX_SIEVE)) { int_linesieve = i; break; }
   }

for (i=0; i<size_alg_base; i++)
   {
   if (alg_fbase[i] > (MAX_SIEVE >> 3)) {alg_8hits = i; break;}
   }

for (i=alg_8hits; i<size_alg_base; i++)
   {
   if (alg_fbase[i] > (MAX_SIEVE >> 2)) {alg_4hits = i; break;}
   }

for (i=alg_4hits; i<size_alg_base; i++)
   {
   if (alg_fbase[i] > (MAX_SIEVE/3)) {alg_3hits = i; break;}
   }
 
for (i=alg_3hits; i<size_alg_base; i++)
   {
   if (alg_fbase[i] > (MAX_SIEVE >> 1)) {alg_2hits = i; break;}
   }

alg_1hits = size_alg_base;
for (i=alg_2hits; i<size_alg_base; i++)
   {
   if (alg_fbase[i] > (MAX_SIEVE)) {alg_1hits = i; break;}
   }

alg_linesieve = size_alg_base; 
for (i=alg_1hits; i<size_alg_base; i++)
   {
   if (alg_fbase[i] > (FUZZY*MAX_SIEVE)) {alg_linesieve = i; break;}
   }

}   /* end of partition_loops */


/************************************************************************/
/*                                                                      */
/*   Routine to decide whether a reduced lattice is 'bad'. Special q    */
/*  that yield highly skewed parallelograms have low yields. Highly     */
/*  skewed means a vertex angle of less than 8  degrees (8  selected by */
/*   data sampling).   Or if the L2 norm is too  big relative to sqrt(q)*/
/*                                                                      */
/************************************************************************/
reject_this_q()

{   /* start of reject_this_q */
double ang,t1,t2,t3,t4;


#define   BAD_ANGLE   8.0
#define   BAD_RATIO   7.0

t1 = sqrt(v1[0] * v1[0] + v1[1] * v1[1]);
t2 = sqrt(v2[0] * v2[0] + v2[1] * v2[1]);
t3 = atan(-v1[0]/v2[0]);
t4 = atan(-v1[1]/v2[1]);
ang = fabs(180.0 -  (t3-t4)*180.0/3.14159);

if (ang < BAD_ANGLE) return(1);
if (max(t1,t2)/sqrt((double)special_q) > BAD_RATIO) return(1);

if (abs(v1[0]) > (1 << MAX_RED_BITS)) return(1); /* not a good reduced lattice */
if (abs(v1[1]) > (1 << MAX_RED_BITS)) return(1);   
if (abs(v2[0]) > (1 << MAX_RED_BITS)) return(1);
if (abs(v2[1]) > (1 << MAX_RED_BITS)) return(1);

return(0);

}   /* end of reject_this_q */

/************************************************************************/
/*                                                                      */
/*   Routine to computer single precision modular inverses              */
/*                                                                      */
/************************************************************************/

int single_modinv (int a, int modulus)

{ /* start of single_modinv */

int ps1, ps2, parity, dividend, divisor, rem, q, t;
double t1;


if (TIME_STATS) t1 = get_time();
if (a == 1) return 1;
q = modulus / a;
rem = modulus - a * q;
dividend = a;
divisor = rem;
ps1 = q;
ps2 = 1;
parity = ~0;


while (divisor > 1)
   {
   rem = dividend - divisor;
   t = rem - divisor;
   if (t >= 0) {
     q += ps1;
     rem = t;
     t -= divisor;
     if (t >= 0) {
       q += ps1;
       rem = t;
       t -= divisor;
       if (t >= 0) {
         q += ps1;
         rem = t;
         t -= divisor;
         if (t >= 0) {
           q += ps1;
           rem = t;
           t -= divisor;
           if (t >= 0) {
             q += ps1;
             rem = t;
             t -= divisor;
             if (t >= 0) {
               q += ps1;
               rem = t;
               t -= divisor;
               if (t >= 0) {
                 q += ps1;
                 rem = t;
                 t -= divisor;
                 if (t >= 0) {
                   q += ps1;
                   rem = t;
                   if (rem >= divisor) {
                     q = dividend/divisor;
                     rem = dividend - q * divisor;
                     q *= ps1;
                   }}}}}}}}}
   q += ps2;
   parity = ~parity;
   dividend = divisor;
   divisor = rem;
   ps2 = ps1;
   ps1 = q;
   }
 
if (TIME_STATS) invtime += (get_time() - t1);
if (parity == 0) return (ps1);
else return (modulus - ps1);

} /* end of single_modinv */


/************************************************************************/
/*                                                                      */
/*  New routine to do lattice basis reduction via Euclidean algorithm   */
/*  experimental version.....only 3 subts.............................. */
/*                                                                      */
/************************************************************************/

newest_latred_exp(q,r, v1, v2)
int q,r, v1[], v2[];

{   /* start of newest_latred_exp */
int a11, a12, a21, a22, t, rem;
double norm_row1, tmp_norm, dot_prod, norm_row2, test;
int accum1, accum2, t1,t2;
int tmp_a11, tmp_a12;


//  maybe always return biggest norm in v1????

if (r == 1)   
   {          
   v1[0] = -(q >> 1);
   v1[1] = q + v1[0];
   v2[0] = 1;
   v2[1] = 1;
   return(0);
   }

a11 = q;
a21 = r;
a12 = 0;
a22 = 1;

//printf("| %d %d |\n", a11, a12);
//printf("| %d %d |\n\n", a21, a22);

/* we take advantage of the fact that most partial quotients are small  */
/* by doing the division by successive conditional subtraction.         */
/* Note:  last conditional subtraction, or maybe last 2 don't make much */
/* difference.                                                          */

norm_row1 = (double)a11 * (double)a11 + (double)a12 * (double)a12;
norm_row2 = (double)a21 * (double)a21 + (double)a22 * (double)a22;

while(1)
   {
   accum1 = a21;
   accum2 = a22;
   rem = a11 - a21;
   if (rem >= a21)
      {
      accum1 += a21;
      accum2 += a22;
      rem -= a21;
      if (rem >= a21)
         {
         accum1 += a21;
         accum2 += a22;
         rem -= a21;

         if (rem >= a21)
            {
            accum1 += a21;
            accum2 += a22;
            rem -= a21;
            if (rem >= a21)
               {
               t = a11/a21;
               accum1 = t * a21;
               accum2 = t * a22;
               }
            }
         }
      }
   tmp_a11 = a11 - accum1;           // update row 1?
   tmp_a12 = a12 - accum2;

   tmp_norm = (double)tmp_a11 * (double)tmp_a11 + (double)tmp_a12 * (double)tmp_a12;

   if (tmp_norm > norm_row1)  break; // update would increase norm so quit

   t1 = tmp_a11;                     // swap rows
   t2 = tmp_a12;
   a11 = a21;
   a12 = a22;
   a21 = t1;
   a22 = t2;
   norm_row1 = norm_row2;
   norm_row2 = tmp_norm;
   }

/*  Need to do a single Gram-Schmidt row reduction if possible */
/*  check to see if t = 0 before doing the division            */
/*  depending on the dot_prod/row_norm ratio we might adjust   */
/*  either row1 or row2                                        */

dot_prod = (double)a11 * (double)a21 + (double)a12 * (double)a22;

test = fabs(dot_prod + dot_prod);

if (test < norm_row2) 
   {
   if (test < norm_row1)   // no update possible
      {
      v1[0] = a11;
      v1[1] = a12;
      v2[0] = a21;
      v2[1] = a22;
      }
   else
      {
      t = round_double(dot_prod/norm_row1);
      v1[0] = a11;
      v1[1] = a12;
      v2[0] = a21 - t * a11;
      v2[1] = a22 - t * a12;
      }
   }
else 
   {
   t = round_double(dot_prod/norm_row2);
   v1[0] = a11 - t * a21;
   v1[1] = a12 - t * a22;
   v2[0] = a21;
   v2[1] = a22;
   }

/*  make b  (2nd value) positive */

if (v1[1] < 0)
   {
   v1[1] = -v1[1];
   v1[0] = -v1[0];
   }
if (v2[1] < 0)
   {
   v2[1] = -v2[1];
   v2[0] = -v2[0];
   }


//printf("final result:\n");
//printf("| %d %d |\n",   v1[0], v1[1]);
//printf("| %d %d |\n\n", v2[0], v2[1]);

return(0);

}   /* end of newest_latred_exp */


/************************************************************************/
/*																		*/
/*	routine to do single precision divide for numbers with a small		*/
/*  limb count.                                                         */
/*																		*/
/************************************************************************/
 
small_div_single(a,b,result)
unsigned int a[],b,result[];
 
{   /* start of small_div_single */
register unsigned int i,rem,as,first;
unsigned int  x[2];
 
//  note: for 1st div,  r = 0. can elim entirely if b > a[i]
//  do direct single/single division otherwise without calling
//  divrem.  Useful if  size is small.

as = SIZE(a);
first = a[as-1];

if (b > first)		  // no division needed!
   {				  // this assignment might not be needed 
   result[as-1] = 0;  // quotient = 0,  rem is as assigned;
   rem = first;
   SIZE(result) = as-1;
   }
else
   {
   result[as-1] = first/b;
   rem = first - result[as-1] *b;
   SIZE(result) = as;
   } 

for (i=as-2; i>0; i--)
	{
	new_divrem_asm(rem,a[i],b,x);
	result[i] = x[0];
	rem = x[1];
	}


return(rem);

}   /* end of small_div_single */


/************************************************************************/
/*																		*/
/*	compute (a*2^30 + b)/c  and (a*2^30 + b) % c						*/
/*	new assembler version inlined for div_single						*/
/*																		*/
/************************************************************************/


__inline void new_divrem_asm(a,b,c,d)
unsigned int a,b,c,d[2];
 
{	/* start of new_divrem_asm */

/* We could use a double length register fromthe mmx instruction set,	*/
/* however, the emmx instruction must be executedto clean up the		*/
/* FP registers every time we use mmx.  Emmx is a very lengthy			*/
/* instruction compared to what we are doing here.						*/

	_asm {

/*		edx:eax = (a << 30) + tempb;									*/

//	mov		eax,a    old code
//	shl		eax,30
//	mov		edx,a
//	shr		edx,2
//	add		eax,b
//	adc		edx,0

	mov		eax,b
	mov		edx,a
	shl		eax,2
	shrd	eax,edx,2
	shr		edx,2

/*  Now divide a * (2^30) which resides in the register pair: edx:eax	*/
/*  eax = edx:eax / c													*/
/*	edx = edx:eax % c													*/

	mov		ecx,c

/*  d[x] is a pointer (moved ahead here for a little pentium optimization*/
	mov		edi,d
    div		ecx

/*   d[0] = ((a << 30) + b) / c;										*/
	mov		DWORD PTR[edi],eax

/*	 d[1] = ((a << 30) + b) % c;										*/ 
 	mov		DWORD PTR[edi]+4,edx

	} // end _asm

}	/* end of new_divrem_asm */

/********************************/
/* auxiliary routine for latred */ // obsolete;  use round_double on Pentiums.
/********************************/

mrint(double d)
{
if (d<0) return(-(long)(-d+0.5));
return((long)(d+0.5));
}


/********************************/
/* auxiliary routine for debug  */
/********************************/

void check_start_points()
{

int i,err,j;
int c,d,rem;
int alg_value[MPDIM], int_value[MPDIM], temp[MPDIM];


(void) printf("V1, V2 = (%d,%d), (%d,%d)\n",v1[0],v1[1], v2[0], v2[1]);

err = 0;
for (i=NUMSMALLP; i<int_linesieve; i++)
   {
   compute_coeffs(int_startpts[i], 1, &c, &d);
   int_norm(d, c, int_value);
   rem = div_single(int_value, int_fbase[i], temp);
   if (rem != 0)
      {
      for (j=1; j<=proj_int_primes[0]; j++) if (proj_int_primes[j] == i) goto cont1;
      (void) printf("int start pt error at %d: %d %d %d %d\n", i, c, d, int_fbase[i], int_startpts[i]);
      err = 1;
      }
cont1: ;
   }

for (i=NUMSMALLP; i<alg_linesieve; i++)
   {
   compute_coeffs(alg_startpts[i], 1, &c, &d);
   alg_norm(d, c, alg_value);
   rem = div_single(alg_value, alg_fbase[i], temp);
   if (rem != 0)
      {
      for (j=1; j<=proj_alg_primes[0]; j++) if (proj_alg_primes[j] == i) goto cont;
      (void) printf("alg start pt error at %d: %d %d %d %d\n", i, c, d, alg_fbase[i], alg_startpts[i]);
      err = 1;
      }
cont: ;
   }

if (err == 1) exit(0);
}


/****************************/
/* aux routine: show vector   */
/****************************/

show_int_vector(i)
int i;
{
(void) printf("%d %d : %d %d\n",int_v1[i][0],int_v1[i][1], int_v2[i][0], int_v2[i][1]);
return(0);
}
show_alg_vector(i)
int i;
{
(void) printf("%d %d : %d %d\n",alg_v1[i][0],alg_v1[i][1], alg_v2[i][0], alg_v2[i][1]);
return(0);
}


/************************************/
/* aux routine: show sieve region   */
/************************************/

void show_region(p)
int p;

{   /* start of show_region */

(void) printf("Region for %d\n",p);
(void) printf("emin,emax = %g %g\n",emin,emax);
(void) printf("Upper:  %g*x + %g :  %g*x + %g, midpt = %g\n",um1,uc1,um2,uc2,upper_midpt);
(void) printf("Lower:  %g*x + %g :  %g*x + %g, midpt = %g\n",lm1,lc1,lm2,lc2,lower_midpt);
(void) printf("=========================================================================\n");

}   /* end of show_region */

/************************/
/*  show sieve areas   */
/************************/

show_areas()
{
int i;
int int_upper_tot, int_lower_tot, alg_upper_tot, alg_lower_tot;

int_upper_tot = 0;
int_lower_tot = 0;
alg_upper_tot = 0;
alg_lower_tot = 0;


//for (i=int_1hits; i<size_int_base; i++)
for (i=int_linesieve; i<300; i++)
   {
   (void) printf("int i, upper, lower: %d %d %d\n",i,int_upper_area[i],int_lower_area[i]);
   int_upper_tot += int_upper_area[i];
   int_lower_tot += int_lower_area[i];
   }

//for (i=alg_linesieve; i<size_alg_base; i++)
for (i=alg_linesieve; i<300; i++)
   {
   (void) printf("alg i, upper, lower: %d %d %d\n",i,alg_upper_area[i],alg_lower_area[i]);
   alg_upper_tot += alg_upper_area[i];
   alg_lower_tot += alg_lower_area[i];
   }

(void) printf("int upper,lower total = %d %d\n",int_upper_tot, int_lower_tot);
(void) printf("alg upper,lower total = %d %d\n",alg_upper_tot, alg_lower_tot);

return(0);
}

/************************/
/* show hash array      */
/************************/

show_hash(base)
int base;
{
int i;

(void) printf("showing hash arrays\n");
for (i=base; i<base+15; i++)
   {
   (void) printf("%d %d %d\n",int_faclist[i][0], int_faclist[i][1], int_faclist[i][2]);
   }
printf("\n******************\n");
for (i=base; i<base+15; i++)
   {
   (void) printf("%d %d %d\n",alg_faclist[i][0], alg_faclist[i][1], alg_faclist[i][2]);
   }
printf("\n");
return(0);
}


/**************************************/
/* show quality of lattice reductions */
/**************************************/

void display_reductions()
{
int q;
double t1,t2,t3,t4,ang;
int special_q, special_q_root,det;

printf("special q stats\n");
printf("  q      root      v1       v2      norm ratio    Angle\n");
printf("=============================================================\n");

for (q=1000; q<size_int_base; q++)
   {
   special_q = int_fbase[q];
   special_q_root = int_fbase[q] - int_roots[q]; 
   reduce_lattice(special_q, special_q_root, v1, v2);
   t1 = sqrt(v1[0] * v1[0] + v1[1] * v1[1]);
   t2 = sqrt(v2[0] * v2[0] + v2[1] * v2[1]);
   t3 = atan(-v1[0]/v2[0]);
    t4 = atan(-v1[1]/v2[1]);
    ang = fabs(180.0 -  (t3-t4)*180.0/3.14159);
    det = v1[0] * v2[1] - v2[0] * v1[1];
    (void) printf("%9d %9d (%5d,%5d) (%5d,%5d) %12.6lf %12.6lf\n", special_q, special_q_root,
         v1[0],v1[1], v2[0],v2[1], max(t1,t2)/sqrt((double)special_q), ang);
//   (void) printf("%9d %12.6lf %12.6lf\n",special_q,  max(t1,t2)/(double)special_q, ang);

   }
}

/***************************************/
/*   Checks/displays skewed lattices   */
/***************************************/

void check_lattice(str, index, iv1, iv2)
int index, iv1[2], iv2[2];
char str[];
{
double t1,t2,t3,t4,ang,t5;

t1 = sqrt(iv1[0] * iv1[0] + iv1[1] * iv1[1]);
t2 = sqrt(iv2[0] * iv2[0] + iv2[1] * iv2[1]);
t3 = atan(-iv1[0]/iv2[0]);
t4 = atan(-iv1[1]/iv2[1]);
ang = fabs(180.0 -  (t3-t4)*180.0/3.14159);
t5 = max(t1,t2)/sqrt((double)int_fbase[index]);

if (ang < 10.0 || t5 > 5.0)
   {
   (void) printf("bad %s sub_lattice at %d (%d,%d) (%d,%d): %12.6f %12.6f\n",str,
               int_fbase[index], iv1[0],iv1[1],iv2[0],iv2[1], ang, t5);
   show_region(int_fbase[index]);
   }
}


/*******************************************************/
/* Routine to compare two different lattice reductions */
/*******************************************************/


cmp_latred(q,r,v1, v2,  v1a, v2a)
int q,r,v1[], v2[], v1a[], v2a[];
{
double n1, n2;

if (v1[0] == v1a[0]  && v1[1] == v1a[1] && v2[0] == v2a[0] && v2[1] == v2a[1]) return(0);
if (v1[0] == v2a[0]  && v1[1] == v2a[1] && v2[0] == v1a[0] && v2[1] == v1a[1])
   {
   printf("swapped at %d %d\n",q,r);
   n1 = (double) v1[0] * v1[0] + (double) v1[1] * v1[1];
   n2 = (double) v2[0] * v2[0] + (double) v2[1] * v2[1];
   if (n1 < n2) printf("2nd norm greater\n");
//   printf("v1, v2 norms = %lf %lf\n", n1, n2);
   return(0);
   }

printf("mismatch %d %d\n",q,r);
printf("%d %d\n",v1[0],v1[1]);
printf("%d %d\n",v2[0],v2[1]);
printf("==========================\n");
printf("%d %d\n",v1a[0],v1a[1]);
printf("%d %d\n",v2a[0],v2a[1]);
return(0);
}



/**********************************/
/*  small, inline assembler       */
/*  routine to round dble to      */
/*  nearest int                   */
/**********************************/

__inline int round_double(double a)

{   /* start of round_double */
int d;
_asm 
	{
	fld a
	fistp d
	};
return d;
}   /* end of round_double */


/************************************************************************/
/*                                                                      */
/*  New routine to do lattice basis reduction via Euclidean algorithm   */
/*  experimental version.....only 3 subts.............................. */
/*                                                                      */
/************************************************************************/

latred(q,r, v1, v2)
int q,r, v1[], v2[];

{   /* start of newest_latred_exp */
int a11, a12, a21, a22, t, rem;
double norm_row1, tmp_norm, dot_prod, norm_row2, test;
int accum1, accum2;
int tmp_a11, tmp_a12;

if (r == 1)   
   {          
   v1[0] = -(q >> 1);
   v1[1] = q + v1[0];
   v2[0] = 1;
   v2[1] = 1;
   return(0);
   }
a11 = q;
a21 = r;
a12 = 0;
a22 = 1;

/* we take advantage of the fact that most partial quotients are small  */
/* by doing the division by successive conditional subtraction.         */
/* Note:  last conditional subtraction, or maybe last 2 don't make much */
/* difference.                                                          */

norm_row2 = (double)a21 * a21 + (double)a22 * a22;

while(1)
   {
   accum1 = a21;
   accum2 = a22;
   rem = a11 - a21;
   if (rem >= a21)
      {
      accum1 += a21;
      accum2 += a22;
      rem -= a21;
      if (rem >= a21)
         {
         accum1 += a21;
         accum2 += a22;
         rem -= a21;
         if (rem >= a21)
            {
            accum1 += a21;
            accum2 += a22;
            rem -= a21;
            if (rem >= a21)
               {
               accum1 += a21;
               accum2 += a22;
               rem -= a21;
               if (rem >= a21)
                  {
                  accum1 += a21;
                  accum2 += a22;
                  rem -= a21;
                  if (rem >= a21)
                     {
                     t = a11/a21;
                     accum1 = t * a21;
                     accum2 = t * a22;
                     }
                  }
               }
            }
         }
      }
   tmp_a11 = a11 - accum1;           // update row 1?
   tmp_a12 = a12 - accum2;

   tmp_norm = (double)tmp_a11 * tmp_a11 + (double)tmp_a12 * tmp_a12;

   if (tmp_norm > norm_row2)  break; // update would increase norm so quit

   a11 = a21;
   a12 = a22;
   a21 = tmp_a11;
   a22 = tmp_a12;
   norm_row2 = tmp_norm;
   }

/*  Need to do a single Gram-Schmidt row reduction if possible */
/*  check to see if t = 0 before doing the division            */
/*  depending on the dot_prod/row_norm ratio we might adjust   */
/*  either row1 or row2                                        */

dot_prod = (double)a11 * a21 + (double)a12 * a22;
norm_row1 = (double) a11 * a11 + (double) a12 * a12;

test = 2 * dot_prod;
if (test < 0) test = -test;

if (test > norm_row2) 
   { 
   t = round_double(dot_prod/norm_row2);
   v1[0] = a11 - t * a21;
   v1[1] = a12 - t * a22;
   v2[0] = a21;
   v2[1] = a22;
   }
else   
   {
   if (test < norm_row1)   
      { 
      v1[0] = a11;
      v1[1] = a12;
      v2[0] = a21;
      v2[1] = a22;
      }
   else
      { 
      t = round_double(dot_prod/norm_row1);
      v1[0] = a11;
      v1[1] = a12;
      v2[0] = a21 - t * a11;
      v2[1] = a22 - t * a12;
      }
   }

/*  make b  (2nd value) positive */

if (v1[1] < 0)
   {
   v1[1] = -v1[1];
   v1[0] = -v1[0];
   }
if (v2[1] < 0)
   {
   v2[1] = -v2[1];
   v2[0] = -v2[0];
   }

return(0);

}   /* end of newest_latred_exp */
