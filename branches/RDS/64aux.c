#include <stdio.h>

#define	RADIX			1073741824
#define	BASE			(double)1073741824.0
#define MPDIM			150
#define	QUAD			5
#define	TRIPLE			4
#define	DOUBLE			3
#define	SINGLE			2
#define	SIZE(x)			x[0]
#define	FIRST(x)		x[x[0] - 1]
#define	LASTDIGIT(x)	x[1]
#define	DEBUG			0
#define	DEBUGS			0
char numstr[300];
 
void mult_single_prec(),mul_dble_by_single();
extern void write_num(),mpcopy(),mul_single(),subt(),mult(),divv(),exit();
extern int div_single(),mpcmp();

/************************************************************************/
/*																		*/
/*	Low level primitives for 64 bit arithmetic; portable version		*/
/*																		*/
/************************************************************************/
 
 
/************************************************************************/
/*																		*/
/*	mod_mul(a,b,c);  returns a * b mod c; a,b,c, unsigned 32 bit		*/
/*	integers															*/
/*																		*/
/************************************************************************/
mod_mult(a,b,c)
unsigned int a,b,c;

{
register unsigned int quo,rem;

quo = (unsigned int)((double)a * (double)b/(double) c + .5);
rem = a*b - quo*c;
if ((int)rem <= -1) rem += c;
return(rem);
}

/************************************************************************/
/*																		*/
/*	compute (a * b + c) / 2^30     and  (a * b + c) % 2^30				*/
/*																		*/
/************************************************************************/
 
void muladd(a,b,c,d)
unsigned int a,b,c,d[2];

{
register double temp;

temp = (double)a * (double)b + (double)c;
d[1] = (unsigned int) (temp/BASE + .5);
d[0] = a*b + c - d[1]*RADIX;
if ((int)d[0] <= -1)
   {
   d[0] += RADIX;
   d[1]--;
   }
}

/************************************************************************/
/*																		*/
/*	compute (a * b ) / 2^30     and  (a * b ) % 2^30					*/
/*																		*/
/************************************************************************/
 
void radmul(a,b,d)
unsigned int a,b,d[2];

{
register unsigned int rem,quo;
register double temp;

temp = (double)a * (double)b;
quo = (unsigned int) (temp/BASE + .5);
rem = a*b  - quo*RADIX;
if ((int)rem <= -1)
   {
   rem += RADIX;
   quo--;
   }
d[0] = rem;
d[1] = quo;
}
 
/************************************************************************/
/*																		*/
/*	compute (a*2^30 + b)/c  and (a*2^30 + b) % c						*/
/*																		*/
/************************************************************************/

void divrem(a,b,c,d)
unsigned int a,b,c,d[2];

{
register double temp;

temp = (double)a * BASE + (double)b;
d[0] = (unsigned int) (temp/(double) c + .5);
d[1] = a*RADIX + b - d[0]*c;
if ((int)d[1] <= -1)
   {
   d[1] += c;
   d[0]--;
   }
}

/************************************************************************/
/*																		*/
/*	divide quad prec number by double prec number						*/
/*																		*/
/************************************************************************/


void div4by2(a,b,q,r)
int a[],b[],q[],r[];

{   /* start of div4by2 */
int quo,rem,temp[MPDIM],temp1[MPDIM],upper,lower,result_case;
double dividend,divisor,test;

#define	BASE2	1152921504606846976.00
#define	BASE3	BASE*BASE2

if (DEBUG)
   {
   (void) printf("in div4, size = %d %d\n",SIZE(a),SIZE(b));
   (void) printf("dividing "); write_num(a,numstr);
   (void) printf("by ");
   write_num(b,numstr);
   if (SIZE(a) > QUAD) exit(0); 
   }

switch(SIZE(a))							/* numerator computation	*/
   {
   case QUAD:							/* quad precision			*/
      dividend = a[1] + a[2]*BASE + a[3]*BASE2 + a[4]*BASE3;
      break;

   case TRIPLE:							/* triple precision			*/
	  dividend = a[1] + a[2]*BASE + a[3]*BASE2;
	  break;

   case DOUBLE:							/* double precision			*/
   	  dividend = a[1] + a[2]*BASE;
	  break;

   case SINGLE:							/* single precision			*/
   	  dividend = a[1];
	  break;
   }

switch(SIZE(b))							/* denominator computation	*/
   {
   case TRIPLE:							/* triple precision			*/
      divisor = b[1] + b[2]*BASE + b[3]*BASE2;
      break;

   case DOUBLE:							/* double precision			*/
	  divisor = b[1] + b[2]*BASE;
	  break;

   case SINGLE:							/* single precision			*/
      rem = div_single(a, LASTDIGIT(b),q);
	  SIZE(r) = SINGLE;
	  FIRST(r) = rem;
	  return;
	  break;
   }

if (DEBUG) (void) printf("dividend,divisor = %20.5f %20.5f\n",dividend,divisor);

test = (dividend/divisor + .01);

if (DEBUG) { (void) printf("test val = %20.5f\n",test); fflush(stdout); }

result_case = SINGLE;
if (test > BASE)  result_case = DOUBLE;
if (test > BASE2) result_case = TRIPLE;

if (DEBUG) { (void) printf("result_case = %d\n",result_case); fflush(stdout); }

switch(result_case)
   {
   case SINGLE:							/* quotient is single precision	*/
      quo = (int)test;
	  if (quo == 0)
	     {
		 mpcopy(a,r);
		 SIZE(q) = SINGLE;
		 FIRST(q) = 0;
		 return;
		 }
	  mul_dble_by_single(b,quo,temp);
	  while(mpcmp(temp,a) == 1)		/* should be at most 1 too big	*/
	     {
		 quo--;
		 if (quo == 0)
		    {
		    mpcopy(a,r);
		    SIZE(q) = SINGLE;
		    FIRST(q) = 0;
		    return;
		    }
		 subt(temp,b,temp);
		 }
	  subt(a,temp,r);
	  SIZE(q) = SINGLE;
	  FIRST(q) = quo;
	  return;
      break;

   case DOUBLE:							/* quotient is double precision	*/
   	  upper = (int)(test/BASE);
	  if (DEBUG) (void) printf("upper = %d\n",upper);
	  q[2] = upper;
	  q[1] = 0;
	  SIZE(q) = DOUBLE;
	  mult(b,q,temp);
	  while (mpcmp(temp,a) == 1)
		{
		q[2]--;
		subt(temp,b,temp);
		}
	  subt(a,temp,temp1);
	  if (SIZE(temp1) == TRIPLE)
		  lower = (int) ((temp1[3]*BASE2 + temp1[2]*BASE + temp1[1])/divisor + 1.0);
	  else if (SIZE(temp1) == DOUBLE)
		  lower = (int) ((temp1[2]*BASE + temp1[1])/divisor + 1.0);
	  else
		  lower = (int) (temp1[1]/divisor + 1.0);
	  q[1] = lower;
	  mult(b,q,temp);
	  while (mpcmp(temp,a) == 1)
	     {
		 q[1]--;
		 subt(temp,b,temp);
		 }
	  subt(a,temp,r);
	  return;
      break;

   case TRIPLE:							/* quotient is triple precision	*/
	  divv(a,b,q,r);					/* don't bother; use full div	*/
	  if (DEBUG) { (void) printf("case 4: rem = "); write_num(r,numstr); }
	  return;
	  break;
   }

}   /* end of div4by 2 */


/************************************************************************/
/*																		*/
/*	square (at most) double precision number							*/
/*																		*/
/************************************************************************/

void square_dble_prec(a, answer)
int a[],answer[];

{   /* start of square_dble_prec */


int x1[2],x21[2],x22[2];

if (DEBUG)
   {
   (void) printf("Squaring, size = %d\n",SIZE(a));
   (void) printf("number = "); write_num(a,numstr);      
   if (SIZE(a) > DOUBLE) 
      {
      (void) printf("attempt to square non double precision number\n");
	  exit(0);
	  }
   }

mult_single_prec(a[1],a[1],x1);
if (SIZE(a) == SINGLE)
   {
   answer[1] = x1[0];
   answer[2] = x1[1];
   if (x1[1] == 0) SIZE(answer) = SINGLE;
   else SIZE(answer) = DOUBLE;
   return;
   }

SIZE(answer) = QUAD;
mult_single_prec(a[1],a[2],x21);
mult_single_prec(a[2],a[2],x22);

answer[1] = x1[0];
answer[2] = x1[1];
answer[3] = x22[0];
answer[4] = x22[1];

if (DEBUG)
   {
   (void) printf("a = %d %d %d\n",a[0],a[1],a[2]);
   (void) printf("squaring double: x1, x21, x22 = %d %d %d %d %d %d\n",x1[0],x1[1],x21[0],x21[1],x22[0],x22[1]);
   }
answer[2] += x21[0];
if (answer[2] >= RADIX)
	{
	answer[2] -= RADIX;
	answer[3] += 1;							/* carry		*/
	}
answer[2] += x21[0];
if (answer[2] >= RADIX)
	{
	answer[2] -= RADIX;
	answer[3] += 1;							/* carry		*/
	}

answer[3] += x21[1];
if (answer[3] >= RADIX)
	{
	answer[3] -= RADIX;
	answer[4] += 1;							/* carry		*/
	}
answer[3] += x21[1];
if (answer[3] >= RADIX)
	{
	answer[3] -= RADIX;
	answer[4] += 1;							/* carry		*/
	}

if (answer[4] == 0) SIZE(answer)--;

if (DEBUG)
   {
   printf("sq dble: size, answer = %d ",SIZE(answer));
   write_num(answer,numstr);
   }

}   /* end of square_dble_prec */

/************************************************************************/
/*																		*/
/*	double length product of two single precision numbers				*/
/*  note: same as muladd without the add								*/
/*																		*/
/************************************************************************/

void mult_single_prec(a,b,answer)
unsigned int a,b, answer[2];

{   /* start of mult_single_prec */

register double temp;

temp = (double)a * (double)b;
answer[1] = (unsigned int) (temp/BASE + .5);
answer[0] = a*b - answer[1]*RADIX;
if ((int)answer[0] <= -1)
   {
   answer[0] += RADIX;
   answer[1]--;
   }
}


/************************************************************************/
/*																		*/
/*	multiply (at most) double prec number by single						*/
/*																		*/
/************************************************************************/

void mul_dble_by_single(a, b, answer)
int a[],b,answer[];

{   /* start of mul_dble_by_single */

int x1[2],x21[2];

mult_single_prec(a[1],b,x1);
if (SIZE(a) == SINGLE)
   {
   answer[1] = x1[0];
   answer[2] = x1[1];
   if (x1[1] == 0) SIZE(answer) = SINGLE;
   else SIZE(answer) = DOUBLE;
   return;
   }

mult_single_prec(a[2],b,x21);

answer[1] = x1[0];
answer[2] = x1[1];
answer[3] = x21[1];
SIZE(answer) = SIZE(a) + 1;

answer[2] += x21[0];
if (answer[2] >= RADIX)
	{
	answer[2] -= RADIX;
	answer[3] += 1;							/* carry		*/
	}

if (answer[3] == 0) SIZE(answer)--;

}   /* end of mul_dble_by_single */
