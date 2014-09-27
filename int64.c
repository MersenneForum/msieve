/********************************************************************************/
/*																				*/
/*	This file contains low level 64-bit arithmetic routines.					*/
/*	There are two versions: a C version, using __int64's as the data type		*/
/*	And an assembler version.	C version: routine name ends in I; asm ends 	*/
/*	in _asm.  Note that there are also floating point versions of these same	*/
/*	routines.  (in 64aux.c)														*/
/********************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mpdim.h"
 
void mul_single_prec(unsigned int a, unsigned int b, unsigned int *d),
	 mul_dble_by_single_asm(unsigned int *n1, unsigned int b, unsigned int *n2),
	 mul_single_prec_asm(unsigned int a, unsigned int b, unsigned int *d),
	 mul_dble_by_single_asm(unsigned int *n1, unsigned int b, unsigned int *n2),
	 square_dble_prec_asm(unsigned int *n1, unsigned int *n2),
	 divrem_asm(unsigned int a, unsigned int b, unsigned int c, unsigned int *d),
	 muladd_asm(unsigned int a, unsigned int b, unsigned int c, unsigned int *d);

int mod_multI(unsigned int a, unsigned int b, unsigned int c),
	mod_mult_asm(unsigned int a, unsigned int b, unsigned int c),
	mod_mult_add_asm(unsigned int a, unsigned int b, unsigned int c, unsigned int d);


void muladdI(unsigned int a, unsigned int b, unsigned int c, unsigned int *d),
	 divremI(unsigned int a, unsigned int b, unsigned int c, unsigned int *d),
	 mul_single_precI(unsigned int a, unsigned int b, unsigned int *d),
	 mul_dble_by_singleI(unsigned int *n1, unsigned int b, unsigned int *n2);

extern void write_num(int *num, char *str),
			mpcopy(int *n1, int *n2),
			mul_single(int *num, int b, int *res),
			subt(int *n1, int *n2, int *n3),
			mult(int *n1, int *n2, int *n3),
			divv(int *n1, int *n2, int *q, int *r);

extern int div_single(int *n1, int b, int *n2),
		   mpcmp(int *n1, int *n2);

/************************************************************************/
/*																		*/
/*	mod_mul(a,b,c);  returns a * b mod c; a,b,c, unsigned 32 bit		*/
/*	integers; assumed  c> a  and c > b									*/
/*																		*/
/************************************************************************/
mod_multI(a,b,c)
unsigned int a,b,c;

{
register unsigned __int64 ab;
unsigned int rem;

ab = (unsigned __int64) a  *  (unsigned __int64) b;
rem = (unsigned int)ab % c;
return(rem);
}

/************************************************************************/
/*																		*/
/*	compute (a * b + c) / 2^30     and  (a * b + c) % 2^30				*/
/*																		*/
/************************************************************************/
 
void muladdI(a,b,c,d)
unsigned int a,b,c,d[2];

{
register unsigned __int64 temp;

temp = (unsigned __int64)a * (unsigned __int64)b + (unsigned __int64)c;
d[1] = (unsigned int)temp >> 30;
d[0] = (unsigned int)temp & 0x3fffffff;
}


 
/************************************************************************/
/*																		*/
/*	compute (a*2^30 + b)/c  and (a*2^30 + b) % c						*/
/*																		*/
/************************************************************************/

void divremI(a,b,c,d)
unsigned int a,b,c,d[2];

{
unsigned __int64 temp;

temp = ((unsigned __int64)a << 30) + (unsigned __int64)b;
d[1] = (unsigned int)temp % c;
d[0] = (unsigned int)temp/c;


}

/************************************************************************/
/*																		*/
/*	multiply single prec number by single; same as muladd w/o the add	*/
/*																		*/
/************************************************************************/

void mul_single_precI(a,b,answer)
unsigned int a,b, answer[2];

{   /* start of mul_single_precI */

register unsigned __int64 temp;

temp = (unsigned __int64)a * (unsigned __int64)b;
answer[1] = (unsigned int)temp >> 30;
answer[0] = (unsigned int)temp & 0x3fffffff;

}


/************************************************************************/
/*																		*/
/*	multiply (at most) double prec number by single						*/
/*																		*/
/************************************************************************/

void mul_dble_by_singleI(a, b, answer)
unsigned int a[],b,answer[];

{   /* start of mul_dble_by_single */

int x1[2],x21[2];

mul_single_precI(a[1],b,x1);
if (SIZE(a) == SINGLE)
   {
   answer[1] = x1[0];
   answer[2] = x1[1];
   if (x1[1] == 0) SIZE(answer) = SINGLE;
   else SIZE(answer) = DOUBLE;
   return;
   }

mul_single_precI(a[2],b,x21);

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


/************************************************************************/
/*																		*/
/*	compute (a*2^30 + b)/c  and (a*2^30 + b) % c						*/
/*	assembler version													*/
/*																		*/
/************************************************************************/


void divrem_asm(a,b,c,d)
unsigned int a,b,c,d[2];
 
{	/* start of divrem_asm */

/* We could use a double length register fromthe mmx instruction set,	*/
/* however, the emmx instruction must be executedto clean up the		*/
/* FP registers every time we use mmx.  Emmx is a very lengthy			*/
/* instruction compared to what we are doing here.						*/

	_asm {

/*		edx:eax = (a << 30) + tempb;									*/

	mov		eax,a
	shl		eax,30
	mov		edx,a
	shr		edx,2
	and		edx,0x3fffffff

	add		eax,b
	adc		edx,0

/*  Now divide a * (2^30) which resides in the register pair: edx:eax	*/
/*  eax = edx:eax / c													*/
/*	edx = edx:eax % c													*/

	mov		ecx,c

/*  d[x] is a pointer (moved ahead here for a pentium optimization      */
	mov		edi,d
    div		ecx

/*   d[0] = ((a << 30) + b) / c;										*/
	mov		DWORD PTR[edi],eax

/*	 d[1] = ((a << 30) + b) % c;										*/ 
 	mov		DWORD PTR[edi]+4,edx

	} // end _asm

}	/* end of divrem_asm */
 

/************************************************************************/
/*																		*/
/*	muladd(a,b,c,d);  returns (a * b + c)/2^30  and						*/
/*	(a * b + c) % 2^30; a,b,c, unsigned 32 bit							*/
/*	integers															*/
/*																		*/
/*	assembler version													*/
/************************************************************************/
void muladd_asm(a,b,c,d)
unsigned int a,b,c,d[2];
 
{   /* start of muladd_asm */

_asm 
	{
	mov eax, a			/* compute a b + c				*/
    mul	b				/* edx:eax  holds product		*/
    add eax, c			/* add c to low word			*/
	adc	edx, 0			/* carry bit					*/

	mov edi, d			/* get address of d				*/
	mov ecx, eax		/* get copy of low word			*/
	
	shl edx, 2			/* get top 30 bits of edx:eax	*/
	shr	ecx, 30
	add edx, ecx
	mov DWORD PTR [edi]+4, edx	/* top 30 in d[1]		*/

	and eax, 1073741823	/* low 30 bits of edx:eax		*/
	mov DWORD PTR [edi], eax	/* store in d[0]		*/
	}

}   /* end of muladd_asm */


/************************************************************************/
/*																		*/
/*	mod_mul(a,b,c);  returns a * b mod c; a,b,c, unsigned 32 bit		*/
/*	integers															*/
/*	assembler version													*/
/*																		*/
/************************************************************************/
mod_mult_asm(a,b,c)
unsigned int a,b,c;

{   /* start of mod_mult_asm */

_asm
	{
	mov eax, a			/* compute a*b  in edx:eax		*/
	mul b
    
	div c				/* div by c, remainder in edx	*/
	mov eax , edx		/* return the remainder in eax	*/
	}
//  lint complains because of no return statement

}   /* end of mod_mul_asm */


/************************************************************************/
/*																		*/
/*	multiply single prec number by single; same as muladd w/o the add	*/
/*																		*/
/************************************************************************/
void mul_single_prec_asm(a,b,d)
unsigned int a,b,d[];

{   /* start of mul_single_prec_asm */

_asm 
	{
	mov eax, a			/* compute a b					*/
    mul	b				/* edx:eax  holds product		*/

	mov edi, d			/* get address of d				*/
	mov ecx, eax		/* get copy of low word			*/
	
	shl edx, 2			/* get top 30 bits of edx:eax	*/
	shr	ecx, 30
	add edx, ecx
	mov DWORD PTR [edi]+4, edx	/* top 30 in d[1]		*/

	and eax, 1073741823	/* low 30 bits of edx:eax		*/
	mov DWORD PTR [edi], eax	/* store in d[0]		*/
	}

}   /* start of mul_single_prec_asm */



/************************************************************************/
/*																		*/
/*	multiply (at most) double prec number by single						*/
/*	This is in C, but it calls the _asm version of mul_single_prec		*/
/*																		*/
/************************************************************************/

void mul_dble_by_single_asm(a, b, answer)
unsigned int a[],b,answer[];

{   /* start of mul_dble_by_single_asm */

int x1[2],x21[2];

mul_single_prec_asm(a[1],b,x1);
if (SIZE(a) == SINGLE)
   {
   answer[1] = x1[0];
   answer[2] = x1[1];
   if (x1[1] == 0) SIZE(answer) = SINGLE;
   else SIZE(answer) = DOUBLE;
   return;
   }

mul_single_prec_asm(a[2],b,x21);

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

}   /* end of mul_dble_by_single_asm */


/************************************************************************/
/*																		*/
/*	square (at most) double precision number							*/
/*	assembler version;  in C but calls mul_single_prec_asm				*/
/*																		*/
/************************************************************************/

void square_dble_prec_asm(a, answer)
unsigned int a[],answer[];

{   /* start of square_dble_prec_asm */


int x1[2],x21[2],x22[2];

mul_single_prec_asm(a[1],a[1],x1);
if (SIZE(a) == SINGLE)
   {
   answer[1] = x1[0];
   answer[2] = x1[1];
   if (x1[1] == 0) SIZE(answer) = SINGLE;
   else SIZE(answer) = DOUBLE;
   return;
   }

SIZE(answer) = QUAD;
mul_single_prec_asm(a[1],a[2],x21);
mul_single_prec_asm(a[2],a[2],x22);

answer[1] = x1[0];
answer[2] = x1[1];
answer[3] = x22[0];
answer[4] = x22[1];

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

}   /* end of square_dble_prec_asm */


/************************************************************************/
/*																		*/
/*	divide quad prec number by double prec number; old version			*/
/*	--> this is slow and badly designed									*/
/*	use new version instead												*/
/************************************************************************/


void old_div4by2_asm(a,b,q,r)
int a[],b[],q[],r[];

{   /* start of old_div4by2_asm */
int quo,rem,temp[MPDIM],temp1[MPDIM],upper,lower,result_case;
double dividend,divisor,test;

#define	BASE2	1152921504606846976.00
#define	BASE3	BASE*BASE2

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

test = (dividend/divisor + .01);

result_case = SINGLE;
if (test > BASE)  result_case = DOUBLE;
if (test > BASE2) result_case = TRIPLE;

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
	  mul_dble_by_single_asm(b,quo,temp);
	  while(mpcmp(temp,a) == 1)			/* should be at most 1 too big	*/
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
	  return;
	  break;
   }

}   /* end of old_div4by2_asm */


/********************************************************************/
/*																	*/
/*	Routine to perform division of (at most) a quad precision		*/
/*	integer by (at most) a triple precision integer;  Most of the	*/
/*	time,  this will be a quad or triple divided by a double		*/
/*	yielding a double, so the routine focusses on optimizing		*/
/*	these cases. b can never be more than double precision			*/
/*	only the remainder gets returned; the quotient is  no interest	*/
/*	The smallest remainder might not get returned.					*/
/*																	*/
/********************************************************************/

void new_div4by2_asm(a,b,rem)
unsigned int a[],b[],rem[];

{   /* start of new_div4by2_asm */

#define	SECOND(x)	x[SIZE(x) - 2]
#define	THIRD(x)	x[SIZE(x) - 3]
#define	BASE2	1152921504606846976.00			/* 2^60			*/
#define	BASE3	BASE*BASE2						/* 2^90			*/

int qsize,r;
double numer,denom,quotient,upper;
unsigned int ans[QUAD],temp[MPDIM];

qsize = SIZE(a) - SIZE(b) + 2;
if (FIRST(b) > FIRST(a)) qsize--;
if (SIZE(b) == SINGLE)					/* b is single precision	*/
	{									/* so call div_single (rare)*/
	r = div_single(a,FIRST(b), ans);
	SIZE(rem) = SINGLE;
	FIRST(rem) = r;
	return;
	}

/*		OK.  b, The divisor is double precsion						*/

switch(SIZE(a))							/* a is dble, triple or quad*/
	{
	case QUAD:							/* should not be possible	*/
		if (qsize == TRIPLE)
	      {
		  divv(a,b,temp,rem);
		  return;
		  }

/*	OK,  we have QUAD with a DOUBLE prec quotient; it doesn't		*/
/*  matter if the remainder isn't as small as possible, as long as	*/
/*  it is in the right residue class and it is reasonably small		*/
/*  thus,  the upper half of the quotient must be right, but the	*/
/*  lower half only needs to be a good estimate: use doubles		*/
	    
		numer = (double)FIRST(a) * BASE3 + (double)SECOND(a) * BASE2
			  + (double)THIRD(a) * BASE + LASTDIGIT(a);
		denom = (double)FIRST(b) * BASE + LASTDIGIT(b);
		quotient = numer/denom;
		upper = quotient/BASE;			/* split quo into upper		*/
		SIZE(ans) = DOUBLE;				/* and lower 30 bits		*/
		FIRST(ans) = (int)upper;
		LASTDIGIT(ans) = (int)( quotient - BASE * FIRST(ans)); 
		mult(b,ans,temp);
		while (mpcmp(temp,a) == 1) subt(temp,b,temp);  
		subt(a,temp,rem);
		while (mpcmp(rem,b) == 1) subt(rem,b,rem);
		return;
	    break;


	case TRIPLE:
		if (qsize == SINGLE)			/* quotient is single prec	*/
		   {
		   numer= (double)FIRST(a) * BASE2 + (double)SECOND(a) * BASE + (double)THIRD(a);
		   denom = (double)FIRST(b) * BASE + LASTDIGIT(b);
		   quotient = numer/denom;
		   mul_single(b, (int)quotient, temp);
		   while (mpcmp(temp,a) == 1) subt(temp,b,temp);
		   subt(a,temp,rem);
		   return;
		   }

/*	OK,  a is triple precision,  the quotient is double precision	*/
/*	This is the most frequent case									*/

		numer = (double)FIRST(a) * BASE2 + (double)SECOND(a) * BASE + (double)THIRD(a);
		denom = (double)FIRST(b) * BASE + LASTDIGIT(b);
		quotient = numer/denom;
		upper = quotient/BASE;			/* split quo into upper		*/
		SIZE(ans) = DOUBLE;				/* and lower 30 bits		*/
		FIRST(ans) = (int)upper;
		LASTDIGIT(ans) = (int)( quotient - BASE * FIRST(ans)); 
		mult(b,ans,temp);
		while (mpcmp(temp,a) == 1) subt(temp,b,temp);
		subt(a,temp,rem);
		break;

/*	a  is double precision,  quotient is single precision			*/

	case DOUBLE:							/* very rare			*/
		if (qsize == 1)
			{
			mpcopy(a,rem);
			return;
			}
		numer = (double)FIRST(a) * BASE + (double)SECOND(a);
		denom = (double)FIRST(b) * BASE + (double)SECOND(b);
		r = (int)(numer/denom);
		mul_single(b, r, temp);
		while (mpcmp(temp,a) == 1) subt(temp,b,temp);
		subt(a, temp, rem);
		break;

	case SINGLE:							/* should not happen	*/
		mpcopy(a,rem);
		return;
		break;
	}

}   /* end of new_div4by2_asm */



/************************************************************************/
/*																		*/
/*	mod_mult_add(a,b,c,d);  returns (ab + c) mod d; a,b,c,d  32 bit		*/
/*	integers;  returns positive value									*/
/*																		*/
/************************************************************************/

modmultadd_asm(a,b,c,d)
int a,b,c,d;

{   /* start of modmultadd_asm */

_asm
	{
	mov eax, a			/* compute a*b  in edx:eax		*/
	imul b
//	
	mov ecx,c
	mov ebx,c
	sar ebx,31			/* sign extend, then add c		*/
	add eax,ecx
	adc edx,ebx			/* carry bit to high order word	*/
//	vs:
//	mov ecx, eax
//	mov ebx, edx
//	mov eax, c
//	cdq
//	add ecx, eax
//	adc ebx, edx	
//	mov eax, ecx
//	mov edx, ebx
//  add eax, c			/* add c to low order word		*/
//  adc edx, 0			/* carry bit					*/

	idiv d				/* div by d, remainder in edx	*/
	cmp edx, 0			/* is it > 0?					*/
	jge SHORT $LABEL1	
	add edx, d			/* if < 0 add the divisor		*/
$LABEL1:
	mov eax , edx		/* return the remainder in eax	*/
	}

}   /* end of modmultadd_asm */



