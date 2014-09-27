 
/************************************************************************/
/*																		*/
/*	This program solves the polynomial congruences:						*/
/*	f(x) = x^d - c  ==  0 mod p  for primes p in the NFS factor 		*/
/*	base. It uses Berlekamp's method. If there is a single root			*/
/*	then computing GCD(f(x), x^p-1 - 1) will find it. Else,				*/
/*	GCD(f(x), (x-s)^(p-1/2) -1 ) has a 50-50 chance of splitting		*/
/*	f(x) into smaller factors. [all done mod p of course]				*/
/*																		*/
/*																		*/
/************************************************************************/
 
#include <stdio.h>
#include <math.h>
 
#define	DEBUG		0
#define	RADIX		1073741824
#define	FBASE		2000000
#define	NBITS		30
#define	HALF_RADIX	32768
#define	MPDIM		150
#define LOGTWO		0.69314
#define NEGCHAR		-5
#define SIEVEMPDIM	500000
#define SMALLPRIME	10
#define	SIGN(x)		(x < 0 ? -1 : 1)

/************************************************************************/
/*																		*/
/*		G L O B A L     V A R I A B L E S								*/
/*																		*/
/************************************************************************/

typedef int *row;

FILE *pcpairs;

int smallprimes[] = {256,243,625,343,121,169,289,361,529,841,961 };

int int_fbase[FBASE],
     int_roots[FBASE],
     alg_roots[FBASE],
     alg_fbase[FBASE],
     pcoeff[10];
 
int expon,
     sfb1,
     sfb2,
     base,
     degree,
     poly_constant,
     mixed_flag,
     sym,
     sign,
     inv_flag;
 
int temp[MPDIM], BIGM[MPDIM];

int store[MPDIM], mem[MPDIM], shift, ppp, pp2, df, notfound, f[MPDIM];
int bbbb[2], bppb[2], aa[MPDIM], ee[MPDIM], nn[MPDIM],power,pirme,pirmem;
int exex[MPDIM];
 
char filetype[12],
     invert[12],
     mixed[12],
     numstr[300];

 
static int single_modinv(int a, int modulus)
{   /* start of single_modinv */
 
register int ps,ps1,ps2=1,parity,dividend,divisor,rem,q;
 
q = modulus/a;
rem = modulus - a*q;
dividend = a;
divisor = rem;
ps1 = q;
parity = 0;
 
while (rem != 0)
	{
	q = 1;
	rem = dividend - divisor;
	if (rem >= divisor)
	   {
	   q = 2;
	   rem -= divisor;
	   if (rem >= divisor)
	      {
	      q = 3;
	      rem -= divisor;
	      if (rem >= divisor)
		{
		q = 4;
		rem -= divisor;
		if (rem >= divisor)
		   {
		   q = dividend/divisor;
		   rem = dividend - q*divisor;
		   }
		}
	      }
	   }
	ps = q*ps1 + ps2;
	ps2 = ps1;
	ps1 = ps;
	dividend = divisor;
	divisor = rem;
	parity = 1 - parity;
	}
 
if (parity == 0) return(ps2);
else return(modulus-ps2);

}   /* end of single_modinv */


/************************************************************************/
/*																		*/
/*			MAIN														*/
/*																		*/
/************************************************************************/

main() 
 
{   /* start of main */
int i;
 
setup();
makepcpairs();
 
pcpairs = fopen("PAIRS","r");
for (i=0; i<sfb2; i++)
   {
   (void) fscanf(pcpairs,"%ld %ld",&alg_fbase[i],&alg_roots[i]);
   }
 
if (strcmp(filetype,"binary") == 0) dump_binary();
else dump_ascii();

(void) printf("done\n");
 
}   /* end of main */
 
/************************************************************************/
/*																		*/
/*	Set up necessary constants, compute primes, call other stuff		*/
/*																		*/
/************************************************************************/
 
setup() 
 
{   /* start of setup */
register int i,tempinv;
int Amodp, Bmodp, prime, Acoeff[MPDIM], Bcoeff[MPDIM];

//zstart();  Why is this needed?
 
(void) printf("give base and lhs exponent\n"); 
(void) scanf("%ld %ld",&base,&expon);  
if (base < 0)				/* general root; not base^expon	*/
   get_num(BIGM, numstr);
 
(void) printf("give integer factor base size\n"); 
(void) scanf("%ld",&sfb1);  
 
(void) printf("give degree, polynomial constant\n");
(void) scanf("%ld %ld",&degree, &poly_constant);
if (degree < 0)
   {
   degree = -degree;
   for (i=degree; i>=0; i--) scanf("%d ",&pcoeff[i]);
   }
else
   {
   if (poly_constant < 0)
	{
	sign = -1;
	poly_constant = -poly_constant;
	}
   else sign = 1;
   }
 
(void) printf("input rhs factor base size\n");
(void) scanf("%ld",&sfb2);
 
 
(void) printf("give max a (as multiple of %ld)\n",SIEVEMPDIM); 
(void) scanf("%ld",&sym);  
(void) printf("sieve ranges from %ld to %ld\n",-sym*SIEVEMPDIM,sym*SIEVEMPDIM);
 
(void) printf("input file type\n");
(void) scanf("%s",filetype);
 
inv_flag = 0;
(void) printf("Invert roots?\n");
(void) scanf("%s",invert);
if (strcmp(invert,"invert") == 0) inv_flag = 1;

mixed_flag = 0;
(void) printf("Mixed LHS Root?\n");
(void) scanf("%s",mixed);
if (strcmp(mixed,"mixedroot") == 0)
   {
   mixed_flag = 1;
   get_num(Acoeff,numstr);
   get_num(Bcoeff,numstr);
   }
 
/*   The left hand size of the congurence is base^expon. Compute this		*/
 
if (base > 0)
   cpower(base,expon,BIGM);		/* BIGM = base^expon		*/
 
/*	Now get the integer factor base and the roots [starting		*/
/*	positions for sieving].						*/


printf("get int fbase\n");
for (i=0; i<SMALLPRIME; i++) 
   {
   int_fbase[i] = zpnext();
   prime = smallprimes[i];
   if (mixed_flag == 0)
      int_roots[i] = div_single(BIGM,prime,temp);
   else
      {
	  Amodp = div_single(Acoeff,prime,temp);
	  Bmodp = div_single(Bcoeff,prime,temp);
	  if (Amodp != 0)
	     {
	     tempinv = single_modinv(Amodp,prime);
	     int_roots[i] = mod_mult(tempinv,Bmodp,prime);
		 }
	  else
	     int_roots[i] = 0;
	  }
   }
 
for (; i<sfb1; i++) 
   {
   int_fbase[i] = zpnext();
   prime = int_fbase[i];
   if (mixed_flag == 0)
      int_roots[i] = div_single(BIGM,prime,temp);
   else
      {
	  Amodp = div_single(Acoeff,prime,temp);
	  Bmodp = div_single(Bcoeff,prime,temp);
	  if (Amodp != 0)
	     {
	     tempinv = single_modinv(Amodp,prime);
	     int_roots[i] = mod_mult(tempinv,Bmodp,prime);
		 }
	  else
	     int_roots[i] = 0;
	  }
   }	
 
(void) printf("first prime is %ld\n",int_fbase[0]); 
(void) printf("last prime is %ld\n",int_fbase[sfb1-1]); 
(void) fflush(stdout); 
 
(void) printf("start doing (p,c) pairs\n"); 
(void) fflush(stdout); 
 
}   /* end of setup */

/************************************************************************/
/*																		*/
/*	Auxiliary routines; simple tasks or call other routines				*/
/*																		*/
/************************************************************************/
//int zzsexp(a,e)
//int a,e; 
//{
//aa[1]=a;
//ee[1]=e;
//zexp(aa,ee,nn,exex);
//return(exex[1]);
//}


int modprod(a,b,p)
int a,b,p;
{
int temp;
int sign;

temp = mod_mult(abs(a),abs(b),p);
sign = SIGN(a)*SIGN(b);
return(sign*temp);
}
/************************************************************************/
/*																		*/
/*	Routine to make monic a polynomial. This is done by multiplying		*/
/*	the coefficients by the inverse of the lead coefficient mod p		*/
/*																		*/
/************************************************************************/
 
char make_monic(f,df,p)
int f[],df,p;
{
register int i, inv;
 
if (f[df] < 0) for (i=df; i>=0; i--) f[i]*= -1; 
if (f[df] == 1) return(0);
if (f[df] != 0) 
   {
   inv = sinv(f[df],p);
   if (modprod(inv,f[df],p) != 1) 
	{
	(void) printf("leading coef %ld non invertible mod %ld\n", f[df],p);
	return(1);
	}
   f[df] = 1;
   for (i=df-1; i>=0; i--) 
   	if (f[i] > 0) f[i] = modprod(f[i],inv,p);
   	else if (f[i] < 0)
	   if ((f[i] = -modprod(-f[i],inv,p)) < 0) f[i]+=p;
   return(0);
   }
 
(void) printf("wrong call to make_monic\n");
return(1);					 /* zero leading coefficient */
}

/************************************************************************/
/*																		*/
/*																		*/
/*	Routine to multiply two polynomials together modulo a third			*/
/*	and all done mod p.													*/
/*																		*/
/************************************************************************/
 
polprod(f,g,rres,mp,dmp,p)
int f[],g[],rres[],mp[],dmp,p; 
{
register int i,j,gi,memdm,dm = dmp-1;
int res[MPDIM], gtop = dm, top = dm;
 
while ((gtop >= 0) && (g[gtop] == 0)) gtop --;
 
if (gtop < 0) 
   {
   for (i=dm; i>=0; i--) rres[i] = 0; 
   return(0);
   }
 
if (mp[dmp] != 1) (void) make_monic(mp,dmp,p);
while ((top > 0) && (mp[top] == 0)) top--;
if (p < HALF_RADIX) 
   {
   if (gi = g[0]) 
	for (i=dm; i>=0; i--) res[i] = (gi*f[i]) % p;
   else
	for (i=dm; i>=0; i--) res[i] = 0;
 
   for (i=dm; i>=0; i--) mem[i] = f[i];
   for (i=1; i<=gtop; i++) 
	{
	memdm = mem[dm];
	for (j=dm; j>top; j--) mem[j] = mem[j-1];
	for (; j>0; j--)
	   if ((mem[j] = (mem[j-1]-memdm*mp[j]) % p) < 0) mem[j] += p;
	if ((mem[0] = (-memdm*mp[0]) % p) < 0) mem[0] += p;
	if (gi = g[i]) for (j=dm; j>=0; j--) res[j] = (res[j]+gi*mem[j]) % p;
	}
   }
else 
   {
   if (gi = g[0]) 
	for (i=dm; i>=0; i--) res[i] = modprod(gi,f[i],p);
   else
	for (i=dm; i>=0; i--) res[i] = 0;
   for (i=dm; i>=0; i--) mem[i] = f[i];
   for (i=1; i <= gtop; i++) 
	{
	memdm = mem[dm];
	for (j=dm; j>top; j--) mem[j] = mem[j-1];
	for (; j>0; j--)
		{
	     if ((mem[j]=(mem[j-1]-modprod(memdm,mp[j],p)) % p) < 0) mem[j]+=p;
		}
	if ((mem[0] = (-modprod(memdm,mp[0],p)) % p) < 0) mem[0] += p;
	if (gi =  g[i]) for (j=dm; j>=0; j--)
		{
		res[j] = (res[j]+modprod(gi,mem[j],p)) % p;
		}
	}
   }
 
for (i=dm; i>=0; i--) rres[i] = res[i]; 
return(0);
}

/************************************************************************/
/*																		*/
/*	Routine to exponentiate a polynomial (x) mod p						*/
/*	Actually used to compute x^p-1 or (x-s)^(p-1)/2						*/
/*	all done mod p; the answer is a polynomial							*/
/*																		*/
/************************************************************************/
xhighmp(result,offset,f,xdf,power,prime)
int f[],xdf,power,prime,result[],offset; 
{
int sq[MPDIM];
register int i;
 
for (i=1;i<=xdf;i++) result[i]=0;
result[0]=1;
 
for (i=2;i<=xdf;i++) sq[i]=0;
sq[1]=1;
sq[0]=offset;
 
while (power > 0) 
   {
   if (power & 1)
	{
	(void) polprod(result,sq,result,f,xdf,prime);
	}
   if (power>>=1)
	{
	 (void) polprod(sq,sq,sq,f,xdf,prime);
	}
   }
}

/************************************************************************/
/*																		*/
/*	Routine to check whether a polynomial is zero						*/
/*																		*/
/************************************************************************/
char polnotzero(f,xdf,p)
int f[],xdf,p;
{
register int i;
 
for (i=0; i<=xdf; i++) if (f[i] % p) return(1); 
return(0);
}

/************************************************************************/
/*																		*/
/*	Routine to determine the degree of a polynomial						*/
/*																		*/
/************************************************************************/
int dgree(f,xdf)
int f[],xdf;
{
register int i=xdf;
while ((i > 0) && (f[i] == 0)) i--;
return(i);
}

/************************************************************************/
/*																		*/
/*	Routine to return GCD of two polynomials mod p						*/
/*																		*/
/************************************************************************/
 
int poly_gcd(f,g,xdf,dg,res,p)
int f[],g[],xdf,dg,res[],p; 
{
register int i,flc,dgg, dff, dif;
int ff[MPDIM], gg[MPDIM];
 
if (xdf >= dg) 
   {
   dgg=dg; dff=xdf;
   for (i=xdf; i>=0; i--) ff[i] = f[i]; 
   for (i=dg; i>=0; i--) gg[i] = g[i]; 
   }
else
   {
   dgg=xdf; dff=dg;
   for (i=dg; i>=0; i--) ff[i]=g[i]; 
   for (i=xdf; i>=0; i--) gg[i]=f[i]; 
   }
 
while (polnotzero(gg,dgg,p)) 
   {
   if ( make_monic(gg,dgg,p) == 1) return(0);
   for (dif=dff-dgg; dif>=0; dif--) 
	{
	if (flc = ff[dff]) 
	   {
	   for (i=dgg-1; i>=0; i--)
		if ((ff[i+dif]-=modprod(flc,gg[i],p)) < 0) ff[i+dif]+=p;
	   }
	dff--;
	}
   for (i=dff; i>=0; i--) mem[i] = ff[i]; 
   for (i=dgg; i>=0; i--) ff[i] = gg[i]; 
   dif=dgree(mem,dff); dff=dgg; dgg=dif;
   for (i=dgg; i>=0; i--) gg[i] = mem[i]; 
   }
 
for (i=dff; i>=0; i--) res[i] = ff[i];  return(dff); 
 
}

/************************************************************************/
/*																		*/
/*	Routine to search for a root when poly splits completely			*/
/*	mod p. Uses Berlekamp's method.										*/
/*																		*/
/************************************************************************/
 
search(shi,ff,dff,lev,ppp)
int shi,ff[],dff,lev,ppp; 
{
int xpow[MPDIM], g[MPDIM], dxpow, dg, ndg, locshi, notsplit=1;

for (locshi=shi; notsplit>0; locshi++) 
   {
   shift = locshi;
   xhighmp(xpow,shift,ff,dff,pp2,ppp);
   if (DEBUG) printf("search: ppp,pp2, shift = %d %d %d\n",ppp,pp2,locshi);
   xpow[0] -= 1;
   if (xpow[0] < 0) xpow[0] += ppp;
   dxpow=dgree(xpow,dff);
   dg = poly_gcd(ff,xpow,dff,dxpow,g,ppp);
   if ((dg > 0) && (dg < dff)) 
	{
	notsplit = 0;
	if (dg == 1) 
	   {
	   notfound--;
	   store[notfound] = ppp-g[0];
	   }
	else 
	   {
	   search(locshi+1,g,dg,lev+1,ppp);
	   }
	xpow[0]+=2;
	if (xpow[0] >= ppp) xpow[0] -= ppp;
	dxpow = dgree(xpow,dff);
	ndg = poly_gcd(ff,xpow,dff,dxpow,g,ppp);
	if (ndg == 1) 
	   {
	   notfound--;
	   store[notfound] = ppp-g[0];
	   }
	else if (ndg > 1)
	   {
	   search(locshi+1,g,ndg,lev+1,ppp);
	   }
	if (ndg+dg-dff) 
	   {
	   notfound--;
	   if (ndg = locshi % ppp) store[notfound]=ppp-ndg;
	   else store[notfound] = 0;
	   }
	}
   }
}

/************************************************************************/
/*																		*/
/*	Routine to compute Integer side factor base and solve for			*/
/*	polynomial roots mod those primes									*/
/*																		*/
/************************************************************************/
makepcpairs() 
 
{   /* start of makepcpairs */
int j, rem, modrem, inverse, result,flucnt=100, counter = 0;
int xpow[MPDIM], g[MPDIM], dxpow, dg, ndg, count, notsplit=1;
 
zpstart();
 
if ((pcpairs = fopen("PAIRS","r")) == NULL) 
   {
   (void) printf("PAIRS must be made, may take a while\n");
   (void) fflush(stdout);
   pcpairs = fopen("PAIRS","w");
 
   aa[0]=1;
   bbbb[0]=1;
   bppb[0]=1;
   nn[0]=1;
   ee[0]=1;
   power = degree;
   df = power;
   pirme = 2;
   for (j=power; j>=0; j--) { f[j] = pcoeff[j]; printf("%d ",f[j]); }
   printf("\n");
   fflush(stdout);
   while (counter < sfb2) 
		{
		for (j=power; j>=0; j--) f[j] = pcoeff[j];
		if ((pcoeff[power] % pirme) == 0)
			{
			(void) fprintf(pcpairs,"%d 0\n",pirme);
			pirme = zpnext();  printf("proj; next prime %d\n",pirme);
			continue;
			}
		nn[1] = pirme;
		pp2 = (pirme-1)/2;
		if (pirme == 2) dg = root2(f,df);
		else
			{ 
			xhighmp(xpow, 0, f, df, (pirme-1), pirme);
			xpow[0]--; 
			if (polnotzero(xpow,df,pirme))
				{ 
				if (xpow[0] < 0) xpow[0] += pirme;
				dxpow = dgree(xpow,df);
				dg = poly_gcd(f,xpow,df,dxpow,g,pirme); 
				}
		else
			{
			dg = df;
			for (j=0; j<=df; j++) g[j] = f[j];
			}
		notfound = dg;

		if (dg == 0)
			{
			if ((f[0] % pirme) == 0)
				{
				dg = 1;
				store[0] = 0;
				}
			}
		else if (dg == 1)
			{
			store[0] = -(g[0] % pirme);
			if (store[0] < 0) store[0] += pirme;
			}
		else search(0L, g, dg, 0L, pirme);
		}
		counter += dg;
		for (j=0; j<dg; j++)
			{
			if (inv_flag == 1)
				{
				if (store[j] != 0 && store[j] != pirme)
				store[j] = single_modinv(store[j],pirme);
				(void) fprintf(pcpairs,"%ld %ld\n",pirme,store[j]);
				(void) fflush(pcpairs);
				}
			else
				{
				(void) fprintf(pcpairs,"%ld %ld\n",pirme,store[j]);
				(void) fflush(pcpairs);
				}
			}
		pirme = zpnext(); 
		if (pirme == 3) pirme = zpnext(); 
		}
   }
else 
   {
   (void) printf("existing file PAIRS will be used, should go fast\n");
   (void) fflush(stdout);
   }
 
}   /* end of makepcpairs */

/************************************************************************/
/*																		*/
/*		Routine to write data back out in ascii							*/
/*																		*/
/************************************************************************/
 
dump_ascii()

{   /* start of dump_ascii */
int i;
FILE *binpc;
 
binpc  =  fopen("ASCII_BASE","w");
(void) fprintf(binpc,"%d\n",base);
(void) fprintf(binpc,"%d\n",poly_constant*sign);
(void) fprintf(binpc,"%d\n",expon);
(void) fprintf(binpc,"%d\n",degree);
(void) fprintf(binpc,"%d\n",sfb1);
(void) fprintf(binpc,"%d\n",sfb2);
 
for (i=0; i<sfb1; i++)
   (void) fprintf(binpc,"%d %d\n",int_fbase[i],int_roots[i]);
 
for (i=0; i<sfb2; i++)
   (void) fprintf(binpc,"%d %d\n",alg_fbase[i],alg_roots[i]);
 
(void) fclose(binpc);
 

}   /* end of dump_ascii */
 
/************************************************************************/
/*																		*/
/*	Routine to create binary version of factor base file				*/
/*																		*/
/************************************************************************/
dump_binary()

{   /* start of dump_binary */
FILE *binpc;
binpc  =  fopen("NFS_BASE","w+b");
(void) fwrite((char*)&base,sizeof(int),1,binpc);
poly_constant = poly_constant*sign;
(void) fwrite((char*)&poly_constant,sizeof(int),1,binpc);
(void) fwrite((char*)&expon,sizeof(int),1,binpc);
(void) fwrite((char*)&degree,sizeof(int),1,binpc);
(void) fwrite((char*)&sfb1,sizeof(int),1,binpc);
(void) fwrite((char*)&sfb2,sizeof(int),1,binpc);
(void) fwrite((char*)int_fbase,sizeof(int),(int)sfb1,binpc);
(void) fwrite((char*)int_roots,sizeof(int),(int)sfb1,binpc);
(void) fwrite((char*)alg_fbase,sizeof(int),(int)sfb2,binpc);
(void) fwrite((char*)alg_roots,sizeof(int),(int)sfb2,binpc);
(void) fclose(binpc);
 
}   /* end of dump_binary */
 
root2(f,df)
int f[],df;
{
int i,count = 0;
int parity = 0;

count = 0;
parity = (f[0] & 1);
if (parity == 0) {store[0] = 0; count++;}
parity = 0;
for (i=0; i<=df; i++) parity += (f[i] & 1);
if ((parity & 1) == 0) {store[count] = 1; count++;}

return(count);
}

dump(str,poly,deg)
char str[];
int poly[],deg;
{
int i;
printf("%s ",str);
for (i=deg; i>=0; i--) printf("%d ",poly[i]);
printf("\n");
fflush(stdout);
}
