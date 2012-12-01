/*--------------------------------------------------------------------
This source distribution is placed in the public domain by its author,
Jason Papadopoulos. You may use it for any purpose, free of charge,
without having to notify anyone. I disclaim any responsibility for any
errors.

Optionally, please be nice and tell me if you find this source to be
useful. Again optionally, if you add to the functionality present here
please consider making those additions public too, so that others may 
benefit from your work.	

$Id$
--------------------------------------------------------------------*/

#if defined(__CUDACC__) && !defined(CUDA_INTRINSICS_H)
#define CUDA_INTRINSICS_H

#ifdef __cplusplus
extern "C"
{
#endif

typedef int int32;
typedef unsigned int uint32;
typedef unsigned long long uint64;
typedef long long int64;

/* 128-bit integers */

typedef struct {
	uint32 w[4];
} uint128;

/*------------------- Low-level functions ------------------------------*/

__device__ void
accum3(uint32 &a0, uint32 &a1, uint32 &a2,
	uint32 b0, uint32 b1) {

	asm("add.cc.u32 %0, %0, %3;   /* inline */   \n\t"
	    "addc.cc.u32 %1, %1, %4;   /* inline */   \n\t"
	    "addc.u32 %2, %2, %5;   /* inline */   \n\t"
		: "+r"(a0), "+r"(a1), "+r"(a2)
		: "r"(b0), "r"(b1), "r"(0) );
}

__device__ void
accum3_shift(uint32 &a0, uint32 &a1, uint32 &a2,
	uint32 b0, uint32 b1) {

	asm("add.cc.u32 %0, %1, %3;   /* inline */   \n\t"
	    "addc.cc.u32 %1, %2, %4;   /* inline */   \n\t"
	    "addc.u32 %2, %5, %5;   /* inline */   \n\t"
		: "=r"(a0), "+r"(a1), "+r"(a2)
		: "r"(b0), "r"(b1), "r"(0) );
}

/*----------------- Squaring ----------------------------------------*/

__device__ uint64 
wide_sqr32(uint32 a)
{
	uint32 a0, a1;

	asm("{ .reg .u64 %dprod; \n\t"
	    "mul.wide.u32 %dprod, %2, %2; \n\t"
	    "cvt.u32.u64 %0, %dprod;      \n\t"
	    "shr.u64 %dprod, %dprod, 32;  \n\t"
	    "cvt.u32.u64 %1, %dprod;      \n\t"
	    "}                   \n\t"
	    : "=r"(a0), "=r"(a1)
	    : "r"(a));

	return (uint64)a1 << 32 | a0;
}

/* -------------------- Modular subtraction ------------------------*/

__device__ uint32 
modsub32(uint32 a, uint32 b, uint32 p) 
{
	uint32 r;

	asm("{  \n\t"
	    ".reg .pred %pborrow;           \n\t"
	    ".reg .u32 %borrow;           \n\t"
	    "mov.b32 %borrow, 0;           \n\t"
	    "sub.cc.u32 %0, %1, %2;        \n\t"
	    "subc.u32 %borrow, %borrow, 0; \n\t"
	    "setp.ne.u32 %pborrow, %borrow, 0;  \n\t"
	    "@%pborrow add.u32 %0, %0, %3; \n\t"
	    "} \n\t"
	    : "=r"(r) : "r"(a), "r"(b), "r"(p) );

	return r;
}

__device__ uint64 
modsub64(uint64 a, uint64 b, uint64 p) 
{
	uint32 r0, r1;
	uint32 a0 = (uint32)a;
	uint32 a1 = (uint32)(a >> 32);
	uint32 b0 = (uint32)b;
	uint32 b1 = (uint32)(b >> 32);
	uint32 p0 = (uint32)p;
	uint32 p1 = (uint32)(p >> 32);

	asm("{  \n\t"
	    ".reg .pred %pborrow;           \n\t"
	    ".reg .u32 %borrow;           \n\t"
	    "mov.b32 %borrow, 0;           \n\t"
	    "sub.cc.u32 %0, %2, %4;        \n\t"
	    "subc.cc.u32 %1, %3, %5;        \n\t"
	    "subc.u32 %borrow, %borrow, 0; \n\t"
	    "setp.ne.u32 %pborrow, %borrow, 0;  \n\t"
	    "@%pborrow add.cc.u32 %0, %0, %6; \n\t"
	    "@%pborrow addc.u32 %1, %1, %7; \n\t"
	    "} \n\t"
	    : "=r"(r0), "=r"(r1)
	    : "r"(a0), "r"(a1), 
	      "r"(b0), "r"(b1), 
	      "r"(p0), "r"(p1));

	return (uint64)r1 << 32 | r0;
}

/*------------------------------- GCD --------------------------------*/
__device__  uint32
gcd64(uint64 x, uint64 y) {

	/* assumes x and y are odd and nonzero */

	uint64 u = x;
	uint64 v = y;

	do {
		v >>= __ffsll(v) - 1;

		x = ullmin(u, v);
		y = ullmax(u, v);
		u = x;
		v = y - x;
	} while (v != 0);

	return u;
}

/*-------------------------- Modular inverse -------------------------*/

__device__ uint32 
modinv32(uint32 a, uint32 p) {

	uint32 ps1, ps2, dividend, divisor, rem, q, t;
	uint32 parity;

	q = 1; rem = a; dividend = p; divisor = a;
	ps1 = 1; ps2 = 0; parity = 0;

	while (divisor > 1) {
		rem = dividend - divisor;
		t = rem - divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t;
		if (rem >= divisor) {
			q = dividend / divisor;
			rem = dividend - q * divisor;
			q *= ps1;
		} } } } } } } } }

		q += ps2;
		parity = ~parity;
		dividend = divisor;
		divisor = rem;
		ps2 = ps1;
		ps1 = q;
	}
	
	if (parity == 0)
		return ps1;
	else
		return p - ps1;
}

__device__ uint64 
modinv64(uint64 a, uint64 p) {

	uint64 ps1, ps2, dividend, divisor, rem, q, t;
	uint32 parity;

	q = 1; rem = a; dividend = p; divisor = a;
	ps1 = 1; ps2 = 0; parity = 0;

	while (divisor > 1) {
		rem = dividend - divisor;
		t = rem - divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t; t -= divisor;
		if (rem >= divisor) { q += ps1; rem = t;
		if (rem >= divisor) {
			q = dividend / divisor;
			rem = dividend - q * divisor;
			q *= ps1;
		} } } } } } } } }

		q += ps2;
		parity = ~parity;
		dividend = divisor;
		divisor = rem;
		ps2 = ps1;
		ps1 = q;
	}
	
	if (parity == 0)
		return ps1;
	else
		return p - ps1;
}

/*------------------- Montgomery arithmetic --------------------------*/
__device__ uint32 
montmul32(uint32 a, uint32 b,
		uint32 n, uint32 w) {

	uint32 acc0, acc1, acc2 = 0;
	uint32 q;
	uint32 prod_lo, prod_hi;

	acc0 = a * b;
	acc1 = __umulhi(a, b);

	q = acc0 * w;

	prod_lo = q * n;
	prod_hi = __umulhi(q, n);

	accum3(acc0, acc1, acc2, prod_lo, prod_hi);

	if (acc2 || acc1 >= n)
		return acc1 - n;
	else
		return acc1;
}

__device__ uint64 
montmul64(uint64 a, uint64 b,
		uint64 n, uint32 w) {

	uint32 a0 = (uint32)a;
	uint32 a1 = (uint32)(a >> 32);
	uint32 b0 = (uint32)b;
	uint32 b1 = (uint32)(b >> 32);
	uint32 n0 = (uint32)n;
	uint32 n1 = (uint32)(n >> 32);
	uint32 acc0, acc1, acc2 = 0;
	uint32 q0, q1;
	uint32 prod_lo, prod_hi;
	uint64 r;

	acc0 = a0 * b0;
	acc1 = __umulhi(a0, b0);
	q0 = acc0 * w;
	prod_lo = q0 * n0;
	prod_hi = __umulhi(q0, n0);
	accum3(acc0, acc1, acc2, prod_lo, prod_hi);

	prod_lo = a0 * b1;
	prod_hi = __umulhi(a0, b1);
	accum3_shift(acc0, acc1, acc2, prod_lo, prod_hi);
	prod_lo = a1 * b0;
	prod_hi = __umulhi(a1, b0);
	accum3(acc0, acc1, acc2, prod_lo, prod_hi);
	prod_lo = q0 * n1;
	prod_hi = __umulhi(q0, n1);
	accum3(acc0, acc1, acc2, prod_lo, prod_hi);
	q1 = acc0 * w;
	prod_lo = q1 * n0;
	prod_hi = __umulhi(q1, n0);
	accum3(acc0, acc1, acc2, prod_lo, prod_hi);

	prod_lo = a1 * b1;
	prod_hi = __umulhi(a1, b1);
	accum3_shift(acc0, acc1, acc2, prod_lo, prod_hi);
	prod_lo = q1 * n1;
	prod_hi = __umulhi(q1, n1);
	accum3(acc0, acc1, acc2, prod_lo, prod_hi);

	r = (uint64)acc1 << 32 | acc0;
	if (acc2 || r >= n)
		return r - n;
	else
		return r;
}

/*------------------ Initializing Montgomery arithmetic -----------------*/
__device__ uint32 
montmul32_w(uint32 n) {

	uint32 res = 2 + n;
	res = res * (2 + n * res);
	res = res * (2 + n * res);
	res = res * (2 + n * res);
	return res * (2 + n * res);
}

__device__ uint32 
montmul32_r(uint32 n) {

	uint32 r0 = ((uint64)1 << 63) % n;
	uint32 r1;

	r1 = r0 + r0;

	if (r1 < r0)
		r1 -= n;

	return modsub32(r1, n, n);
}

__device__ uint64 
montmul64_r(uint64 n, uint32 w) {

	uint32 shift;
	uint32 i;
	uint64 shifted_n;
	uint64 res;

	shift = __clzll(n);
	shifted_n = n << shift;
	res = -shifted_n;

	for (i = 64 - shift; i < 72; i++) {
		if (res >> 63)
			res = res + res - shifted_n;
		else
			res = res + res;

		if (res >= shifted_n)
			res -= shifted_n;
	}

	res = res >> shift;
	res = montmul64(res, res, n, w);
	res = montmul64(res, res, n, w);
	return montmul64(res, res, n, w);
}

/*-------------------------- Mod -------------------------------------*/

	/* the following routines implement algorithms
	   from Ernst W. Mayer's paper "Effecient Long
	   Division via Montgomery Multiply" */

__device__ uint32
modwinv_32(uint32 n) { /* compute 1/n mod 2**32 */

	uint32 v = 3 * n ^ 2;
	v = v * ((uint32)2 - n * v);
	v = v * ((uint32)2 - n * v);
	v = v * ((uint32)2 - n * v);

	return v;
}

__device__ uint64
modwinv_64(uint64 n) { /* compute 1/n mod 2**64 */

	uint64 v = 3 * n ^ 2;
	v = v * ((uint64)2 - n * v);
	v = v * ((uint64)2 - n * v);
	v = v * ((uint64)2 - n * v);
	v = v * ((uint64)2 - n * v);

	return v;
}

__device__ uint32
mod160t_32(uint32 n) { /* compute 2**160 mod n */

	uint32 t = -n;
	t = ((uint64)t * t) % n;
	t = ((uint64)t * t) % n;
	t = ((uint64)t * -n) % n;

	return t;
}

__device__ uint64
mod192t_64(uint64 n) { /* compute 2**192 mod n */

	/* assumes n < 2**63 */

	uint64 a, b, t = 0;

	a = b = -n % n;

	while (a != 0) {
		if ((a & 1) == 1) {
			t += b;
			if (t >= n)
				t -= n;
		}

		a >>= 1;
		b <<= 1;
		if (b >= n)
			b -= n;
	}

	a = t;
	b = -n % n;
	t = 0;

	while (a != 0) {
		if ((a & 1) == 1) {
			t += b;
			if (t >= n)
				t -= n;
		}

		a >>= 1;
		b <<= 1;
		if (b >= n)
			b -= n;
	}

	return t;
}

__device__ uint32
mod128_32(uint128 num, uint32 v, uint32 t,
		uint32 n, uint32 w) {

	uint32 x, rem, tmp;

	x = num.w[0];
	tmp = x * v;
	rem = __umulhi(tmp, n);
	x = num.w[1];
	tmp = x - rem;
	tmp = tmp * v;
	if (rem > x)
		tmp++;
	rem = __umulhi(tmp, n);
	x = num.w[2];
	tmp = x - rem;
	tmp = tmp * v;
	if (rem > x)
		tmp++;
	rem = __umulhi(tmp, n);
	x = num.w[3];
	tmp = x - rem;
	tmp = tmp * v;
	if (rem > x)
		tmp++;
	rem = __umulhi(tmp, n);
	rem = montmul32(rem, t, n, w);

	if (rem == 0)
		return 0;
	else
		return n - rem;
}

__device__  uint64
mod128_64(uint128 num, uint64 v, uint64 t,
		uint64 n, uint32 w) {

	uint64 x, rem, tmp;

	x = (uint64)num.w[1] << 32 | num.w[0];
	tmp = x * v;
	rem = mul64hi(tmp, n);
	x = (uint64)num.w[3] << 32 | num.w[2];
	tmp = x - rem;
	tmp = tmp * v;
	if (rem > x)
		tmp++;
	rem = mul64hi(tmp, n);
	rem = montmul64(rem, t, n, w);

	if (rem == 0)
		return 0;
	else
		return n - rem;
}

#ifdef __cplusplus
}
#endif

#endif /* defined(__CUDACC__) && !defined(CUDA_INTRINSICS_H) */

