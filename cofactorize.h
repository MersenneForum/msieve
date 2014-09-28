#ifndef _COFACTORIZE_H_
#define _COFACTORIZE_H_

#include <gmp.h>

typedef signed char s8;
typedef unsigned char u8;
typedef signed short s16;
typedef unsigned short u16;
typedef signed int s32;
typedef unsigned int u32;

#ifdef _MSC_VER
typedef signed __int64 s64;
typedef unsigned __int64 u64;
#else
typedef long long s64;
typedef unsigned long long u64;
#endif

u32 tinyqs(mpz_t n, mpz_t factor1, mpz_t factor2);
u32 squfof(mpz_t n);

#endif /* !_COFACTORIZE_H_ */
