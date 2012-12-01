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

#ifndef _STAGE1_CORE_H_
#define _STAGE1_CORE_H_

#ifdef __CUDACC__
#include "cuda_intrinsics.h"
#endif

#ifdef __cplusplus
extern "C" {
#endif

#define FOUND_ARRAY_SIZE 1000

typedef struct {
	uint64 p;
	uint128 pp;
	uint128 root;
} specialq_t;

typedef struct {
	uint32 p1;
	uint32 p2;
	int64 offset;
	specialq_t q;
} found_t;

#ifdef __cplusplus
}

#endif

#endif /* !_STAGE1_CORE_H_ */
