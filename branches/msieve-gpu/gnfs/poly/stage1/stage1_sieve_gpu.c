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

#include "stage1.h"
#include "stage1_core.h"

#define HOST_BATCH_SIZE 30000

/*------------------------------------------------------------------------*/
typedef struct {
	uint32 num_p;
	uint32 last_p;

	uint32 p[HOST_BATCH_SIZE];
	uint32 lattice_size[HOST_BATCH_SIZE];
	uint64 roots[POLY_BATCH_SIZE][HOST_BATCH_SIZE];
} p_soa_var_t;

static void
p_soa_var_reset(p_soa_var_t *soa)
{
	soa->num_p = 0;
	soa->last_p = 0;
}

static void 
store_p_soa(uint64 p, uint32 num_roots, uint32 which_poly, 
		mpz_t *roots, void *extra)
{
	lattice_fb_t *L = (lattice_fb_t *)extra;
	p_soa_var_t *soa;
	uint32 num;

	if (num_roots != 1) {
		printf("error: num_roots > 1\n");
		exit(-1);
	}

	if (L->fill_p)
		soa = (p_soa_var_t *)L->p_array;
	else
		soa = (p_soa_var_t *)L->q_array;

	num = soa->num_p;
	if (p != soa->last_p) {
		soa->p[num] = (uint32)p;
		soa->lattice_size[num] = (uint32)(L->poly->batch[
				which_poly].sieve_size / ((double)p * p));
		soa->num_p++;
		soa->last_p = (uint32)p;
		soa->roots[which_poly][num] = gmp2uint64(roots[0]);
	}
	else {
		soa->roots[which_poly][num - 1] = gmp2uint64(roots[0]);
	}
}

/*------------------------------------------------------------------------*/
uint32
sieve_lattice_batch(msieve_obj *obj, lattice_fb_t *L,
			uint32 threads_per_block,
			gpu_info_t *gpu_info, CUfunction gpu_kernel)
{
	uint32 i;
	p_soa_var_t * p_array = (p_soa_var_t *)L->p_array;
	p_soa_var_t * q_array = (p_soa_var_t *)L->q_array;
	uint32 num_blocks;
	uint32 num_p_offset;
	uint32 num_q_offset;
	void *gpu_ptr;
	uint32 num_poly = L->poly->num_poly;

	p_soa_t *p_marshall = L->p_marshall;
	p_soa_t *q_marshall = L->q_marshall;
	uint32 num_q_done = 0;

	i = 0;
	gpu_ptr = (void *)(size_t)L->gpu_p_array;
	CUDA_ALIGN_PARAM(i, __alignof(gpu_ptr));
	CUDA_TRY(cuParamSetv(gpu_kernel, (int)i, 
			&gpu_ptr, sizeof(gpu_ptr)))
	i += sizeof(gpu_ptr);

	CUDA_ALIGN_PARAM(i, __alignof(uint32));
	num_p_offset = i;
	i += sizeof(uint32);

	gpu_ptr = (void *)(size_t)L->gpu_q_array;
	CUDA_ALIGN_PARAM(i, __alignof(gpu_ptr));
	CUDA_TRY(cuParamSetv(gpu_kernel, (int)i, 
			&gpu_ptr, sizeof(gpu_ptr)))
	i += sizeof(gpu_ptr);

	CUDA_ALIGN_PARAM(i, __alignof(uint32));
	num_q_offset = i;
	i += sizeof(uint32);

	CUDA_ALIGN_PARAM(i, __alignof(uint32));
	CUDA_TRY(cuParamSeti(gpu_kernel, i, num_poly))
	i += sizeof(uint32);

	gpu_ptr = (void *)(size_t)L->gpu_found_array;
	CUDA_ALIGN_PARAM(i, __alignof(gpu_ptr));
	CUDA_TRY(cuParamSetv(gpu_kernel, (int)i, 
			&gpu_ptr, sizeof(gpu_ptr)))
	i += sizeof(gpu_ptr);

	CUDA_TRY(cuParamSetSize(gpu_kernel, i))

	while (num_q_done < q_array->num_p) {

		uint32 q_left;
		uint32 num_p_done = 0;
		time_t curr_time;
		double elapsed;

		uint32 curr_num_q = MIN(L->found_array_size,
					q_array->num_p - num_q_done);

		curr_num_q = MIN(curr_num_q, P_SOA_BATCH_SIZE);
		q_left = q_array->num_p - (num_q_done + curr_num_q);
		if (q_left > 0 && q_left < L->found_array_size / 4)
			curr_num_q /= 2;

		memcpy(q_marshall->p, 
			q_array->p + num_q_done,
			curr_num_q * sizeof(uint32));

		for (i = 0; i < num_poly; i++) {
			memcpy(q_marshall->roots[i],
				q_array->roots[i] + num_q_done,
				curr_num_q * sizeof(uint64));
		}

		CUDA_TRY(cuMemcpyHtoD(L->gpu_q_array, q_marshall,
				P_SOA_BATCH_SIZE * (2 * sizeof(uint32) +
					num_poly * sizeof(uint64))))
		CUDA_TRY(cuParamSeti(gpu_kernel, num_q_offset, curr_num_q))

		while (num_p_done < p_array->num_p) {

			uint32 p_left;
			uint32 curr_num_p = MIN(L->found_array_size,
						p_array->num_p - num_p_done);

			curr_num_p = MIN(curr_num_p, P_SOA_BATCH_SIZE);
			p_left = p_array->num_p - (num_p_done + curr_num_p);
			if (p_left > 0 && p_left < L->found_array_size / 4)
				curr_num_p /= 2;

			memcpy(p_marshall->p, 
				p_array->p + num_p_done,
				curr_num_p * sizeof(uint32));
			memcpy(p_marshall->lattice_size, 
				p_array->lattice_size + num_p_done,
				curr_num_p * sizeof(uint32));

			for (i = 0; i < num_poly; i++) {
				memcpy(p_marshall->roots[i],
					p_array->roots[i] + num_p_done,
					curr_num_p * sizeof(uint64));
			}

			CUDA_TRY(cuMemcpyHtoD(L->gpu_p_array, p_marshall,
				P_SOA_BATCH_SIZE * (2 * sizeof(uint32) +
					num_poly * sizeof(uint64))))

			CUDA_TRY(cuParamSeti(gpu_kernel, num_p_offset, 
						curr_num_p))
#if 0
			printf("qnum %u pnum %u\n", curr_num_q, curr_num_p);
#endif

			num_blocks = gpu_info->num_compute_units;
			if (curr_num_q < L->found_array_size) {
				num_blocks = (curr_num_q + 
					threads_per_block - 1) /
					threads_per_block;
			}

			CUDA_TRY(cuLaunchGrid(gpu_kernel, 
						num_blocks, 1))

			CUDA_TRY(cuMemcpyDtoH(L->found_array, 
						L->gpu_found_array, 
						num_blocks * 
						threads_per_block *
							sizeof(found_t)))

			for (i = 0; i < threads_per_block *
					num_blocks; i++) {
				found_t *f = L->found_array + i;

				if (f->p > 0) {
					handle_collision(L->poly, 
						f->which_poly,
						f->p, f->proot, 
						f->offset, f->q);
				}
			}

			num_p_done += curr_num_p;
		}

		if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			return 1;

		curr_time = time(NULL);
		elapsed = curr_time - L->start_time;
		if (elapsed > L->deadline)
			return 1;

		num_q_done += curr_num_q;
	}

	return 0;
}

/*------------------------------------------------------------------------*/
uint32
sieve_lattice_gpu(msieve_obj *obj, lattice_fb_t *L, 
		sieve_fb_t *sieve_small, sieve_fb_t *sieve_large, 
		uint32 small_p_min, uint32 small_p_max, 
		uint32 large_p_min, uint32 large_p_max,
		gpu_info_t *gpu_info, CUmodule gpu_module)
{
	uint32 i;
	uint32 min_small, min_large;
	uint32 quit = 0;
	p_soa_var_t * p_array;
	p_soa_var_t * q_array;
	clock_t clock_start;

	CUfunction gpu_kernel;
	uint32 threads_per_block;

	L->p_marshall = (p_soa_t *)xmalloc(sizeof(p_soa_t));
	L->q_marshall = (p_soa_t *)xmalloc(sizeof(p_soa_t));
	p_array = L->p_array = (p_soa_var_t *)xmalloc(
					sizeof(p_soa_var_t));
	q_array = L->q_array = (p_soa_var_t *)xmalloc(
					sizeof(p_soa_var_t));

	CUDA_TRY(cuModuleGetFunction(&gpu_kernel, gpu_module, 
				"sieve_kernel"))
	CUDA_TRY(cuMemAlloc(&L->gpu_p_array, sizeof(p_soa_t)))
	CUDA_TRY(cuMemAlloc(&L->gpu_q_array, sizeof(p_soa_t)))

	CUDA_TRY(cuFuncGetAttribute((int *)&threads_per_block, 
			CU_FUNC_ATTRIBUTE_MAX_THREADS_PER_BLOCK,
			gpu_kernel))

	CUDA_TRY(cuFuncSetBlockShape(gpu_kernel, 
				threads_per_block, 1, 1))

	L->found_array_size = threads_per_block *
				gpu_info->num_compute_units;
	L->found_array = (found_t *)xmalloc(L->found_array_size *
					sizeof(found_t));
	CUDA_TRY(cuMemAlloc(&L->gpu_found_array, 
			L->found_array_size * sizeof(found_t)))

	printf("------- %u-%u %u-%u\n",
			small_p_min, small_p_max,
			large_p_min, large_p_max);

	min_large = large_p_min;
	sieve_fb_reset(sieve_small, (uint64)large_p_min, 
			(uint64)large_p_max, 1, 1);

	clock_start = clock();

	while (min_large < large_p_max) {

		L->fill_p = 0;
		p_soa_var_reset(q_array);
		for (i = 0; i < HOST_BATCH_SIZE && 
				min_large < large_p_max; i++) {
			min_large = sieve_fb_next(sieve_small, L->poly,
						store_p_soa, L);
		}
		if (q_array->num_p == 0)
			goto finished;

		min_small = small_p_min;
		sieve_fb_reset(sieve_large, 
				(uint64)small_p_min, (uint64)small_p_max,
				1, 1);

		while (min_small <= small_p_max) {

			L->fill_p = 1;
			p_soa_var_reset(p_array);
			for (i = 0; i < HOST_BATCH_SIZE && 
					min_small < small_p_max; i++) {
				min_small = sieve_fb_next(sieve_large, L->poly,
							store_p_soa, L);
			}
			if (p_array->num_p == 0)
				goto finished;

			if (sieve_lattice_batch(obj, L, threads_per_block,
						gpu_info, gpu_kernel)) {
				quit = 1;
				goto finished;
			}
		}
	}

	printf("%lf\n", (double)(clock() - clock_start) / CLOCKS_PER_SEC);

finished:
	CUDA_TRY(cuMemFree(L->gpu_p_array))
	CUDA_TRY(cuMemFree(L->gpu_q_array))
	CUDA_TRY(cuMemFree(L->gpu_found_array))
	free(p_array);
	free(q_array);
	free(L->p_marshall);
	free(L->q_marshall);
	free(L->found_array);
	return quit;
}


