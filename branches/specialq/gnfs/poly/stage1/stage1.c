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

#include <stage1.h>

/*------------------------------------------------------------------------*/
void 
stage1_bounds_init(bounds_t *bounds, poly_stage1_t *data)
{
	mpz_init_set(bounds->gmp_high_coeff_begin, 
			data->gmp_high_coeff_begin);
	mpz_init_set(bounds->gmp_high_coeff_end, 
			data->gmp_high_coeff_end);
	bounds->norm_max = data->norm_max;
}

/*------------------------------------------------------------------------*/
void 
stage1_bounds_free(bounds_t *bounds)
{
	mpz_clear(bounds->gmp_high_coeff_begin);
	mpz_clear(bounds->gmp_high_coeff_end);
}

/*------------------------------------------------------------------------*/
void 
stage1_bounds_update(bounds_t *bounds, double N, 
			double high_coeff, uint32 degree)
{
	double skewness_min = 1.0;

	switch (degree) {
	case 4:
		skewness_min = sqrt(pow(N / high_coeff, 1./4.) / 
					bounds->norm_max);
		bounds->coeff_max = bounds->norm_max;
		break;

	case 5:
		skewness_min = pow(pow(N / high_coeff, 1./5.) / 
					bounds->norm_max, 2./3.);
		bounds->coeff_max = bounds->norm_max / sqrt(skewness_min);
		break;

	case 6:
		skewness_min = sqrt(pow(N / high_coeff, 1./6.) / 
					bounds->norm_max);
		bounds->coeff_max = bounds->norm_max / skewness_min;
		break;

	default:
		printf("error: unexpected poly degree %d\n", degree);
		exit(-1);
	}

	bounds->p_size_max = bounds->coeff_max / skewness_min;
}

/*------------------------------------------------------------------------*/
void
poly_search_init(poly_search_t *poly, poly_stage1_t *data)
{
	uint32 i;

	for (i = 0; i < POLY_BATCH_SIZE; i++) {
		curr_poly_t *c = poly->batch + i;
		mpz_init(c->high_coeff);
		mpz_init(c->trans_N);
		mpz_init(c->trans_m0);
		mpz_init(c->mp_sieve_size);
	}
	mpz_init_set(poly->N, data->gmp_N);
	mpz_init(poly->m0);
	mpz_init(poly->p);
	mpz_init(poly->tmp1);
	mpz_init(poly->tmp2);
	mpz_init(poly->tmp3);
	mpz_init(poly->tmp4);
	mpz_init(poly->tmp5);

	poly->degree = data->degree;
	poly->callback = data->callback;
	poly->callback_data = data->callback_data;

#ifdef HAVE_CUDA
	CUDA_TRY(cuCtxCreate(&poly->gpu_context, 
			CU_CTX_BLOCKING_SYNC,
			poly->gpu_info->device_handle))

	switch (poly->degree) {
	case 4:
	case 6:
		CUDA_TRY(cuModuleLoad(&poly->gpu_module48, 
				"stage1_core_deg46_48.ptx"))
		CUDA_TRY(cuModuleLoad(&poly->gpu_module64, 
				"stage1_core_deg46_64.ptx"))
		break;

	case 5:
		CUDA_TRY(cuModuleLoad(&poly->gpu_module48, 
				"stage1_core_deg5_48.ptx"))
		CUDA_TRY(cuModuleLoad(&poly->gpu_module64, 
				"stage1_core_deg5_64.ptx"))
		break;
	}
#endif

}

/*------------------------------------------------------------------------*/
void
poly_search_free(poly_search_t *poly)
{
	uint32 i;

	for (i = 0; i < POLY_BATCH_SIZE; i++) {
		curr_poly_t *c = poly->batch + i;
		mpz_clear(c->high_coeff);
		mpz_clear(c->trans_N);
		mpz_clear(c->trans_m0);
		mpz_clear(c->mp_sieve_size);
	}
	mpz_clear(poly->N);
	mpz_clear(poly->m0);
	mpz_clear(poly->p);
	mpz_clear(poly->tmp1);
	mpz_clear(poly->tmp2);
	mpz_clear(poly->tmp3);
	mpz_clear(poly->tmp4);
	mpz_clear(poly->tmp5);
#ifdef HAVE_CUDA
	CUDA_TRY(cuCtxDestroy(poly->gpu_context)) 
#endif
}

/*------------------------------------------------------------------------*/
void
handle_collision(poly_search_t *poly, uint32 which_poly,
			uint32 p1, uint32 p2, uint32 special_q,
			uint64 special_q_root, uint128 res)
{
	curr_poly_t *c = poly->batch + which_poly;

	/* p1, p2, and special_q should always be pairwise coprime
	 * when we get here, but let's be defensive and check anyway. */
	if (mp_gcd_1(special_q, p1) != 1)
		return;
	if (mp_gcd_1(special_q, p2) != 1)
		return;
	if (mp_gcd_1(p1, p2) != 1)
		return;

	mpz_set_ui(poly->p, (unsigned long)p1);
	mpz_mul_ui(poly->p, poly->p, (unsigned long)p2);
	mpz_mul_ui(poly->p, poly->p, (unsigned long)special_q);

	mpz_gcd(poly->tmp3, poly->p, c->high_coeff);
	if (mpz_cmp_ui(poly->tmp3, 1))
		return;

	uint64_2gmp(special_q_root, poly->tmp1);
	mpz_import(poly->tmp2, 4, -1, sizeof(uint32), 0, 0, &res);
	mpz_set_ui(poly->tmp3, (unsigned long)special_q);

	mpz_mul(poly->tmp3, poly->tmp3, poly->tmp3);
	mpz_addmul(poly->tmp1, poly->tmp2, poly->tmp3);
	mpz_sub(poly->tmp1, poly->tmp1, c->mp_sieve_size);
	mpz_add(poly->m0, c->trans_m0, poly->tmp1);

	/* check */
	mpz_pow_ui(poly->tmp1, poly->m0, (mp_limb_t)poly->degree);
	mpz_mul(poly->tmp2, poly->p, poly->p);
	mpz_sub(poly->tmp1, c->trans_N, poly->tmp1);
	mpz_tdiv_r(poly->tmp3, poly->tmp1, poly->tmp2);
	if (mpz_cmp_ui(poly->tmp3, (mp_limb_t)0)) {
		gmp_printf("poly %u %u %u %u %Zd\n",
				which_poly, special_q, p1, p2, poly->m0);
		printf("crap\n");
		return;
	}

	mpz_mul_ui(poly->tmp1, c->high_coeff, (mp_limb_t)poly->degree);
	mpz_tdiv_qr(poly->m0, poly->tmp2, poly->m0, poly->tmp1);
	mpz_invert(poly->tmp3, poly->tmp1, poly->p);

	mpz_sub(poly->tmp4, poly->tmp3, poly->p);
	if (mpz_cmpabs(poly->tmp3, poly->tmp4) < 0)
		mpz_set(poly->tmp4, poly->tmp3);

	mpz_sub(poly->tmp5, poly->tmp2, poly->tmp1);
	if (mpz_cmpabs(poly->tmp2, poly->tmp5) > 0)
		mpz_add_ui(poly->m0, poly->m0, (mp_limb_t)1);
	else
		mpz_set(poly->tmp5, poly->tmp2);

	mpz_addmul(poly->m0, poly->tmp4, poly->tmp5);

	gmp_printf("poly %2u %u %u %u %Zd\n",
			which_poly, special_q, p1, p2, poly->m0);

	poly->callback(c->high_coeff, poly->p, poly->m0,
			c->coeff_max, poly->callback_data);
}

/*------------------------------------------------------------------------*/
static void
search_coeffs_core(msieve_obj *obj, poly_search_t *poly, 
			uint32 deadline)
{
	uint32 i, j;
	uint32 degree = poly->degree;
	uint32 num_poly = poly->num_poly;
	uint32 mult = 0;

	switch (degree) {
	case 4: mult = 4 * 4 * 4 * 4; break;
	case 5: mult = 5 * 5 * 5 * 5 * 5; break;
	case 6: mult = 6 * 6 * 6 * 6 * 6 * 6; break;
	}

	for (i = 0; i < num_poly; i++) {
		curr_poly_t *c = poly->batch + i;

		mpz_mul_ui(c->trans_N, poly->N, (mp_limb_t)mult);
		for (j = 0; j < degree - 1; j++)
			mpz_mul(c->trans_N, c->trans_N, c->high_coeff);

		mpz_root(c->trans_m0, c->trans_N, (mp_limb_t)degree);

		mpz_tdiv_q(poly->m0, poly->N, c->high_coeff);
		mpz_root(poly->m0, poly->m0, (mp_limb_t)degree);
		c->sieve_size = c->coeff_max / mpz_get_d(poly->m0) * 
				c->p_size_max * c->p_size_max / degree;
		mpz_set_d(c->mp_sieve_size, c->sieve_size);
	}

	sieve_lattice(obj, poly, num_poly * deadline);
}

/*------------------------------------------------------------------------*/
static void
search_coeffs(msieve_obj *obj, poly_search_t *poly, 
		bounds_t *bounds, uint32 deadline)
{
	uint32 i, j, p;
	mpz_t curr_high_coeff, tmp;
	double dn = mpz_get_d(poly->N);
	uint32 digits = mpz_sizeinbase(poly->N, 10);
	double start_time = get_cpu_time();
	uint32 deadline_per_coeff = 800;
#ifdef HAVE_CUDA
	uint32 batch_size = (poly->degree == 5) ? POLY_BATCH_SIZE : 1;
#else
	uint32 batch_size = 1;
#endif

	if (digits <= 100)
		deadline_per_coeff = 5;
	else if (digits <= 105)
		deadline_per_coeff = 20;
	else if (digits <= 110)
		deadline_per_coeff = 30;
	else if (digits <= 120)
		deadline_per_coeff = 50;
	else if (digits <= 130)
		deadline_per_coeff = 100;
	else if (digits <= 140)
		deadline_per_coeff = 200;
	else if (digits <= 150)
		deadline_per_coeff = 400;
	printf("deadline: %u seconds per coefficient\n", deadline_per_coeff);

	mpz_init(curr_high_coeff);
	mpz_init(tmp);
	mpz_sub_ui(curr_high_coeff, bounds->gmp_high_coeff_begin, (mp_limb_t)1);
	mpz_fdiv_q_ui(curr_high_coeff, curr_high_coeff, 
			(mp_limb_t)HIGH_COEFF_MULTIPLIER);
	mpz_mul_ui(curr_high_coeff, curr_high_coeff, 
			(mp_limb_t)HIGH_COEFF_MULTIPLIER);

	poly->num_poly = 0;

	while (1) {
		curr_poly_t *c = poly->batch + poly->num_poly;

		mpz_add_ui(curr_high_coeff, curr_high_coeff,
				(mp_limb_t)HIGH_COEFF_MULTIPLIER);

		if (mpz_cmp(curr_high_coeff, 
					bounds->gmp_high_coeff_end) > 0) {

			if (poly->num_poly > 0) {
				search_coeffs_core(obj, poly, 
					   deadline_per_coeff);
			}
			break;
		}

		mpz_divexact_ui(tmp, curr_high_coeff, 
					(mp_limb_t)HIGH_COEFF_MULTIPLIER);
		for (i = p = 0; i < PRECOMPUTED_NUM_PRIMES; i++) {
			p += prime_delta[i];

			if (p > HIGH_COEFF_PRIME_LIMIT)
				break;

			for (j = 0; j < HIGH_COEFF_POWER_LIMIT; j++) {
				if (mpz_divisible_ui_p(tmp, (mp_limb_t)p))
					mpz_divexact_ui(tmp, tmp, (mp_limb_t)p);
				else
					break;
			}
		}
		if (mpz_cmp_ui(tmp, (mp_limb_t)1))
			continue;

		stage1_bounds_update(bounds, dn, 
					mpz_get_d(curr_high_coeff),
					poly->degree);

		mpz_set(c->high_coeff, curr_high_coeff);
		c->p_size_max = bounds->p_size_max;
		c->coeff_max = bounds->coeff_max;

		if (++poly->num_poly == batch_size) {
			search_coeffs_core(obj, poly, 
					deadline_per_coeff);

			if (obj->flags & MSIEVE_FLAG_STOP_SIEVING)
				break;

			if (deadline) {
				double curr_time = get_cpu_time();
				double elapsed = curr_time - start_time;

				if (elapsed > deadline)
					break;
			}
			poly->num_poly = 0;
		}
	}

	mpz_clear(curr_high_coeff);
}

/*------------------------------------------------------------------------*/
void
poly_stage1_init(poly_stage1_t *data,
		 stage1_callback_t callback, void *callback_data)
{
	memset(data, 0, sizeof(poly_stage1_t));
	mpz_init_set_ui(data->gmp_N, (mp_limb_t)0);
	mpz_init_set_ui(data->gmp_high_coeff_begin, (mp_limb_t)0);
	mpz_init_set_ui(data->gmp_high_coeff_end, (mp_limb_t)0);
	data->callback = callback;
	data->callback_data = callback_data;
}

/*------------------------------------------------------------------------*/
void
poly_stage1_free(poly_stage1_t *data)
{
	mpz_clear(data->gmp_N);
	mpz_clear(data->gmp_high_coeff_begin);
	mpz_clear(data->gmp_high_coeff_end);
}

/*------------------------------------------------------------------------*/
uint32
poly_stage1_run(msieve_obj *obj, poly_stage1_t *data)
{
	bounds_t bounds;
	poly_search_t poly;
#ifdef HAVE_CUDA
	gpu_config_t gpu_config;

	gpu_init(&gpu_config);
	if (gpu_config.num_gpu == 0) {
		printf("error: no CUDA-enabled GPUs found\n");
		exit(-1);
	}
	if (obj->which_gpu >= (uint32)gpu_config.num_gpu) {
		printf("error: GPU %u does not exist "
			"or is not CUDA-enabled\n", obj->which_gpu);
		exit(-1);
	}
	logprintf(obj, "using GPU %u (%s)\n", obj->which_gpu,
			gpu_config.info[obj->which_gpu].name);

	poly.gpu_info = gpu_config.info + obj->which_gpu; 
#endif

	stage1_bounds_init(&bounds, data);
	poly_search_init(&poly, data);

	search_coeffs(obj, &poly, &bounds, data->deadline);

	poly_search_free(&poly);
	stage1_bounds_free(&bounds);
	return 1;
}
