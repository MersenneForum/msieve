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

#include "stage2.h"

typedef struct {
	uint64 x;
	uint64 y;
	uint64 z;
} lattice_t;

#define MAX_FACTORS 8

typedef struct {
	uint32 start;
	uint32 stride_y;
	uint32 stride_z;
} xyprog_t;

typedef struct {
	uint32 p;
	uint32 latsize_mod_p;
	uint32 table_size;
	uint32 num_roots;
	uint32 max_sieve_val;
	xyprog_t *roots;
	uint16 *invtable;
	uint8 *sieve;
} xydata_t;

/*-------------------------------------------------------------------------*/
#define APPLY_ROTATION(dist) {						\
	double x_dist = floor((dist) * direction[0]);			\
	double y_dist = floor((dist) * direction[1]);			\
	double z_dist = floor((dist) * direction[2]);			\
	bpoly.coeff[0] = apoly->coeff[0] - dbl_d * x_dist;		\
	bpoly.coeff[1] = apoly->coeff[1] + dbl_p * x_dist - dbl_d * y_dist; \
	bpoly.coeff[2] = apoly->coeff[2] + dbl_p * y_dist - dbl_d * z_dist; \
	bpoly.coeff[3] = apoly->coeff[3] + dbl_p * z_dist;		\
}

static void
compute_line_size(double max_norm, dpoly_t *apoly,
		  double dbl_p, double dbl_d, double direction[3],
		  double last_line_min_in, double last_line_max_in,
		  double *line_min, double *line_max)
{
	double v0;
	double new_xlate, new_skewness;
	double d0, d1, offset;
	double last_line_min = last_line_min_in;
	double last_line_max = last_line_max_in;
	dpoly_t bpoly = *apoly;

	APPLY_ROTATION(last_line_min);
	v0 = optimize_basic(&bpoly, &new_skewness, &new_xlate);
	offset = 1e-6 * fabs(last_line_min);
	offset = MAX(offset, 10000.0);

	if (v0 > max_norm) {
		d0 = last_line_min;
		d1 = last_line_min + offset;
		while (1) {
			double new_d1 = last_line_min + offset;
			APPLY_ROTATION(new_d1);
			v0 = optimize_basic(&bpoly, &new_skewness, &new_xlate);
			if (v0 <= max_norm || new_d1 >= last_line_max) {
				d1 = new_d1;
				break;
			}
			d0 = new_d1;
			offset *= 2;
		}
	}
	else {
		d0 = last_line_min - offset;
		d1 = last_line_min;
		while (1) {
			double new_d0 = last_line_min - offset;
			APPLY_ROTATION(new_d0);
			v0 = optimize_basic(&bpoly, &new_skewness, &new_xlate);
			if (v0 > max_norm) {
				d0 = new_d0;
				break;
			}
			d1 = new_d0;
			offset *= 2;
		}
	}

	while (d1 - d0 > 500 && d1 - d0 > 1e-6 * fabs(d0)) {
		double dd = (d0 + d1) / 2;
		APPLY_ROTATION(dd);
		v0 = optimize_basic(&bpoly, &new_skewness, &new_xlate);
		if (v0 > max_norm)
			d0 = dd;
		else
			d1 = dd;
	}
	*line_min = d0;

	APPLY_ROTATION(last_line_max);
	v0 = optimize_basic(&bpoly, &new_skewness, &new_xlate);
	offset = 1e-6 * fabs(last_line_max);
	offset = MAX(offset, 10000.0);

	if (v0 > max_norm) {
		d0 = last_line_max - offset;
		d1 = last_line_max;
		while (1) {
			double new_d0 = last_line_max - offset;
			APPLY_ROTATION(new_d0);
			v0 = optimize_basic(&bpoly, &new_skewness, &new_xlate);
			if (v0 <= max_norm || new_d0 <= last_line_min) {
				d0 = new_d0;
				break;
			}
			d1 = new_d0;
			offset *= 2;
		}
	}
	else {
		d0 = last_line_max;
		d1 = last_line_max + offset;
		while (1) {
			double new_d1 = last_line_max + offset;
			APPLY_ROTATION(new_d1);
			v0 = optimize_basic(&bpoly, &new_skewness, &new_xlate);
			if (v0 > max_norm) {
				d1 = new_d1;
				break;
			}
			d0 = new_d1;
			offset *= 2;
		}
	}

	while (d1 - d0 > 500 && d1 - d0 > 1e-6 * fabs(d0)) {
		double dd = (d0 + d1) / 2;
		APPLY_ROTATION(dd);
		v0 = optimize_basic(&bpoly, &new_skewness, &new_xlate);
		if (v0 > max_norm)
			d1 = dd;
		else
			d0 = dd;
	}
	*line_max = d1;
}

/*-------------------------------------------------------------------------*/
static uint64
find_lattice_size_xyz(double line_length)
{
	return (uint64)2*2*2*2*2*3*3*3*5*5*7;
}

/*-------------------------------------------------------------------------*/
static uint32
find_lattice_primes_xyz(sieve_prime_t *primes, uint32 num_primes,
			uint64 lattice_size, sieve_prime_t *lattice_primes)
{
	uint32 i, j;
	uint32 num_lattice_primes = 0;

	for (i = 0; i < num_primes; i++) {

		sieve_prime_t *curr_prime = primes + i;
		uint32 num_powers = curr_prime->num_powers;
		uint32 num_lattice_powers = 0;

		if (lattice_size % curr_prime->prime)
			break;

		for (j = 0; j < num_powers; j++, num_lattice_powers++) {
			sieve_power_t *sp = curr_prime->powers + j;
			if (lattice_size % sp->power)
				break;
		}

		lattice_primes[num_lattice_primes] = *curr_prime;
		lattice_primes[num_lattice_primes].num_powers = 
						num_lattice_powers;
		num_lattice_primes++;
	}

	return num_lattice_primes;
}

/*-------------------------------------------------------------------------*/
static uint32
find_lattice_primes_xy(sieve_prime_t *primes, uint32 num_primes,
			uint64 lattice_size_xyz, sieve_prime_t *lattice_primes,
			double line_length)
{
	uint32 i;
	uint32 num_lattice_primes = 0;
	double target_size = line_length / (double)lattice_size_xyz / 1000;
	double curr_size = 1.0;

	for (i = 0; i < num_primes; i++) {

		sieve_prime_t *curr_prime = primes + i;

		if (num_lattice_primes == MAX_FACTORS)
			break;

		if (lattice_size_xyz % curr_prime->prime == 0)
			continue;

		if (curr_size * curr_prime->prime >= target_size)
			break;

		lattice_primes[num_lattice_primes] = *curr_prime;
		lattice_primes[num_lattice_primes].num_powers = 1;
		curr_size *= curr_prime->prime;
		num_lattice_primes++;
	}

	return num_lattice_primes;
}

/*-------------------------------------------------------------------------*/
#define MAX_DIM 64

#define MAX_ROOTS MAX_DIM

typedef struct {
	uint16 score;
	uint8 power;
	uint8 num_roots;
	uint8 roots[MAX_ROOTS][3];
} hit_t;

/*-------------------------------------------------------------------------*/
static void 
do_sieving_xyz(sieve_root_t *r, uint16 *sieve,
		uint32 contrib, uint32 dim)
{
	uint32 i, j;
	uint32 start = r->start;
	uint32 step = r->step;
	uint32 resclass = r->resclass;
	uint32 resclass2 = mp_modmul_1(resclass, resclass, step);

	for (i = 0; i < dim; i++) {

		uint32 ri = start;

		for (j = 0; j < dim; j++) {

			uint32 rj = ri;

			while (rj < dim) {
				sieve[rj] += contrib;
				rj += step;
			}
			sieve += dim;
			ri = mp_modsub_1(ri, resclass, step);
		}

		start = mp_modsub_1(start, resclass2, step);
	}
}

/*-------------------------------------------------------------------------*/
static void 
find_hits_xyz(sieve_prime_t *lattice_primes,
	uint32 num_lattice_primes, hit_t *hitlist)
{
	uint32 i, j, k;
	uint16 *sieve;

	sieve = (uint16 *)xmalloc(MAX_DIM * MAX_DIM *
				MAX_DIM * sizeof(uint16));

	for (i = 0; i < num_lattice_primes; i++) {

		sieve_prime_t *curr_prime = lattice_primes + i;
		uint32 num_powers = curr_prime->num_powers;
		uint32 dim = curr_prime->powers[num_powers-1].power;
		uint32 size = dim * dim * dim;
		uint32 max_score;
		hit_t *hits = hitlist + i;

		memset(sieve, 0, size * sizeof(uint16));

		for (j = 0; j < num_powers; j++) {

			sieve_power_t *sp = curr_prime->powers + j;
			uint32 num_roots = sp->num_roots;
			uint32 contrib = sp->sieve_contrib;

			for (k = 0; k < num_roots; k++) {
				do_sieving_xyz(sp->roots + k, sieve, 
						contrib, dim);
			}
		}

		for (j = max_score = 0; j < size; j++) {
			uint32 curr_score = sieve[j];
			max_score = MAX(max_score, curr_score);
		}

		for (j = k = 0; j < size; j++) {
			uint32 curr_score = sieve[j];

			if (curr_score == max_score) {
				if (k == MAX_ROOTS)
					break;

				hits->roots[k][0] = j % dim;
				hits->roots[k][1] = (j / dim) % dim;
				hits->roots[k][2] = (j / dim) / dim;
				k++;
			}
		}
		hits->score = max_score;
		hits->power = dim;
		hits->num_roots = k;
	}

	free(sieve);
}

/*-------------------------------------------------------------------------*/
#define CRT_ACCUM(idx)						\
	crt_accum[idx][0] = crt_accum[(idx)+1][0] + 		\
				crt_prod[idx] * 		\
				hitlist[idx].roots[i##idx][0];	\
	crt_accum[idx][1] = crt_accum[(idx)+1][1] + 		\
				crt_prod[idx] * 		\
				hitlist[idx].roots[i##idx][1];	\
	crt_accum[idx][2] = crt_accum[(idx)+1][2] + 		\
				crt_prod[idx] * 		\
				hitlist[idx].roots[i##idx][2];	

static uint32 
compute_lattices_xyz(hit_t *hitlist, uint32 num_lattice_primes,
		   lattice_t **lattices_out, uint64 lattice_size_xyz)
{
	uint32 i;
	int32 i0, i1, i2, i3, i4, i5, i6, i7;
	uint64 crt_prod[MAX_FACTORS];
	uint64 crt_accum[MAX_FACTORS + 1][3];
	uint32 num_lattices;
	lattice_t *lattices;

	for (i = 0, num_lattices = 1; i < num_lattice_primes; i++) {
		hit_t *hits = hitlist + i;
		uint32 p = hits->power;

		crt_prod[i] = lattice_size_xyz / p;
		crt_prod[i] *= mp_modinv_1((uint32)(crt_prod[i] % p), p);
		num_lattices *= hits->num_roots;
	}

	lattices = (lattice_t *)xmalloc(num_lattices * sizeof(lattice_t));
	i0 = i1 = i2 = i3 = i4 = i5 = i6 = i7 = i = 0;
	crt_accum[num_lattice_primes][0] = 0;
	crt_accum[num_lattice_primes][1] = 0;
	crt_accum[num_lattice_primes][2] = 0;

	switch (num_lattice_primes) {
	case 8:
		for (i7 = hitlist[7].num_roots - 1; i7 >= 0; i7--) {
			CRT_ACCUM(7)
	case 7:
		for (i6 = hitlist[6].num_roots - 1; i6 >= 0; i6--) {
			CRT_ACCUM(6)
	case 6:
		for (i5 = hitlist[5].num_roots - 1; i5 >= 0; i5--) {
			CRT_ACCUM(5)
	case 5:
		for (i4 = hitlist[4].num_roots - 1; i4 >= 0; i4--) {
			CRT_ACCUM(4)
	case 4:
		for (i3 = hitlist[3].num_roots - 1; i3 >= 0; i3--) {
			CRT_ACCUM(3)
	case 3:
		for (i2 = hitlist[2].num_roots - 1; i2 >= 0; i2--) {
			CRT_ACCUM(2)
	case 2:
		for (i1 = hitlist[1].num_roots - 1; i1 >= 0; i1--) {
			CRT_ACCUM(1)
	default:
		for (i0 = hitlist[0].num_roots - 1; i0 >= 0; i0--) {
			CRT_ACCUM(0)
			lattices[i].x = crt_accum[0][0] % lattice_size_xyz;
			lattices[i].y = crt_accum[0][1] % lattice_size_xyz;
			lattices[i].z = crt_accum[0][2] % lattice_size_xyz;
			i++;
		}}}}}}}}
	}

	*lattices_out = lattices;
	return num_lattices;
}

/*-------------------------------------------------------------------------*/
static void
xydata_alloc(sieve_prime_t *lattice_primes, 
		uint32 num_lattice_primes, 
		uint64 lattice_size_xyz,
		xydata_t *xydata)
{
	uint32 i;

	for (i = 0; i < num_lattice_primes; i++) {

		sieve_prime_t *curr_prime = lattice_primes + i;
		xydata_t *curr_xydata = xydata + i;
		uint32 p = curr_prime->prime;
		uint32 table_size = p * p;
		uint32 num_roots = curr_prime->powers[0].num_roots;

		curr_xydata->p = p;
		curr_xydata->latsize_mod_p = lattice_size_xyz % p;
		curr_xydata->table_size = table_size;
		curr_xydata->num_roots = num_roots;
		curr_xydata->max_sieve_val = MIN(num_roots, 6);
		curr_xydata->roots = (xyprog_t *)xmalloc(num_roots * 
							sizeof(xyprog_t));
		curr_xydata->invtable = (uint16 *)xmalloc(p * sizeof(uint16));
		curr_xydata->sieve = (uint8 *)xmalloc(table_size * 
							sizeof(uint8));
	}
}

/*-------------------------------------------------------------------------*/
static void
xydata_free(xydata_t *xydata, uint32 num_lattice_primes)
{
	uint32 i;

	for (i = 0; i < num_lattice_primes; i++) {

		xydata_t *curr_xydata = xydata + i;

		free(curr_xydata->roots);
		free(curr_xydata->invtable);
		free(curr_xydata->sieve);
	}
}

/*-------------------------------------------------------------------------*/
static void 
xydata_init(sieve_prime_t *lattice_primes, xydata_t *xydata,
		uint32 num_lattice_primes, 
		lattice_t *lattice_xyz, int64 z_base)
{
	uint32 i, j;

	for (i = 0; i < num_lattice_primes; i++) {

		sieve_prime_t *curr_prime = lattice_primes + i;
		sieve_power_t *sp = curr_prime->powers + 0;
		xydata_t *curr_xydata = xydata + i;

		uint32 p = curr_xydata->p;
		uint32 num_roots = curr_xydata->num_roots;
		uint16 *invtable = curr_xydata->invtable;
		xyprog_t *roots = curr_xydata->roots;

		uint32 latsize_mod_p = curr_xydata->latsize_mod_p;
		uint32 x_mod_p = lattice_xyz->x % p;
		uint32 y_mod_p = lattice_xyz->y % p;
		int64 z_start = z_base + lattice_xyz->z;
		int32 z_start_mod = z_start % p;
		uint32 z_mod_p = (z_start_mod < 0) ? 
				(z_start_mod + (int32)p) : z_start_mod;

		for (j = 0; j < num_roots; j++) {
			sieve_root_t *r = sp->roots + j;
			xyprog_t *curr_xyprog = roots + j;
			uint32 start = r->start;
			uint32 resclass = r->resclass;
			uint32 resclass2 = mp_modmul_1(resclass, resclass, p);

			curr_xyprog->stride_y = mp_modmul_1(resclass, 
							latsize_mod_p, p);
			curr_xyprog->stride_z = mp_modmul_1(resclass2, 
							latsize_mod_p, p);

			start = mp_modsub_1(start, 
					mp_modmul_1(resclass, y_mod_p, p), 
					p);
			curr_xyprog->start = mp_modsub_1(start, 
					mp_modmul_1(resclass2, z_mod_p, p), 
					p);
		}

		for (j = 0; j < p; j++) {
			invtable[x_mod_p] = j;
			x_mod_p = mp_modadd_1(x_mod_p, latsize_mod_p, p);
		}
	}
}

/*-------------------------------------------------------------------------*/
uint32 g_num_max_score;

static void 
find_hits_xy(xydata_t *xydata, uint32 num_lattice_primes, 
		lattice_t *lattice_xyz, int64 z_base, 
		uint32 z_blocks, uint64 lattice_size_xyz)
{
	uint32 i, j, k, m;

	for (i = 0; i < z_blocks; i++) {

		for(j = 0; j < num_lattice_primes; j++) {

			xydata_t *curr_xydata = xydata + j;
			uint32 p = curr_xydata->p;
			uint32 table_size = curr_xydata->table_size;
			uint32 num_roots = curr_xydata->num_roots;
			uint32 max_sieve_val = curr_xydata->max_sieve_val;
			xyprog_t *roots = curr_xydata->roots;
			uint8 *sieve = curr_xydata->sieve;
			uint16 *invtable = curr_xydata->invtable;

			memset(sieve, 0, table_size * sizeof(uint8));

			for (k = 0; k < num_roots; k++) {

				uint8 *row = sieve;
				xyprog_t *curr_prog = roots + k;
				uint32 start = curr_prog->start;
				uint32 stride_y = curr_prog->stride_y;

				for (m = 0; m < p; m++) {
					row[invtable[start]]++;
					row += p;
					start = mp_modsub_1(start, stride_y, p);
				}
			}

			for (k = 0; k < table_size; k++) {
				if (sieve[k] == max_sieve_val)
					break;
			}
			if (k == table_size)
				break;
		}

		if (j == num_lattice_primes) {
			/* yay! */
			g_num_max_score++;
		}

		for(j = 0; j < num_lattice_primes; j++) {

			xydata_t *curr_xydata = xydata + j;
			uint32 p = curr_xydata->p;
			uint32 num_roots = curr_xydata->num_roots;
			xyprog_t *roots = curr_xydata->roots;

			for (k = 0; k < num_roots; k++) {

				xyprog_t *curr_prog = roots + k;
				curr_prog->start = mp_modsub_1(
							curr_prog->start,
							curr_prog->stride_z, 
							p);
			}
		}
	}
}

/*-------------------------------------------------------------------------*/
void
root_sieve_run_deg6(poly_stage2_t *data, double alpha_proj)
{
	uint32 i;
	stage2_curr_data_t *s = (stage2_curr_data_t *)data->internal;
	root_sieve_t *rs = &s->root_sieve;
	curr_poly_t *c = &s->curr_poly;

	uint64 lattice_size_xyz;
	sieve_prime_t lattice_primes[MAX_FACTORS];
	uint32 num_lattice_primes;
	hit_t hitlist[MAX_FACTORS];
	uint32 num_lattices;
	lattice_t *lattices;

	double max_norm = data->max_norm * exp(-alpha_proj);
	double dbl_p = mpz_get_d(c->gmp_p);
	double dbl_d = mpz_get_d(c->gmp_d);
	dpoly_t apoly;
	double direction[3];
	double line_min, line_max;
	int64 z_base;
	uint32 z_blocks;
	xydata_t xydata[MAX_FACTORS];

	apoly.degree = data->degree;
	for (i = 0; i <= data->degree; i++)
		apoly.coeff[i] = mpz_get_d(c->gmp_a[i]);

	direction[0] = 0;
	direction[1] = 0;
	direction[2] = 1;
	compute_line_size(max_norm, &apoly,
			dbl_p, dbl_d, direction,
			-10000, 10000, &line_min, &line_max);

	lattice_size_xyz = find_lattice_size_xyz((line_max - line_min) / 100);

	num_lattice_primes = find_lattice_primes_xyz(rs->primes, 
					rs->num_primes, lattice_size_xyz, 
					lattice_primes);

	find_hits_xyz(lattice_primes, num_lattice_primes, hitlist);

	num_lattices = compute_lattices_xyz(hitlist, 
			num_lattice_primes, &lattices,
			lattice_size_xyz);

	z_base = lattice_size_xyz * (line_min / lattice_size_xyz + 1);
	z_blocks = (line_max - line_min) / lattice_size_xyz;

	direction[0] = 0;
	direction[1] = 1;
	direction[2] = 0;
	compute_line_size(max_norm, &apoly,
			dbl_p, dbl_d, direction,
			-10000, 10000, &line_min, &line_max);

	num_lattice_primes = find_lattice_primes_xy(rs->primes, 
					rs->num_primes, lattice_size_xyz, 
					lattice_primes, line_max - line_min);

	xydata_alloc(lattice_primes, num_lattice_primes, 
			lattice_size_xyz, xydata);

	for (i = 0; i < num_lattices; i++) {

		xydata_init(lattice_primes, xydata, 
				num_lattice_primes,
				lattices + i, z_base);

		find_hits_xy(xydata, num_lattice_primes, 
				lattices + i, z_base, 
				z_blocks, lattice_size_xyz);

		if (i % 512 == 0) {
			printf("."); fflush(stdout);
		}
	}

	printf("num best lattices: %u\n", g_num_max_score);


	xydata_free(xydata, num_lattice_primes);
	free(lattices);
	exit(-1);
}
