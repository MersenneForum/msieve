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
	uint32 start;
	uint32 stride_z;
} xyprog_t;

typedef struct {
	uint32 p;
	uint32 latsize_mod_p;
	uint32 table_size;
	uint32 num_roots;
	uint32 max_sieve_val;
	uint16 contrib;
	xyprog_t *roots;
	uint8 *invtable_y;
	uint8 *sieve;
} xydata_t;

/*-------------------------------------------------------------------------*/
static uint32
find_lattice_primes(sieve_prime_t *primes, uint32 num_primes,
			uint64 lattice_size_xyz, sieve_prime_t *lattice_primes,
			uint64 *lattice_size_xy, double line_length)
{
	uint32 i;
	uint32 num_lattice_primes = 0;
	uint64 tmp = 1;
	double target_size = line_length / (double)lattice_size_xyz / 1e4;
	double curr_size = 1.0;

	for (i = 0; i < num_primes; i++) {

		sieve_prime_t *curr_prime = primes + i;

		if (num_lattice_primes == MAX_CRT_FACTORS)
			break;

		if (curr_prime->powers[0].num_roots == 0 ||
		    lattice_size_xyz % curr_prime->prime == 0)
			continue;

		if (curr_size * curr_prime->prime >= target_size)
			break;

		tmp *= curr_prime->prime;
		lattice_primes[num_lattice_primes] = *curr_prime;
		lattice_primes[num_lattice_primes].num_powers = 1;
		curr_size *= curr_prime->prime;
		num_lattice_primes++;
	}

	*lattice_size_xy = tmp;
	return num_lattice_primes;
}

/*-------------------------------------------------------------------------*/
#define MAX_XY_LATTICES 10

static void 
root_sieve_xy(root_sieve_t *rs, xydata_t *xydata, 
		uint32 num_lattice_primes, 
		lattice_t *lattice_xyz, int64 z_base,
		double line_min, double line_max)
{
	uint32 i, j, k;
	hit_t hitlist[MAX_CRT_FACTORS];
	lattice_t lattices_xy[MAX_XY_LATTICES];
	sieve_xy_t *xy = &rs->xydata;
	uint32 num_lattices = 1;

	for (i = 0; i < num_lattice_primes; i++) {

		xydata_t *curr_xydata = xydata + i;
		uint32 p = curr_xydata->p;
		uint32 table_size = curr_xydata->table_size;
		uint8 *sieve = curr_xydata->sieve;
		uint32 max_sieve_val = curr_xydata->max_sieve_val;
		uint16 contrib = curr_xydata->contrib;
		hit_t *hits = hitlist + i;

		for (j = k = 0; j < table_size; j++) {
			if (sieve[j] >= max_sieve_val) {
				hits->score[k] = sieve[j] * contrib;
				hits->roots[k][0] = j % p;
				hits->roots[k][1] = j / p;
				k++;
			}
		}

		hits->power = p;
		hits->num_roots = k;
		num_lattices *= k;
	}

	num_lattices = MIN(num_lattices, MAX_XY_LATTICES);

	compute_lattices(hitlist, num_lattice_primes,
			lattices_xy, xy->lattice_size, num_lattices, 2);

	xy->apoly = rs->apoly;
	xy->apoly.coeff[3] += z_base * rs->dbl_p;
	xy->apoly.coeff[2] -= z_base * rs->dbl_d;

	mpz_set_d(xy->tmp1, line_min);
	mpz_tdiv_q(xy->y_base, xy->tmp1, xy->mp_lattice_size);
	mpz_mul(xy->y_base, xy->y_base, xy->mp_lattice_size);
	xy->y_blocks = (line_max - line_min) / xy->dbl_lattice_size;

	for (i = 0; i < num_lattices; i++) {

		lattice_t *lattice_xy = lattices_xy + i;

		uint64_2gmp(lattice_xy->x, xy->tmp1);
		uint64_2gmp(lattice_xyz->x, xy->tmp2);
		mpz_mul(xy->resclass_x, xy->tmp1, xy->crt0);
		mpz_addmul(xy->resclass_x, xy->tmp2, xy->crt1);

		uint64_2gmp(lattice_xy->y, xy->tmp1);
		uint64_2gmp(lattice_xyz->y, xy->tmp2);
		mpz_mul(xy->resclass_y, xy->tmp1, xy->crt0);
		mpz_addmul(xy->resclass_y, xy->tmp2, xy->crt1);

		mpz_tdiv_r(xy->resclass_x, xy->resclass_x, 
				xy->mp_lattice_size);
		mpz_tdiv_r(xy->resclass_y, xy->resclass_y, 
				xy->mp_lattice_size);

		xy->curr_score = lattice_xyz->score + lattice_xy->score;
		rs->curr_z = z_base + lattice_xyz->z; 

		sieve_x_run(rs);
	}
}

/*-------------------------------------------------------------------------*/
static void
xydata_alloc(sieve_prime_t *lattice_primes, 
		uint32 num_lattice_primes, 
		uint64 lattice_size_xyz,
		xydata_t *xydata,
		uint32 num_roots_min)
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
		curr_xydata->max_sieve_val = MIN(num_roots - 1, num_roots_min);
		curr_xydata->contrib = curr_prime->powers[0].sieve_contrib;
		curr_xydata->roots = (xyprog_t *)xmalloc(num_roots * 
							sizeof(xyprog_t));
		curr_xydata->sieve = (uint8 *)xmalloc((2 * table_size) * 
							sizeof(uint8));
		curr_xydata->invtable_y = curr_xydata->sieve + table_size;
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
		free(curr_xydata->sieve);
	}
}

/*-------------------------------------------------------------------------*/
static void 
xydata_init(sieve_prime_t *lattice_primes, xydata_t *xydata,
		uint32 num_lattice_primes, 
		lattice_t *lattice_xyz, int64 z_base)
{
	uint32 i, j, k, m;

	for (i = 0; i < num_lattice_primes; i++) {

		sieve_prime_t *curr_prime = lattice_primes + i;
		sieve_power_t *sp = curr_prime->powers + 0;
		xydata_t *curr_xydata = xydata + i;

		uint32 p = curr_xydata->p;
		uint32 num_roots = curr_xydata->num_roots;
		uint8 *invtable_y = curr_xydata->invtable_y;
		xyprog_t *roots = curr_xydata->roots;

		uint32 latsize_mod_p = curr_xydata->latsize_mod_p;
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
			uint32 ytmp = y_mod_p;
			uint32 stride_y = mp_modmul_1(resclass, 
							latsize_mod_p, p);

			curr_xyprog->stride_z = mp_modmul_1(resclass2, 
							latsize_mod_p, p);

			start = mp_modsub_1(start, 
					mp_modmul_1(resclass, y_mod_p, p), 
					p);
			curr_xyprog->start = mp_modsub_1(start, 
					mp_modmul_1(resclass2, z_mod_p, p), 
					p);

			for (k = m = 0; k < p; k++) {
				invtable_y[ytmp] = m;
				ytmp = mp_modadd_1(ytmp, latsize_mod_p, p);
				m = mp_modadd_1(m, stride_y, p);
			}
			invtable_y += p;
		}
	}
}

/*-------------------------------------------------------------------------*/
static void 
find_hits(root_sieve_t *rs, xydata_t *xydata, 
		uint32 num_lattice_primes, 
		lattice_t *lattice_xyz)
{
	uint32 i, j, k, m;
	sieve_xyz_t *xyz = &rs->xyzdata;
	uint64 lattice_size_xyz = xyz->lattice_size;
	int64 z_base = xyz->z_base;
	uint32 z_blocks = xyz->z_blocks;

	for (i = 0; i < z_blocks; i++) {

		for(j = 0; j < num_lattice_primes; j++) {

			xydata_t *curr_xydata = xydata + j;
			uint32 p = curr_xydata->p;
			uint32 table_size = curr_xydata->table_size;
			uint32 num_roots = curr_xydata->num_roots;
			uint32 max_sieve_val = curr_xydata->max_sieve_val;
			xyprog_t *roots = curr_xydata->roots;
			uint8 *sieve = curr_xydata->sieve;
			uint8 *invtable_y = curr_xydata->invtable_y;

			memset(sieve, 0, table_size * sizeof(uint8));

			for (k = 0; k < num_roots; k++) {

				uint8 *row = sieve;
				xyprog_t *curr_prog = roots + k;
				uint32 start = curr_prog->start;

				for (m = 0; m < p; m++) {
					uint32 curr_start = mp_modsub_1(start, 
							invtable_y[m], p);
					row[curr_start]++;
					row += p;
				}
				invtable_y += p;
			}

			for (k = 0; k < table_size; k++) {
				if (sieve[k] >= max_sieve_val)
					break;
			}
			if (k == table_size)
				break;
		}

		if (j == num_lattice_primes) {
			root_sieve_xy(rs, xydata, 
					num_lattice_primes,
					lattice_xyz, 
					z_base + i * lattice_size_xyz,
					xyz->y_line_min[i],
					xyz->y_line_max[i]);
		}

		for (j = 0; j < num_lattice_primes; j++) {

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
sieve_xy_alloc(sieve_xy_t *xy)
{
	mpz_init(xy->y_base);
	mpz_init(xy->mp_lattice_size);
	mpz_init(xy->resclass_x);
	mpz_init(xy->resclass_y);
	mpz_init(xy->crt0);
	mpz_init(xy->crt1);
	mpz_init(xy->tmp1);
	mpz_init(xy->tmp2);
	mpz_init(xy->tmp3);
	mpz_init(xy->tmp4);
}

/*-------------------------------------------------------------------------*/
void
sieve_xy_free(sieve_xy_t *xy)
{
	mpz_clear(xy->y_base);
	mpz_clear(xy->mp_lattice_size);
	mpz_clear(xy->resclass_x);
	mpz_clear(xy->resclass_y);
	mpz_clear(xy->crt0);
	mpz_clear(xy->crt1);
	mpz_clear(xy->tmp1);
	mpz_clear(xy->tmp2);
	mpz_clear(xy->tmp3);
	mpz_clear(xy->tmp4);
}

/*-------------------------------------------------------------------------*/
void
sieve_xy_run(root_sieve_t *rs)
{
	uint32 i;

	poly_stage2_t *data = (poly_stage2_t *)rs->root_heap.extra;
	sieve_xyz_t *xyz = &rs->xyzdata;
	int64 z_base = xyz->z_base;

	sieve_xy_t *xy = &rs->xydata;
	sieve_prime_t *lattice_primes = xy->lattice_primes;
	uint32 num_lattice_primes;

	double direction[3] = {0, 1, 0};
	double line_min, line_max;
	xydata_t xydata[MAX_CRT_FACTORS];

	uint64 inv_xy;
	uint64 inv_xyz;

	compute_line_size_deg6(rs->max_norm, &rs->apoly,
			rs->dbl_p, rs->dbl_d, direction,
			-10000, 10000, &line_min, &line_max);
	if (line_min > line_max)
		return;

	num_lattice_primes = xy->num_lattice_primes = 
				find_lattice_primes(rs->primes, 
					rs->num_primes, xyz->lattice_size, 
					lattice_primes, &xy->lattice_size,
					line_max - line_min);

	inv_xy = mp_modinv_2(xyz->lattice_size, xy->lattice_size);
	inv_xyz = mp_modinv_2(xy->lattice_size, xyz->lattice_size);
	uint64_2gmp(xy->lattice_size, xy->tmp1);
	uint64_2gmp(inv_xy, xy->tmp2);
	uint64_2gmp(xyz->lattice_size, xy->tmp3);
	uint64_2gmp(inv_xyz, xy->tmp4);
	mpz_mul(xy->mp_lattice_size, xy->tmp1, xy->tmp3);
	mpz_mul(xy->crt0, xy->tmp2, xy->tmp3);
	mpz_mul(xy->crt1, xy->tmp1, xy->tmp4);
	xy->dbl_lattice_size = mpz_get_d(xy->mp_lattice_size);

	xydata_alloc(lattice_primes, num_lattice_primes, 
			xyz->lattice_size, xydata, xy->num_roots_min);

	for (i = 0; i < xyz->num_lattices; i++) {

		lattice_t *curr_lattice_xyz = xyz->lattices + i;

		xydata_init(lattice_primes, xydata, 
				num_lattice_primes,
				curr_lattice_xyz, z_base);

		find_hits(rs, xydata, num_lattice_primes, 
				curr_lattice_xyz);

		if (data->obj->flags & MSIEVE_FLAG_STOP_SIEVING)
			break;
	}

	xydata_free(xydata, num_lattice_primes);
}
