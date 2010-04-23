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
	uint32 stride_y;
} xprog_t;

typedef struct {
	uint32 p;
	uint32 latsize_mod_p;
	uint32 num_roots;
	uint32 max_sieve_val;
	uint16 contrib;
	xprog_t *roots;
	uint8 *sieve;
} xdata_t;

/*-------------------------------------------------------------------------*/
static uint32
find_lattice_primes(sieve_prime_t *primes, uint32 num_primes,
			mpz_t mp_lattice_size, sieve_prime_t *lattice_primes,
			uint64 *lattice_size_x, double line_length)
{
	uint32 i;
	uint32 num_lattice_primes = 0;
	uint64 tmp = 1;
	double target_size = line_length / mpz_get_d(mp_lattice_size) / 1e8;
	double curr_size = 1.0;

	for (i = 0; i < num_primes; i++) {

		sieve_prime_t *curr_prime = primes + i;

		if (num_lattice_primes == MAX_CRT_FACTORS)
			break;

		if (mpz_tdiv_ui(mp_lattice_size, curr_prime->prime) == 0)
			continue;

		if (curr_size * curr_prime->prime >= target_size)
			break;

		tmp *= curr_prime->prime;
		lattice_primes[num_lattice_primes] = *curr_prime;
		lattice_primes[num_lattice_primes].num_powers = 1;
		curr_size *= curr_prime->prime;
		num_lattice_primes++;
	}

	*lattice_size_x = tmp;
	return num_lattice_primes;
}

/*-------------------------------------------------------------------------*/
#define MAX_X_LATTICES 30

static void 
root_sieve_x(root_sieve_t *rs, xdata_t *xdata, 
		uint32 num_lattice_primes, uint32 which_yblock)
{
	uint32 i, j, k;
	hit_t hitlist[MAX_CRT_FACTORS];
	lattice_t lattices_x[MAX_X_LATTICES];
	sieve_xy_t *xy = &rs->xydata;
	sieve_x_t *x = &rs->xdata;
	uint32 num_lattices = 1;

	double line_min, line_max;
	double direction[3] = {1, 0, 0};

	for (i = 0; i < num_lattice_primes; i++) {

		xdata_t *curr_xdata = xdata + i;
		uint32 p = curr_xdata->p;
		uint8 *sieve = curr_xdata->sieve;
		uint32 max_sieve_val = curr_xdata->max_sieve_val;
		uint16 contrib = curr_xdata->contrib;
		hit_t *hits = hitlist + i;

		for (j = k = 0; j < p; j++) {
			if (sieve[j] >= max_sieve_val) {
				hits->score[k] = sieve[j] * contrib;
				hits->roots[k][0] = j;
				hits->roots[k][1] = 0;
				hits->roots[k][2] = 0;
				k++;
			}
		}

		hits->power = p;
		hits->num_roots = k;
		num_lattices *= k;
	}

	num_lattices = MIN(num_lattices, MAX_X_LATTICES);

	compute_lattices(hitlist, num_lattice_primes,
			lattices_x, x->lattice_size, num_lattices);

	mpz_add(rs->curr_y, xy->y_base, xy->resclass_y);
	mpz_addmul_ui(rs->curr_y, xy->mp_lattice_size, which_yblock);

	x->apoly = xy->apoly;
	x->apoly.coeff[2] += mpz_get_d(rs->curr_y) * rs->dbl_p;
	x->apoly.coeff[1] -= mpz_get_d(rs->curr_y) * rs->dbl_d;
	compute_line_size_deg6(rs->max_norm, &x->apoly,
			rs->dbl_p, rs->dbl_d, direction,
			x->last_line_min, x->last_line_max,
			&line_min, &line_max);

	mpz_set_d(x->tmp1, 0.01 * line_min);
	mpz_tdiv_q(x->x_base, x->tmp1, x->mp_lattice_size);
	mpz_mul(x->x_base, x->x_base, x->mp_lattice_size);
	x->x_blocks = 0.01 * (line_max - line_min) / x->dbl_lattice_size;

	for (i = 0; i < num_lattices; i++) {

		lattice_t *lattice_x = lattices_x + i;

		uint64_2gmp(lattice_x->x, x->tmp1);
		mpz_mul(x->resclass, x->tmp1, x->crt0);
		mpz_addmul(x->resclass, xy->resclass_x, x->crt1);

		mpz_tdiv_r(x->resclass, x->resclass, x->mp_lattice_size);

		x->curr_score = xy->curr_score + lattice_x->score;

		root_sieve_line(rs);
	}
}

/*-------------------------------------------------------------------------*/
static void
xdata_alloc(sieve_prime_t *lattice_primes, 
		uint32 num_lattice_primes, 
		mpz_t mp_lattice_size,
		xdata_t *xdata)
{
	uint32 i;

	for (i = 0; i < num_lattice_primes; i++) {

		sieve_prime_t *curr_prime = lattice_primes + i;
		xdata_t *curr_xdata = xdata + i;
		uint32 p = curr_prime->prime;
		uint32 num_roots = curr_prime->powers[0].num_roots;

		curr_xdata->p = p;
		curr_xdata->latsize_mod_p = mpz_tdiv_ui(mp_lattice_size, p);
		curr_xdata->num_roots = num_roots;
		curr_xdata->max_sieve_val = MIN(num_roots, 5);
		curr_xdata->contrib = curr_prime->powers[0].sieve_contrib;
		curr_xdata->roots = (xprog_t *)xmalloc(num_roots * 
							sizeof(xprog_t));
		curr_xdata->sieve = (uint8 *)xmalloc(p * sizeof(uint8));
	}
}

/*-------------------------------------------------------------------------*/
static void
xdata_free(xdata_t *xdata, uint32 num_lattice_primes)
{
	uint32 i;

	for (i = 0; i < num_lattice_primes; i++) {

		xdata_t *curr_xdata = xdata + i;

		free(curr_xdata->roots);
		free(curr_xdata->sieve);
	}
}

/*-------------------------------------------------------------------------*/
static void 
xdata_init(sieve_prime_t *lattice_primes, xdata_t *xdata,
		uint32 num_lattice_primes, 
		mpz_t y_base, mpz_t resclass_y, int64 curr_z)
{
	uint32 i, j;

	for (i = 0; i < num_lattice_primes; i++) {

		sieve_prime_t *curr_prime = lattice_primes + i;
		sieve_power_t *sp = curr_prime->powers + 0;
		xdata_t *curr_xdata = xdata + i;

		uint32 p = curr_xdata->p;
		uint32 num_roots = curr_xdata->num_roots;
		xprog_t *roots = curr_xdata->roots;

		uint32 latsize_mod_p = curr_xdata->latsize_mod_p;
		uint32 yres_mod_p = mpz_tdiv_ui(resclass_y, p);
		uint32 ybase_mod_p = mpz_fdiv_ui(y_base, p);
		uint32 y_mod_p = mp_modadd_1(yres_mod_p, ybase_mod_p, p);
		int32 z_mod = curr_z % p;
		uint32 z_mod_p = (z_mod < 0) ?  (z_mod + (int32)p) : z_mod;

		for (j = 0; j < num_roots; j++) {
			sieve_root_t *r = sp->roots + j;
			xprog_t *curr_xprog = roots + j;
			uint32 start = r->start;
			uint32 resclass = r->resclass;
			uint32 resclass2 = mp_modmul_1(resclass, resclass, p);

			curr_xprog->stride_y = mp_modmul_1(resclass, 
							latsize_mod_p, p);
			start = mp_modsub_1(start, 
					mp_modmul_1(resclass, y_mod_p, p), 
					p);
			curr_xprog->start = mp_modsub_1(start, 
					mp_modmul_1(resclass2, z_mod_p, p), 
					p);
		}
	}
}

/*-------------------------------------------------------------------------*/
static void 
find_hits(root_sieve_t *rs, xdata_t *xdata, 
		uint32 num_lattice_primes, uint32 y_blocks)
{
	uint32 i, j, k;

	for (i = 0; i < y_blocks; i++) {

		for(j = 0; j < num_lattice_primes; j++) {

			xdata_t *curr_xdata = xdata + j;
			uint32 p = curr_xdata->p;
			uint32 num_roots = curr_xdata->num_roots;
			uint32 max_sieve_val = curr_xdata->max_sieve_val;
			xprog_t *roots = curr_xdata->roots;
			uint8 *sieve = curr_xdata->sieve;

			memset(sieve, 0, p * sizeof(uint8));

			for (k = 0; k < num_roots; k++) {

				xprog_t *curr_prog = roots + k;

				sieve[curr_prog->start]++;
			}

			for (k = 0; k < p; k++) {
				if (sieve[k] >= max_sieve_val)
					break;
			}
			if (k == p)
				break;
		}

		if (j == num_lattice_primes) {
			root_sieve_x(rs, xdata, num_lattice_primes, i);
		}

		for (j = 0; j < num_lattice_primes; j++) {

			xdata_t *curr_xdata = xdata + j;
			uint32 p = curr_xdata->p;
			uint32 num_roots = curr_xdata->num_roots;
			xprog_t *roots = curr_xdata->roots;

			for (k = 0; k < num_roots; k++) {

				xprog_t *curr_prog = roots + k;
				curr_prog->start = mp_modsub_1(
							curr_prog->start,
							curr_prog->stride_y, p);
			}
		}
	}
}

/*-------------------------------------------------------------------------*/
void
sieve_x_alloc(sieve_x_t *x)
{
	mpz_init(x->mp_lattice_size);
	mpz_init(x->resclass);
	mpz_init(x->crt0);
	mpz_init(x->crt1);
	mpz_init(x->tmp1);
	mpz_init(x->tmp2);
	mpz_init(x->tmp3);
	mpz_init(x->tmp4);
}

/*-------------------------------------------------------------------------*/
void
sieve_x_free(sieve_x_t *x)
{
	mpz_clear(x->mp_lattice_size);
	mpz_clear(x->resclass);
	mpz_clear(x->crt0);
	mpz_clear(x->crt1);
	mpz_clear(x->tmp1);
	mpz_clear(x->tmp2);
	mpz_clear(x->tmp3);
	mpz_clear(x->tmp4);
}

/*-------------------------------------------------------------------------*/
void
sieve_x_run(root_sieve_t *rs)
{
	sieve_xy_t *xy = &rs->xydata;

	sieve_x_t *x = &rs->xdata;
	sieve_prime_t *lattice_primes = x->lattice_primes;
	uint32 num_lattice_primes;

	double direction[3] = {1, 0, 0};
	double line_min, line_max;
	xdata_t xdata[MAX_CRT_FACTORS];

	compute_line_size_deg6(rs->max_norm, &xy->apoly,
			rs->dbl_p, rs->dbl_d, direction,
			-10000, 10000, &line_min, &line_max);

	num_lattice_primes = x->num_lattice_primes = 
				find_lattice_primes(rs->primes, 
					rs->num_primes, xy->mp_lattice_size, 
					lattice_primes, &x->lattice_size,
					line_max - line_min);
	x->last_line_min = line_min;
	x->last_line_max = line_max;

	uint64_2gmp(x->lattice_size, x->tmp1);
	mpz_invert(x->crt0, xy->mp_lattice_size, x->tmp1);
	mpz_invert(x->crt1, x->tmp1, xy->mp_lattice_size);
	mpz_mul(x->crt0, x->crt0, xy->mp_lattice_size);
	mpz_mul(x->crt1, x->crt1, x->tmp1);
	mpz_mul(x->mp_lattice_size, xy->mp_lattice_size, x->tmp1);
	x->dbl_lattice_size = mpz_get_d(x->mp_lattice_size);

	xdata_alloc(lattice_primes, num_lattice_primes, 
			xy->mp_lattice_size, xdata);

	xdata_init(lattice_primes, xdata, num_lattice_primes,
			xy->y_base, xy->resclass_y, 
			rs->curr_z);

	find_hits(rs, xdata, num_lattice_primes, xy->y_blocks);

	xdata_free(xdata, num_lattice_primes);
}
