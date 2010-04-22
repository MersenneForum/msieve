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

/*-------------------------------------------------------------------------*/
static uint64
find_lattice_size(double line_length)
{
	return (uint64)2*2*2*2*2*3*3*3*5*5*7;
}

/*-------------------------------------------------------------------------*/
static uint32
find_lattice_primes(sieve_prime_t *primes, uint32 num_primes,
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
static void 
do_sieving(sieve_root_t *r, uint16 *sieve,
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
#define MAX_DIM 64

static void 
find_hits(sieve_prime_t *lattice_primes,
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
				do_sieving(sp->roots + k, sieve, 
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

				hits->score[k] = curr_score;
				hits->roots[k][0] = j % dim;
				hits->roots[k][1] = (j / dim) % dim;
				hits->roots[k][2] = (j / dim) / dim;
				k++;
			}
		}
		hits->power = dim;
		hits->num_roots = k;
	}

	free(sieve);
}

/*-------------------------------------------------------------------------*/
void
sieve_xyz_alloc(sieve_xyz_t *xyz)
{
	memset(xyz, 0, sizeof(sieve_xyz_t));
}

/*-------------------------------------------------------------------------*/
void 
sieve_xyz_free(sieve_xyz_t *xyz)
{
	free(xyz->lattices);
}

/*-------------------------------------------------------------------------*/
void
sieve_xyz_run(root_sieve_t *rs)
{
	uint32 i;
	sieve_xyz_t *xyz = &rs->xyzdata;
	hit_t hitlist[MAX_CRT_FACTORS];
	uint64 lattice_size;
	uint32 num_lattice_primes;
	uint32 num_lattices;

	double direction[3] = {0, 0, 1};
	double line_min, line_max;

	compute_line_size_deg6(rs->max_norm, &rs->apoly,
			rs->dbl_p, rs->dbl_d, direction,
			-10000, 10000, &line_min, &line_max);

	lattice_size = xyz->lattice_size = 
			find_lattice_size((line_max - line_min) / 100);

	xyz->z_base = line_min / lattice_size - 1;
	xyz->z_base *= lattice_size;
	xyz->z_blocks = (line_max - line_min) / lattice_size;

	num_lattice_primes = xyz->num_lattice_primes = 
			find_lattice_primes(rs->primes, 
					rs->num_primes, lattice_size, 
					xyz->lattice_primes);

	find_hits(xyz->lattice_primes, num_lattice_primes, hitlist);

	for (i = 0, num_lattices = 1; i < num_lattice_primes; i++) {
		num_lattices *= hitlist[i].num_roots;
	}

	if (num_lattices > xyz->num_lattices) {
		xyz->lattices = (lattice_t *)xrealloc(xyz->lattices,
					num_lattices * sizeof(lattice_t));
	}
	xyz->num_lattices = num_lattices;

	compute_lattices(hitlist, num_lattice_primes, xyz->lattices,
			lattice_size, num_lattices);
}