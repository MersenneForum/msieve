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
	uint32 x_off;
	uint32 score;
} line_score_t;

#define LINE_HEAP_SIZE 10

typedef struct {
	uint32 num_entries;
	line_score_t entries[LINE_HEAP_SIZE];
} line_score_heap_t;

/*-------------------------------------------------------------------------*/
/* boilerplate code for managing heaps */

#define HEAP_SWAP(a,b) { tmp = a; a = b; b = tmp; }
#define HEAP_PARENT(i)  (((i)-1) >> 1)
#define HEAP_LEFT(i)    (2 * (i) + 1)
#define HEAP_RIGHT(i)   (2 * (i) + 2)

static void
heapify(line_score_t *h, uint32 index, uint32 size) {

	uint32 c;
	line_score_t tmp;
	for (c = HEAP_LEFT(index); c < (size-1); 
			index = c, c = HEAP_LEFT(index)) {

		if (h[c].score > h[c+1].score)
			c++;

		if (h[index].score > h[c].score) {
			HEAP_SWAP(h[index], h[c]);
		}
		else {
			return;
		}
	}
	if (c == (size-1) && h[index].score > h[c].score) {
		HEAP_SWAP(h[index], h[c]);
	}
}

static void
make_heap(line_score_t *h, uint32 size) {

	int32 i;
	for (i = HEAP_PARENT(size); i >= 0; i--)
		heapify(h, (uint32)i, size);
}

static void
save_rotation(line_score_heap_t *heap, uint32 x_off, uint32 score)
{
	line_score_t *h = heap->entries;

	if (heap->num_entries <= LINE_HEAP_SIZE - 1) {
		line_score_t *s = h + heap->num_entries++;
		s->x_off = x_off;
		s->score = score;
		if (heap->num_entries == LINE_HEAP_SIZE)
			make_heap(h, LINE_HEAP_SIZE);
	}
	else if (h->score < score) {
		h->x_off = x_off;
		h->score = score;
		heapify(h, 0, LINE_HEAP_SIZE);
	}
}

/*-------------------------------------------------------------------------*/
static void
sieve_one_block(uint16 *sieve_block, uint32 sieve_block_size,
		sieve_prime_t *primes, uint32 num_primes)
{
	uint32 i, j;
	uint64 *block_array = (uint64 *)sieve_block;
	uint32 block_words = sieve_block_size / UNROLL;

	memset(block_array, 0, block_words * sizeof(uint64));

	for (i = 0; i < num_primes; i++) {
		sieve_prime_t *sp = primes + i;
		uint64 *contrib_array = (uint64 *)sp->contrib_array;
		uint32 contrib_words = sp->contrib_array_size;
		uint32 contrib_offset = sp->contrib_array_offset;
		uint32 block_offset = 0;

		while (block_offset < block_words) {
			uint32 curr_words = MIN(contrib_words - contrib_offset,
						block_words - block_offset);
			uint64 *b = block_array + block_offset;
			uint64 *c = contrib_array + contrib_offset;

#if defined(GCC_ASM32X) && defined(HAS_MMX) && defined(NDEBUG)
			j = 0;
			ASM_G volatile(
			    "cmpl $0, %3              \n\t"
			    "je 1f                   \n\t"
			    ALIGN_LOOP
			    "0:                       \n\t"
			    "movq (%1,%0,8), %%mm0    \n\t"
			    "movq 8(%1,%0,8), %%mm1   \n\t"
			    "movq 16(%1,%0,8), %%mm2  \n\t"
			    "movq 24(%1,%0,8), %%mm3  \n\t"
			    "paddw (%2,%0,8), %%mm0   \n\t"
			    "paddw 8(%2,%0,8), %%mm1  \n\t"
			    "paddw 16(%2,%0,8), %%mm2 \n\t"
			    "paddw 24(%2,%0,8), %%mm3 \n\t"
			    "movq %%mm0, (%1,%0,8)    \n\t"
			    "movq %%mm1, 8(%1,%0,8)   \n\t"
			    "movq %%mm2, 16(%1,%0,8)  \n\t"
			    "movq %%mm3, 24(%1,%0,8)  \n\t"
			    "addl $4, %0              \n\t"
			    "cmpl %3, %0              \n\t"
			    "jne 0b                   \n\t"
			    "1:                       \n\t"
			    :"+r"(j)
			    :"r"(b), "r"(c), "g"(curr_words & (uint32)(~3))
			    :"%mm0", "%mm1", "%mm2", "%mm3", "memory", "cc");
#else
			for (j = 0; j < (curr_words & (uint32)(~3)); j += 4) {
				b[j+0] += c[j+0];
				b[j+1] += c[j+1];
				b[j+2] += c[j+2];
				b[j+3] += c[j+3];
			}
#endif
			for (; j < curr_words; j++)
				b[j] += c[j];

			block_offset += curr_words;
			contrib_offset = mp_modadd_1(contrib_offset,
						curr_words, contrib_words);
		}

		sp->contrib_array_offset = contrib_offset;
	}

#if defined(GCC_ASM32X) && defined(HAS_MMX)
	ASM_G volatile("emms");
#endif
}

/*-------------------------------------------------------------------------*/
static void
prepare_sieve_lattice(root_sieve_t *rs)
{
	uint32 i, j, k, m;
	sieve_prime_t *primes = rs->primes;
	uint32 num_primes = rs->num_primes;
	uint16 *inv = rs->sieve_block;
	sieve_x_t *x = &rs->xdata;

	for (i = 0; i < num_primes; i++) {

		sieve_prime_t *curr_prime = primes + i;
		uint32 num_powers = curr_prime->num_powers;
		uint16 *contrib_array = curr_prime->contrib_array;
		uint32 contrib_array_size = curr_prime->contrib_array_size;

		memset(contrib_array, 0, contrib_array_size * sizeof(uint16));

		for (j = 0; j < num_powers; j++) {

			sieve_power_t *sp = curr_prime->powers + j;
			uint32 power = sp->power;
			uint32 num_roots = sp->num_roots;
			uint16 contrib = sp->sieve_contrib;
			uint32 mod_lattice = mpz_tdiv_ui(x->mp_lattice_size,
							power);
			uint32 xmin_mod;
			uint32 y_mod;
			uint32 z_mod;
			uint32 resclass_mod;
			int32 tmpval;
			uint32 default_step;

			if (mod_lattice == 0)
				continue;

			resclass_mod = mpz_tdiv_ui(x->resclass, power);
			xmin_mod = mpz_fdiv_ui(x->x_base, power);
			y_mod = mpz_fdiv_ui(rs->curr_y, power);
			tmpval = rs->curr_z % power;
			z_mod = (tmpval < 0) ? tmpval + (int32)power : tmpval;

			default_step = power / mp_gcd_1(mod_lattice, power);

			for (k = 0; k < power; k++)
				inv[k] = power;

			for (k = resclass_mod, m = 0; m < power; m++) {
				if (inv[k] == power)
					inv[k] = m;
				k = mp_modadd_1(k, mod_lattice, power);
			}

			for (k = 0; k < num_roots; k++) {
				sieve_root_t *r = sp->roots + k;
				uint32 start;
				uint32 step = r->step;
				uint32 resclass = r->resclass;

				start = mp_modadd_1(y_mod, 
						    resclass * z_mod % power, 
						    power);
				start = mp_modsub_1(r->start, 
						    resclass * start % power, 
						    power);
				start = mp_modsub_1(start, xmin_mod, power);

				if (step == power) {
					start = inv[start];
					if (start == power)
						continue;
					step = default_step;
				}
				else {
					uint32 mult = power / step;

					start = start % step;
					start = inv[mult * start];
					if (start == power)
						continue;
					step /= mp_gcd_1(mod_lattice, step);
					start = start % step;
				}

				while (start < contrib_array_size) {
					contrib_array[start] += contrib;
					start += step;
				}
			}
		}

		curr_prime->contrib_array_offset = 0;

		for (j = 1; j < UNROLL; j++) {
			memcpy(contrib_array + j * contrib_array_size, 
				contrib_array, 
				contrib_array_size * sizeof(uint16));
		}
	}
}

/*-------------------------------------------------------------------------*/
void
root_sieve_line(root_sieve_t *rs)
{
	uint32 i, j;
	sieve_x_t *x = &rs->xdata;
	uint32 worst_score = 0;
	line_score_heap_t line_heap;
	uint32 num_blocks = x->x_blocks / DEFAULT_BLOCK_SIZE + 1;

	sieve_prime_t *primes = rs->primes;
	uint32 num_primes = rs->num_primes;
	uint16 *block = rs->sieve_block;

	prepare_sieve_lattice(rs);

	line_heap.num_entries = 0;

	for (i = 0; i < num_blocks; i++) {

		sieve_one_block(block, DEFAULT_BLOCK_SIZE,
				primes, num_primes);

		for (j = 0; j < DEFAULT_BLOCK_SIZE; j++) {
			uint32 score = block[j];
			if (score >= worst_score) {

				save_rotation(&line_heap, 
						DEFAULT_BLOCK_SIZE * i + j,
						score);

				if (line_heap.num_entries == LINE_HEAP_SIZE) {
					worst_score = line_heap.entries[
						      0].score;
				}
			}
		}
	}

	for (i = 0; i < line_heap.num_entries; i++) {
		mpz_add(rs->curr_x, x->x_base, x->resclass);
		mpz_addmul_ui(rs->curr_x, x->mp_lattice_size,
				line_heap.entries[i].x_off);
		optimize_final(rs->curr_x, rs->curr_y, rs->curr_z,
				(poly_stage2_t *)rs->root_heap.extra);
	}
}
