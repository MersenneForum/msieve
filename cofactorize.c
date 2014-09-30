#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "alphad.h"
#include "alphag.h"
#include "mpdim.h"
#include "functions.h"
#include "cofactorize.h"

static mpz_t rat_max_1LP;
static mpz_t rat_min_2LP;
static mpz_t rat_max_2LP;
static mpz_t rat_min_3LP;
static mpz_t rat_max_3LP;

static mpz_t alg_max_1LP;
static mpz_t alg_min_2LP;
static mpz_t alg_max_2LP;
static mpz_t alg_min_3LP;
static mpz_t alg_max_3LP;

static mpz_t n, t1, t2, t3, fac1, fac2;

void cofac_setup()
{
  mpz_t tmp;

  mpz_init(tmp);
  mpz_init(n);
  mpz_init(t1);
  mpz_init(t2);
  mpz_init(t3);
  mpz_init(fac1);
  mpz_init(fac2);

  /* rational bounds */
  mpz_init(rat_max_1LP);
  u64_2gmp(BPRIME_LIMIT_R, rat_max_1LP);

  /* min for 2 large primes is the square of the FB limit */
  mpz_init_set_ui(rat_min_2LP, int_pmax);
  mpz_mul(rat_min_2LP, rat_min_2LP, rat_min_2LP);

  /* the max for 2 large primes is (large_prime_bound)^1.9 */
  u64_2gmp(BPRIME_LIMIT_R, tmp);
  mpz_init_set_d(rat_max_2LP, pow((double)BPRIME_LIMIT_R, 0.9));
  mpz_mul(rat_max_2LP, rat_max_2LP, tmp);

  /* if configured for 3LP, the min for 3 large primes 
    is (FB limit)^3.2 and the max is (large_prime_bound)^2.8 */
  if (num_lp_r == 3) 
  {
    mpz_set_ui(tmp, int_pmax);
    mpz_init_set_d(rat_min_3LP, pow((double)int_pmax, 1.2));
    mpz_mul(tmp, tmp, tmp);
    mpz_mul(rat_min_3LP, rat_min_3LP, tmp);

    u64_2gmp(BPRIME_LIMIT_R, tmp);
    mpz_init_set_d(rat_max_3LP, pow((double)BPRIME_LIMIT_R, 0.8));
    mpz_mul(tmp, tmp, tmp);
    mpz_mul(rat_max_3LP, rat_max_3LP, tmp);
  }

  /* algebraic bounds */
  mpz_init(alg_max_1LP);
  u64_2gmp(BPRIME_LIMIT_A, alg_max_1LP);

  /* min for 2 large primes is the square of the FB limit */
  mpz_init_set_ui(alg_min_2LP, alg_pmax);
  mpz_mul(alg_min_2LP, alg_min_2LP, alg_min_2LP);

  /* the max for 2 large primes is (large_prime_bound)^1.9 */
  u64_2gmp(BPRIME_LIMIT_A, tmp);
  mpz_init_set_d(alg_max_2LP, pow((double)BPRIME_LIMIT_A, 0.9));
  mpz_mul(alg_max_2LP, alg_max_2LP, tmp);

  /* if configured for 3LP, the min for 3 large primes 
    is (FB limit)^3.2 and the max is (large_prime_bound)^2.8 */
  if (num_lp_a == 3) 
  {
    mpz_set_ui(tmp, alg_pmax);
    mpz_init_set_d(alg_min_3LP, pow((double)alg_pmax, 1.2));
    mpz_mul(tmp, tmp, tmp);
    mpz_mul(alg_min_3LP, alg_min_3LP, tmp);

    u64_2gmp(BPRIME_LIMIT_A, tmp);
    mpz_init_set_d(alg_max_3LP, pow((double)BPRIME_LIMIT_A, 0.8));
    mpz_mul(tmp, tmp, tmp);
    mpz_mul(alg_max_3LP, alg_max_3LP, tmp);
  }

  mpz_clear(tmp);
}

static void mp2gmp(int *x, mpz_t out)
{
  mpz_import(out, SIZE(x), -1, sizeof(int), 0, 2, &x[1]);
}

static void gmp2mp(mpz_t x, int *out)
{
  mpz_export(&out[1], &out[0], -1, 0, sizeof(int), 2, x);
}

int cofactorize(int *mp_n, int *LP1, int *LP2, int *LP3, int is_rat)
{
  mpz_t *max_1LP = &rat_max_1LP;
  mpz_t *min_2LP = &rat_min_2LP;
  mpz_t *max_2LP = &rat_max_2LP;
  mpz_t *min_3LP = &rat_min_3LP;
  mpz_t *max_3LP = &rat_max_3LP;
  int num_LP = num_lp_r;

  if (!is_rat) 
  {
    mpz_t *max_1LP = &alg_max_1LP;
    mpz_t *min_2LP = &alg_min_2LP;
    mpz_t *max_2LP = &alg_max_2LP;
    mpz_t *min_3LP = &alg_min_3LP;
    mpz_t *max_3LP = &alg_max_3LP;
    num_LP = num_lp_a;
  }

  SIZE(LP1) = SIZE(LP2) = SIZE(LP3) = 0;
  mp2gmp(mp_n, n);

  /* less than the 1LP bound: success */
  if (mpz_cmp(n, *max_1LP) <= 0)
  {
    small_mpcopy(mp_n, LP1);
    return 1;
  }

  /* not within the 2LP or 3LP bounds: failure */
  if (mpz_cmp(n, *min_2LP) < 0 ||
      num_LP == 2 && mpz_cmp(n, *max_2LP) > 0 ||
      num_LP == 3 && ((mpz_cmp(n, *max_2LP) > 0 &&
                       mpz_cmp(n, *min_3LP) < 0) ||
                      mpz_cmp(n, *max_3LP) > 0))
  {
    return 0;
  }

  /* perform a compositeness test on n 
     (check whether 2^(n-1) mod n != 1) */

  mpz_set_ui(t1, 2);
  mpz_sub_ui(t2, n, 1);
  mpz_powm(t3, t1, t2, n);
  if (mpz_cmp_ui(t3, 1) == 0)
  {
    return 0;
  }

  /* attempt to split n */

  if (mpz_sizeinbase(n, 2) <= 58)
  {
    u32 f = squfof(n);
    if (f < 2)
    {
      return 0;
    }
    mpz_set_ui(fac1, f);
    mpz_divexact_ui(fac2, n, f);
  }
  else
  {
    if (tinyqs(n, fac1, fac2) == 0)
    {
      return 0;
    }
  }

  /* split succeeded; test for 2LP relation */

  if (mpz_cmp(fac1, *max_1LP) <= 0 &&
      mpz_cmp(fac2, *max_1LP) <= 0)
  {
    gmp2mp(fac1, LP1);
    gmp2mp(fac2, LP2);
    return 1;
  }

  /* the smallest factor has to be below the LP bound */

  if (mpz_cmp(fac1, fac2) > 0)
  {
    mpz_swap(fac1, fac2);
  }

  if (mpz_cmp(fac1, *max_1LP) <= 0)
  {
    /* the larger factor needs to be in the 2LP range and
       to be composite */

    if (mpz_cmp(fac2, *min_2LP) >= 0 &&
        mpz_cmp(fac2, *max_2LP) <= 0)
    {
      mpz_set_ui(t1, 2);
      mpz_sub_ui(t2, fac2, 1);
      mpz_powm(t3, t1, t2, fac2);
      if (mpz_cmp_ui(t3, 1) != 0)
      {
        /* attempt to split the larger factor */

        if (mpz_sizeinbase(fac2, 2) <= 58)
        {
          u32 f = squfof(fac2);
          if (f < 2)
          {
            return 0;
          }
          mpz_set_ui(t1, f);
          mpz_divexact_ui(t2, fac2, f);
        }
        else
        {
          if (tinyqs(fac2, t1, t2) == 0)
          {
            return 0;
          }
        }

        /* split succeeded; if both factors are less than the
           large prime bound then we found a 3LP relation */

        if (mpz_cmp(t1, *max_1LP) <= 0 &&
            mpz_cmp(t2, *max_1LP) <= 0)
        {
          gmp2mp(fac1, LP1);
          gmp2mp(t1, LP2);
          gmp2mp(t2, LP3);
          return 1;
        }
      }
    }
  } 

  return 0;
}
