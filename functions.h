/****************************************************************/
/*																*/
/*	Function templates for lattice siever						*/
/*																*/
/****************************************************************/
void
	mul_single(int *a, int b, int *c),
	get_file(int *a, FILE *ptr),
	write_num(int *a, char *str),
	add(int *a, int *b, int *c),
	get_logs(int *a, int b, unsigned char *c, int d),
	output_relation(int i1, int i2, int a3, int a4, int s, int a, int b),
	dump_params(int a, char *b, char *c),
	find_start_points(),
	siever(int i, int s),
	int_norm(int d, int c, int *n),
	subt(int *a, int *b, int *c),
	dump_int_factors(int a, int b),
	dump_alg_factors(int a, int b),
	enqu(int a, int *b),
	gcd(int a, int b),
	divv(int *a, int *b, int *q, int *r),
	mpcopy(int *a, int *b),
	trial_int_div(int *n, int s, int *r1, int *r2, int l, int c),
	sieve_scanner(int s),
	mpower(int *b, int *n, int *m, int *t),
	mpctoi(char *s, int *n),
	mult(int *a, int *b, int *c),
	mpsqrt(int *n, int *t),
	get_range(int *a, int *b, int *c),
	mark_range_as_used(int a, int b),
	compute_coeffs(int a, int b, int *v1, int *v2),
	even_int_siever(int c, int s),
	even_alg_siever(int c, int s),
	prepare_bounds(int *v1, int *v2, int s),
	sort_int_factor_list(int *list),
	check_lattice(char *str, int i, int*v1, int *v2),
	prime_power_roots(),
	alg_smallp_roots(),
	int_smallp_roots(),
	factor_by_resieve(int s),
	check_start_points(),
	update_global_count(),		// not used
	find_success_locations(),
	vector_int_sieve(int s),
	vector_alg_sieve(int s),
	display_reductions(),
	init_faclist(),
	partition_loops(),
	sqfprep(),
	get_ascii(),
	get_fbase(),
	restore_logs(),
	reset_startpts(),
	init_count_variables(),
	init_time_variables(),
	show_region(int index),
	update_int_start_points(int s),
	update_alg_start_points(int s),
    update_range_file(int q_value, int current_q);
 
int
	trial_alg_div(int *n, int s, int *r, int l, int c),
	get_ascii_fbase(),
	alg_norm(int a, int b, int *n),
	mpcmp(int *a, int *b),
	div_single(int *n, int d, int *r),
	modmultadd_asm(int a, int b, int c, int d),
	single_modinv(int a, int p),
	single_modinv_my(int a, int p),
	ptest(int *n),
	fast_ptest(int *n),
	relprime(int a, int b),
	squfof(int *n, int *r1, int *r2),
	small_qs(int *n, int *r1, int *r2),
	mod_mult_asm(int a, int b, int c),
	pollardrho2(int *n, int *f1, int *f2),
	reduce_lattice(int a, int b, int *v1, int *v2),
    newest_latred_exp(int a, int b, int *v1, int *v2),
	mrint(double d),
	binary_inv(int a, int b),
	int_hash(int c, int d),
	alg_hash(int c, int d),
	show_hash(int b),
	reject_this_q(),
	show_areas(),
	show_int_vector(),
	show_alg_vector();

__int64 prectime();

double get_time();

unsigned char get_alg_log(int b, int a),
			  refine_int_log(int b, int a),
			  refine_alg_log(int b, int a, int s);

__inline int round_double(double a);
__inline void new_divrem_asm(unsigned int a,  unsigned int b, unsigned int c,  unsigned int *d);