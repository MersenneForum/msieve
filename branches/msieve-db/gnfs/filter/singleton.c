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

#include "filter.h"

#define IDEAL_DELETE 126
#define IDEAL_COUNT_BIAS 127

/*--------------------------------------------------------------------*/
static uint32 pack_ideal(uint8 *array, ideal_t *i)
{
	uint32 count = 0;

	count = compress_p(array, i->p << 1 | i->rat_or_alg, count);
	count = compress_p(array, i->r, count);
	return count;
}

/*--------------------------------------------------------------------*/
static void unpack_ideal(ideal_t *i, uint8 *array)
{
	uint32 count = 0;
	uint64 p = decompress_p(array, &count);
	uint64 r = decompress_p(array, &count);

	i->p = p >> 1;
	i->r = r;
	i->rat_or_alg = p & 1;
}

/*--------------------------------------------------------------------*/
static int compare_ideals(DB *db, const DBT *i1,
			const DBT *i2)
{
	/* Ordering is by prime, then by root of prime,
	   then by rational or algebraic type. This ordering
	   is designed to put the smallest ideals first */

	ideal_t a, b;

	unpack_ideal(&a, (uint8 *)i1->data);
	unpack_ideal(&b, (uint8 *)i2->data);

	if (a.p < b.p)
		return -1;
	if (a.p > b.p)
		return 1;
		
	if (a.r < b.r)
		return -1;
	if (a.r > b.r)
		return 1;

	if (a.rat_or_alg < b.rat_or_alg)
		return -1;
	if (a.rat_or_alg > b.rat_or_alg)
		return 1;
		
	return 0;
}

#if 0
/*--------------------------------------------------------------------*/
static int compress_ideal(DB *db, 
			const DBT *prev_key, const DBT *prev_data, 
			const DBT *key, const DBT *data, 
			DBT *dest)
{
	ideal_t i1, i2;
	uint64 p1, p2;
	uint64 r2;
	uint32 count = 0;
	uint8 p_diff_zero = 0;
	uint8 tmp[24];

	unpack_ideal(&i1, (uint8 *)prev_key->data);
	unpack_ideal(&i2, (uint8 *)key->data);

	p1 = i1.p;
	p2 = i2.p;
	r2 = i2.r;

	if (p1 == p2)
		p_diff_zero = 1;

	count = compress_p(tmp, 
			r2 << 2 | i2.rat_or_alg << 1 | p_diff_zero, 
			count);

	if (!p_diff_zero)
		count = compress_p(tmp, p2 - p1, count);

	count = compress_p(tmp, data->size, count);

	dest->size = count + data->size;
	if (dest->ulen < dest->size)
		return DB_BUFFER_SMALL;

	memcpy(dest->data, tmp, count);
	memcpy((uint8 *)dest->data + count, data->data, data->size);
	return 0;
}

/*--------------------------------------------------------------------*/
static int decompress_ideal(DB *db, 
		const DBT *prev_key, const DBT *prev_data, 
		DBT *compressed, 
		DBT *key, DBT *data)
{
	ideal_t i1, i2;
	uint64 r2;
	uint8 p_diff_zero;
	uint32 count = 0;
	uint32 c_count = 0;
	uint32 data_size;
	uint8 tmp[24];

	r2 = decompress_p((uint8 *)compressed->data, &c_count);
	p_diff_zero = r2 & 1;
	i2.rat_or_alg = (r2 >> 1) & 1;
	i2.r = r2 >> 2;

	unpack_ideal(&i1, (uint8 *)prev_key->data);
	i2.p = i1.p;
	if (!p_diff_zero)
		i2.p += decompress_p((uint8 *)compressed->data, &c_count);

	data_size = (uint32)decompress_p((uint8 *)compressed->data, &c_count);

	count = pack_ideal(tmp, &i2);
	key->size = count;
	data->size = data_size;
	if (key->ulen < key->size || data->ulen < data->size)
		return DB_BUFFER_SMALL;

	memcpy(key->data, tmp, count);
	memcpy(data->data, (uint8 *)compressed->data + c_count, data_size);
	compressed->size = c_count + data_size;
	return 0;
}
#endif

/*--------------------------------------------------------------------*/
static void write_lp_file_final(msieve_obj *obj,
			void *rel_read_stream,
			uint64 mem_size,
			filter_t *filter)
{
	char buf[LINE_BUF_SIZE];
	uint64 num_relations;
	uint64 last_relnum;
	DBT *rel_key;
	DBT *rel_data;
	FILE * out_fp;
	relation_ideal_t tmp_lp;
	uint32 delete_rel;
	uint32 header_words;

	logprintf(obj, "commencing LP file conversion\n");

	/* pass 2: break up the sorted list of large ideals into
	   a sorted list or relations; convert each ideal to an
	   integer and squeeze out relations with singleton ideals */

	sprintf(buf, "%s.lp", obj->savefile.name);
	out_fp = fopen(buf, "wb");
	if (out_fp == NULL) {
		printf("error opening LP file\n");
		exit(-1);
	}

	stream_db_read_init(rel_read_stream, DB_DUP,
				mem_size,
				&rel_key,
				&rel_data);

	num_relations = 0;
	last_relnum = (uint64)(-1);
	delete_rel = 0;
        header_words = (sizeof(relation_ideal_t) -
		 	sizeof(tmp_lp.ideal_list)) / sizeof(uint32);

	while (rel_key != NULL) {

		uint64 curr_relnum = decompress_one_p((uint8 *)rel_key->data);
		uint32 curr_ideal = decompress_one_p((uint8 *)rel_data->data);
		
		if (curr_relnum != last_relnum) {

			if (!delete_rel && last_relnum != (uint64)(-1)) {
				num_relations++;
				fwrite(&tmp_lp, 
					header_words + tmp_lp.ideal_count, 
					sizeof(uint32), out_fp); 
			}

			memset(&tmp_lp, 0, sizeof(tmp_lp));
			tmp_lp.rel_index = curr_relnum;
			last_relnum = curr_relnum;
			delete_rel = 0;
		}

		if (curr_ideal == IDEAL_DELETE)
			delete_rel = 1;
		else if (curr_ideal < IDEAL_COUNT_BIAS)
			tmp_lp.gf2_factors = curr_ideal;
		else
			tmp_lp.ideal_list[tmp_lp.ideal_count++] = 
					curr_ideal - IDEAL_COUNT_BIAS;

		/* get next ideal */

		stream_db_read_next(rel_read_stream, 
				&rel_key, &rel_data);
	}

	if (!delete_rel && last_relnum != (uint64)(-1)) {
		num_relations++;
		fwrite(&tmp_lp, header_words + tmp_lp.ideal_count, 
			sizeof(uint32), out_fp); 
	}

	filter->num_relations = num_relations;
	logprintf(obj, "reduced to %" PRIu64 " relations "
			"and %" PRIu64 " unique ideals\n", 
			filter->num_relations, filter->num_ideals);
	stream_db_read_close(rel_read_stream, 1);
	fclose(out_fp);
}

/*--------------------------------------------------------------------*/
typedef struct {
	uint8 ideal[24];
	uint32 ideal_count;
	uint8 relnum[12];
	uint32 relnum_count;
	uint64 ideal_assigned;
	uint32 is_singleton;
} ideal_list_t;

static void delete_singletons(ideal_list_t *list,
			uint32 list_size,
			uint64 *num_singleton_ideals,
			void *rel_write_stream)
{
	uint32 i;
	uint8 tmp_idealnum[12];
	uint32 tmp_idealnum_count = compress_one_p(tmp_idealnum, 
						IDEAL_DELETE);

	for (i = 0; i < list_size; i++) {
		ideal_list_t *il = list + i;

		if (!il->is_singleton)
			continue;

		(*num_singleton_ideals)++;
		stream_db_write_next(rel_write_stream,
				il->relnum, il->relnum_count,
				tmp_idealnum, tmp_idealnum_count);
	}
}

static void write_singletons(ideal_list_t *list,
			uint32 list_size,
			uint64 *num_ideals,
			uint64 *ideal_assigned,
			void *rel_write_stream)
{
	uint32 i;
	uint8 tmp_idealnum[12];
	uint32 tmp_idealnum_count;

	for (i = 0; i < list_size; i++) {
		ideal_list_t *il = list + i;

		if (!il->is_singleton)
			continue;

		(*num_ideals)++;
		il->ideal_assigned = (*ideal_assigned)++;

		tmp_idealnum_count = compress_one_p(tmp_idealnum, 
						il->ideal_assigned);

		stream_db_write_next(rel_write_stream,
				il->relnum, il->relnum_count,
				tmp_idealnum, tmp_idealnum_count);
	}
}

static uint32 free_rel_allowed(ideal_list_t *list,
			uint32 list_size)
{
	uint32 i;
	ideal_t tmp_ideal;

	for (i = 0; i < list_size; i++) {
		ideal_list_t *il = list + i;

		unpack_ideal(&tmp_ideal, il->ideal);

		if (tmp_ideal.rat_or_alg == ALGEBRAIC_IDEAL &&
		    2 * tmp_ideal.p + 1 == tmp_ideal.r) {
			return 0;
		}
	}

	return 1;
}

static void write_free_rel(ideal_list_t *list,
			uint32 list_size,
			uint64 free_rel_num,
			void *rel_write_stream)
{
	uint32 i;
	uint8 tmp_idealnum[12];
	uint32 tmp_idealnum_count;
	uint8 tmp_relnum[12];
	uint32 tmp_relnum_count = compress_one_p(tmp_relnum,
				0x8000000000 + free_rel_num);

	for (i = 0; i < list_size; i++) {
		ideal_list_t *il = list + i;

		tmp_idealnum_count = compress_one_p(tmp_idealnum, 
						il->ideal_assigned);

		stream_db_write_next(rel_write_stream,
				tmp_relnum, tmp_relnum_count,
				tmp_idealnum, tmp_idealnum_count);
	}
}

static void write_lp_file_pass2(msieve_obj *obj,
			void *ideal_read_stream,
			uint64 mem_size,
			uint32 free_rel_target,
			filter_t *filter)
{
	savefile_t *savefile = &obj->savefile;
	char buf[LINE_BUF_SIZE];
	uint64 num_ideals;
	uint64 num_singleton_ideals;
	uint64 num_free_rels;
	uint64 ideal_assigned;
	DBT *ideal_key;
	DBT *ideal_data;
	ideal_list_t list[MAX_POLY_DEGREE + 1];
	ideal_t tmp_ideal;
	uint64 last_p;
	uint8 tmp_idealnum[12];
	uint32 tmp_idealnum_count;
	void *rel_write_stream;
	uint32 i;

	logprintf(obj, "commencing singleton removal, pass 2\n");

	/* pass 2: break up the sorted list of large ideals into
	   a sorted list or relations; convert each ideal to an
	   integer and squeeze out relations with singleton ideals */

	sprintf(buf, "%s.filter/tmplp", savefile->name);

	rel_write_stream = stream_db_init(obj, NULL, 0,
				compare_relnums,
				NULL,
				NULL,
				buf);
	if (rel_write_stream == NULL) {
		logprintf(obj, "error opening lp write stream\n");
		exit(-1);
	}

	stream_db_write_init(rel_write_stream, DB_DUP,
				0.9 * mem_size);

	stream_db_read_init(ideal_read_stream, DB_DUP,
				0.1 * mem_size,
				&ideal_key,
				&ideal_data);

	num_ideals = 0;
	num_singleton_ideals = 0;
	num_free_rels = 0;
	ideal_assigned = IDEAL_COUNT_BIAS;
	i = (uint32)(-1);
	last_p = (uint64)(-1);

	/* read all the ideals in ascending order. The first big
	   batch of them will be the list of relations, each
	   containing the number of small factors (i.e. factors below
	   the filtering bound). Then we transition to the large
	   ideals */

	while (ideal_key != NULL) {

		ideal_list_t *il = list + i;

		unpack_ideal(&tmp_ideal, (uint8 *)ideal_key->data);

		if (tmp_ideal.p < 1000) {

			/* current ideal is not an ideal at all, but the
			   count of small factors for the relation */

			tmp_idealnum_count = compress_one_p(tmp_idealnum, tmp_ideal.p);

			stream_db_write_next(rel_write_stream,
					ideal_data->data, ideal_data->size,
					tmp_idealnum, tmp_idealnum_count);
		}
		else if (last_p != (uint64)(-1) &&
			 ideal_key->size == il->ideal_count &&
			 memcmp(ideal_key->data, il->ideal, il->ideal_count) == 0) {

			/* current ideal equals the last one, write it out */

			/* if the previous occurrence of this ideal was 
			   listed as a singleton, then we have two ideals
			   to write out */

			if (il->is_singleton) {

				num_ideals++;
				il->ideal_assigned = ideal_assigned++;
				il->is_singleton = 0;

				tmp_idealnum_count = compress_one_p(tmp_idealnum, 
							il->ideal_assigned);

				stream_db_write_next(rel_write_stream,
						il->relnum, il->relnum_count,
						tmp_idealnum, tmp_idealnum_count);
			}

			stream_db_write_next(rel_write_stream,
					ideal_data->data, ideal_data->size,
					tmp_idealnum, tmp_idealnum_count);
		}
		else {
			/* first time ideal has been seen */

			i++;
			if (tmp_ideal.p != last_p) {
				/* ideal is modulo a different prime; throw
				   away all the singleton ideals seen so far,
				   unless the number of ideals modulo last_p is
				   large enough so that we can issue a free 
				   relation. In that case everything can be saved */

				if (i == free_rel_target &&
				    free_rel_allowed(list, i)) {

					write_singletons(list, i, 
							&num_ideals, 
							&ideal_assigned,
							rel_write_stream);
					write_free_rel(list, i, 
							num_free_rels++,
							rel_write_stream);
				}
				else {
					delete_singletons(list, i, 
							&num_singleton_ideals, 
							rel_write_stream);
				}
				last_p = tmp_ideal.p;
				i = 0;
			}

			il = list + i;
			il->is_singleton = 1;
			il->ideal_count = ideal_key->size;
			il->relnum_count = ideal_data->size;
			memcpy(il->ideal, ideal_key->data, ideal_key->size);
			memcpy(il->relnum, ideal_data->data, ideal_data->size);
		}


		/* get next ideal */

		stream_db_read_next(ideal_read_stream, 
				&ideal_key, &ideal_data);
	}

	/* handle the final list of ideals */
	if (++i == free_rel_target &&
	    free_rel_allowed(list, i)) {

		write_singletons(list, i, 
				&num_ideals, 
				&ideal_assigned,
				rel_write_stream);
		write_free_rel(list, i, 
				num_free_rels++,
				rel_write_stream);
	}
	else {
		delete_singletons(list, i, 
				&num_singleton_ideals, 
				rel_write_stream);
	}

	/* free pass 2 stuff; delete the list of input ideals */

	logprintf(obj, "found %" PRIu64 " ideals (plus "
			"%" PRIu64 " singletons)\n",
			num_ideals, num_singleton_ideals);
	logprintf(obj, "found %" PRIu64 
			" free relations\n", num_free_rels);

	stream_db_write_close(rel_write_stream);
	stream_db_read_close(ideal_read_stream, 1);
	filter->num_ideals = num_ideals;

	/* final pass: convert the distributed list of relation
	   files into the conventional LP disk file */

	write_lp_file_final(obj, rel_write_stream,
			mem_size, filter);

	stream_db_free(rel_write_stream);
}

/*--------------------------------------------------------------------*/
void nfs_write_lp_file(msieve_obj *obj, factor_base_t *fb,
			uint64 mem_size,
			uint32 num_line_num_files,
			filter_t *filter)
{
	savefile_t *savefile = &obj->savefile;
	char buf[LINE_BUF_SIZE];
	uint64 curr_relation;
	uint64 num_relations;
	uint64 next_relation;
	uint64 num_skipped_b;
	void *ideal_write_stream;
	void *read_stream;
	DBT *relnum_key;
	DBT *relnum_data;
	int32 status;
	mpz_t scratch;

	uint8 tmp_factors[COMPRESSED_P_MAX_SIZE];
	uint32 array_size;
	relation_t tmp_rel;

	tmp_rel.factors = tmp_factors;

	logprintf(obj, "commencing singleton removal, pass 1\n");

	savefile_open(savefile, SAVEFILE_READ);

	/* pass 1: break up the input dataset into a list
	   of relations and a list of large ideals
	
	   Note that we do not sort duplicates, and do not
	   compress <ideal,relnum> pairs. This lets BDB
	   store the key only once for all the duplicates,
	   and takes up even less space than trying to
	   compress sorted neighbors */

	sprintf(buf, "%s.filter/tmpideal", savefile->name);

	ideal_write_stream = stream_db_init(obj, NULL, 0,
				compare_ideals,
				NULL,
				NULL,
				buf);
	if (ideal_write_stream == NULL) {
		logprintf(obj, "error opening ideal write stream\n");
		exit(-1);
	}

	sprintf(buf, "%s.filter/tmpi", savefile->name);

	read_stream = stream_db_init(obj, NULL, 
				num_line_num_files,
				compare_relnums,
				compress_relnum,
				decompress_relnum,
				buf);
	if (read_stream == NULL) {
		logprintf(obj, "error opening relnum read stream\n");
		exit(-1);
	}

	stream_db_write_init(ideal_write_stream, DB_DUP,
				0.9 * mem_size);

	stream_db_read_init(read_stream, 0,
				0.1 * mem_size,
				&relnum_key,
				&relnum_data);
	curr_relation = 0;
	num_relations = 0;
	next_relation = 0;
	num_skipped_b = 0;
	mpz_init(scratch);
	if (relnum_key != NULL)
		next_relation = decompress_one_p((uint8 *)relnum_key->data);
	savefile_read_line(buf, sizeof(buf), savefile);

	while (!savefile_eof(savefile)) {

		relation_lp_t tmp_ideal;
		ideal_t fake_ideal;
		uint8 packed_num[24] = {0};
		uint8 packed_ideal[24];
		uint32 packed_num_count;
		uint32 packed_ideal_count;
		uint32 i;

		if (buf[0] != '-' && !isdigit(buf[0])) {
			/* no relation on this line */
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}

		curr_relation++;
		if (curr_relation % 10000000 == 0) {
			printf("read %" PRIu64 "M relations\n", 
					curr_relation / 1000000);
		}

		if (curr_relation < next_relation) {
			/* relation needs to be skipped */
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}

		/* parse relation */

		status = nfs_read_relation(buf, fb, &tmp_rel, 
					&array_size, 0, scratch);
		if (status != 0) {

			if (status == -99)
				num_skipped_b++;
			else
			    logprintf(obj, "error %d reading relation %" PRIu64 "\n",
					status, curr_relation);
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}

		/* relation is good; get its ideal decomposition */

		num_relations++;
		find_large_ideals(&tmp_rel, &tmp_ideal, 
					filter->filtmin_r,
					filter->filtmin_a);

		/* issue one DB write for each ideal; to avoid forgetting
		   about relations that have no large ideals, also create
		   a 'fake' ideal to hold the count of small factors */

		fake_ideal.p = tmp_ideal.gf2_factors;
		fake_ideal.rat_or_alg = RATIONAL_IDEAL;
		fake_ideal.r = 0;

		packed_num_count = compress_one_p(packed_num, num_relations);
		packed_ideal_count = pack_ideal(packed_ideal, &fake_ideal);
		stream_db_write_next(ideal_write_stream,
					packed_ideal,
					packed_ideal_count,
					packed_num, 
					packed_num_count);

		for (i = 0; i < tmp_ideal.ideal_count; i++) {
			packed_ideal_count = pack_ideal(packed_ideal,
						tmp_ideal.ideal_list + i);

			stream_db_write_next(ideal_write_stream,
					packed_ideal,
					packed_ideal_count,
					packed_num, 
					packed_num_count);
		}

		/* get next relation number */

		stream_db_read_next(read_stream, 
				&relnum_key, &relnum_data);
		if (relnum_key == NULL)
			break;

		next_relation = decompress_one_p((uint8 *)relnum_key->data);

		/* get the next line */

		savefile_read_line(buf, sizeof(buf), savefile);
	}

	/* free pass 1 stuff; delete the list of input
	   relation numbers */

	mpz_clear(scratch);
	savefile_close(savefile);
	stream_db_write_close(ideal_write_stream);
	stream_db_read_close(read_stream, 1);

	if (num_skipped_b > 0)
		logprintf(obj, "skipped %" PRIu64 " relations with b > 2^32\n",
				num_skipped_b);
	logprintf(obj, "found %" PRIu64 " relations\n", num_relations);
	filter->num_relations = num_relations;

	/* pass 2: merge the temporary DBs into a single output one */

	write_lp_file_pass2(obj, ideal_write_stream,
			mem_size, fb->afb.poly.degree + 1, 
			filter);

	stream_db_free(ideal_write_stream);
	stream_db_free(read_stream);
}
