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

/*--------------------------------------------------------------------*/
static int compare_relations(DB *db, const DBT *rel1,
				const DBT *rel2)
{
	int64 a1, a2;
	uint32 b1, b2;

	memcpy(&a1, rel1->data, sizeof(int64));
	memcpy(&b1, (uint8 *)rel1->data + sizeof(int64), sizeof(uint32));
	memcpy(&a2, rel2->data, sizeof(int64));
	memcpy(&b2, (uint8 *)rel2->data + sizeof(int64), sizeof(uint32));

	if (b1 > b2)
		return 1;
	if (b1 < b2)
		return -1;

	if (a1 > a2)
		return 1;
	if (a1 < a2)
		return -1;

	return 0;
}

/*--------------------------------------------------------------------*/
static int compress_relation(DB *db, 
			const DBT *prev_key, const DBT *prev_data, 
			const DBT *key, const DBT *data, 
			DBT *dest)
{
	int64 a2;
	uint32 b1, b2;
	uint32 count = 0;
	uint8 a_neg = 0;
	uint8 b_diff_zero = 0;
	uint8 tmp[24];

	memcpy(&b1, (uint8 *)prev_key->data + sizeof(int64), sizeof(uint32));
	memcpy(&a2, key->data, sizeof(int64));
	memcpy(&b2, (uint8 *)key->data + sizeof(int64), sizeof(uint32));

	if (a2 < 0) {
		a_neg = 2;
		a2 = -a2;
	}
	if (b1 == b2)
		b_diff_zero = 1;

	count = compress_p(tmp, 
			a2 << 2 | a_neg | b_diff_zero, 
			count);

	if (!b_diff_zero)
		count = compress_p(tmp, b2 - b1, count);

	/* the BDB docs only show this in an example and not in the
	   API reference, but the compressed buffer is going to have 
	   several key/data pairs concatenated, and the decompression
	   routine has to have each relation know how large its 
	   compressed version is because it cannot blindly use the
	   size of the compressed buffer */

	count = compress_p(tmp, data->size, count);

	dest->size = count + data->size;
	if (dest->ulen < dest->size)
		return DB_BUFFER_SMALL;

	memcpy(dest->data, tmp, count);
	memcpy((uint8 *)dest->data + count, data->data, data->size);
	return 0;
}

/*--------------------------------------------------------------------*/
static int decompress_relation(DB *db, 
		const DBT *prev_key, const DBT *prev_data, 
		DBT *compressed, 
		DBT *key, DBT *data)

{
	uint64 a2;
	uint32 b1, b2;
	uint32 b_diff_zero;
	uint32 a_neg;
	uint32 count = 0;
	uint32 data_size;

	a2 = decompress_p((uint8 *)compressed->data, &count);
	a_neg = a2 & 2;
	b_diff_zero = a2 & 1;
	a2 >>= 2;
	if (a_neg)
		a2 = -a2;

	memcpy(&b2, (uint8 *)prev_key->data + sizeof(int64), sizeof(uint32));
	if (!b_diff_zero)
		b1 += (uint32)decompress_p((uint8 *)compressed->data, &count);

	data_size = (uint32)decompress_p((uint8 *)compressed->data, &count);

	key->size = sizeof(uint64) + sizeof(uint32);
	data->size = data_size;
	if (key->ulen < key->size || data->ulen < data->size)
		return DB_BUFFER_SMALL;

	memcpy(data->data, (uint8 *)compressed->data + count, data_size);
	memcpy(key->data, &a2, sizeof(uint64));
	memcpy((uint8 *)key->data + sizeof(uint64), &b2, sizeof(uint32));

	/* the BDB docs only show this in an example and not in the
	   API reference, but the compressed buffer needs to know
	   how many bytes the current relation needed */

	compressed->size = count + data_size;
	return 0;
}

/*--------------------------------------------------------------------*/
static DB * init_relation_db(msieve_obj *obj, char *name, uint32 open_flags)
{
	DB * reldb = NULL;
	int32 status;

	status = db_create(&reldb, NULL, 0);
	if (status != 0) {
		logprintf(obj, "DB creation failed: %s\n",
				db_strerror(status));
		goto error_finished;
	}

	status = reldb->set_cachesize(reldb, 0, 100*1024*1024, 0);
	if (status != 0) {
		logprintf(obj, "DB set cache failed: %s\n",
				db_strerror(status));
		goto error_finished;
	}

	status = reldb->set_bt_compare(reldb, compare_relations);
	if (status != 0) {
		logprintf(obj, "DB set compare fcn failed: %s\n",
				db_strerror(status));
		goto error_finished;
	}

	status = reldb->set_bt_compress(reldb, compress_relation,
				decompress_relation);
	if (status != 0) {
		logprintf(obj, "DB set compress fcn failed: %s\n",
				db_strerror(status));
		goto error_finished;
	}

	status = reldb->open(reldb, NULL, name, NULL, DB_BTREE, 
				open_flags, 
				(open_flags & DB_RDONLY) ? 0444 : 0664);
	if (status != 0) {
		logprintf(obj, "DB open failed: %s\n",
				db_strerror(status));
		goto error_finished;
	}

	return reldb;

error_finished:
	if (reldb)
		reldb->close(reldb, 0);
	return NULL;

}

/*--------------------------------------------------------------------*/
#define LOG2_BIN_SIZE 17
#define BIN_SIZE (1 << (LOG2_BIN_SIZE))
#define TARGET_HITS_PER_PRIME 40.0

uint32 nfs_purge_duplicates(msieve_obj *obj, factor_base_t *fb,
				uint64 max_relations, 
				uint64 *num_relations_out) {

	uint32 i;
	int32 status;
	savefile_t *savefile = &obj->savefile;
	char buf[LINE_BUF_SIZE];
	uint64 curr_relation;
	uint64 num_relations;
	uint64 num_skipped_b;
	uint32 batch_size;
	mpz_t scratch;

	uint32 *prime_bins;
	double bin_max;

	uint8 tmp_factors[COMPRESSED_P_MAX_SIZE];
	relation_t tmp_rel;

	DB * reldb;
	DBT store_buf;
	void *store_ptr;

	logprintf(obj, "commencing duplicate removal\n");

	tmp_rel.factors = tmp_factors + 2;

	savefile_open(savefile, SAVEFILE_READ);

	sprintf(buf, "%s.db", savefile->name);
	reldb = init_relation_db(obj, buf, DB_CREATE | DB_TRUNCATE);
//	reldb = init_relation_db(obj, buf, DB_RDONLY);
	if (reldb == NULL) {
		printf("error opening relation DB\n");
		exit(-1);
	}
#if 0
	{
		DBC *dbc;
		DBT key, data;

		memset(&key, 0, sizeof(DBT));
		memset(&data, 0, sizeof(DBT));

		status = reldb->cursor(reldb, NULL, &dbc, DB_CURSOR_BULK);
		if (status != 0) {
			logprintf(obj, "DBC open failed: %s\n",
						db_strerror(status));
			exit(-1);
		}

		while ((status = dbc->get(dbc, &key, 
				&data, DB_NEXT)) != DB_NOTFOUND) {

			if (status != 0) {
				logprintf(obj, "DBC get failed: %s\n",
							db_strerror(status));
				exit(-1);
			}

			printf("%u %u\t", key.size, data.size);
			if (key.size == 12) {
				int64 a;
				uint32 b;
				memcpy(&a, key.data, 8);
				memcpy(&b, (uint8 *)key.data + 8, 4);
				printf("%u %" PRId64 , b, a);
			}
			printf("\n");
		}
		dbc->close(dbc);
		reldb->close(reldb, 0);
		exit(0);
	}
#endif
	memset(&store_buf, 0, sizeof(DBT));
	store_buf.data = xcalloc(1, 30*1024*1024);
	store_buf.ulen = 30*1024*1024;
	store_buf.flags = DB_DBT_USERMEM | DB_DBT_BULK;
	DB_MULTIPLE_WRITE_INIT(store_ptr, &store_buf);

	prime_bins = (uint32 *)xcalloc((size_t)1 << (32 - LOG2_BIN_SIZE),
					sizeof(uint32));

	num_relations = 0;
	curr_relation = 0;
	num_skipped_b = 0;
	batch_size = 0;
	mpz_init(scratch);
	savefile_read_line(buf, sizeof(buf), savefile);

	while (!savefile_eof(savefile)) {

		uint32 num_r, num_a;
		uint32 array_size;
		uint8 *rel_key;
		uint8 *rel_data;

		if (buf[0] != '-' && !isdigit(buf[0])) {

			/* no relation on this line */

			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}

		/* read and verify the relation */

		curr_relation++;
		if (max_relations && curr_relation >= max_relations)
			break;

		status = nfs_read_relation(buf, fb, &tmp_rel, 
					&array_size, 0, scratch);
		if (status != 0) {
			if (status == -99)
				num_skipped_b++;
			else
			    logprintf(obj, "error %d reading "
					"relation %" PRIu64 "\n",
					status, num_relations);
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}

		if (tmp_rel.b == 0) {
			/* skip free relations; we'll add them in later */
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}

		num_relations++;
		if (num_relations % 10000000 == 0) {
			printf("read %" PRIu64 "M relations\n", 
					curr_relation / 1000000);
		}

		/* relation is good; queue for DB write */

		num_r = tmp_rel.num_factors_r;
		num_a = tmp_rel.num_factors_a;

		DB_MULTIPLE_KEY_RESERVE_NEXT(store_ptr, &store_buf,
				rel_key, sizeof(uint32)+sizeof(int64),
				rel_data, array_size + 2);

		batch_size++;
		if (rel_key == NULL || rel_data == NULL) {

			/* relation doesn't fit, commit the batch and
			   start the next */

			printf("batch size %u\n", batch_size);
			status = reldb->sort_multiple(reldb, &store_buf,
					NULL, DB_MULTIPLE_KEY);
			if (status != 0) {
				logprintf(obj, "DB bulk sort failed: %s\n",
						db_strerror(status));
				exit(-1);
			}

			status = reldb->put(reldb, NULL, &store_buf,
					NULL, DB_MULTIPLE_KEY);
			if (status != 0) {
				logprintf(obj, "DB bulk write failed: %s\n",
						db_strerror(status));
				exit(-1);
			}

			DB_MULTIPLE_WRITE_INIT(store_ptr, &store_buf);

			DB_MULTIPLE_KEY_RESERVE_NEXT(store_ptr, &store_buf,
					rel_key, sizeof(int64)+sizeof(uint32),
					rel_data, array_size + 2);
			batch_size = 1;
		}

		memcpy(rel_key, &tmp_rel.a, sizeof(int64));
		memcpy(rel_key + sizeof(int64), &tmp_rel.b, sizeof(uint32));
		rel_data[0] = num_r;
		rel_data[1] = num_a;
		memcpy(rel_data + 2, tmp_rel.factors, array_size);

		/* add the factors of tmp_rel to the counts of (32-bit) primes */
	   
		for (i = array_size = 0; i < num_r + num_a; i++) {
			uint64 p = decompress_p(tmp_rel.factors, 
						&array_size);

			if (p >= ((uint64)1 << 32))
				continue;

			prime_bins[p / BIN_SIZE]++;
		}

		/* get the next line */

		savefile_read_line(buf, sizeof(buf), savefile);
	}

	/* commit the final batch */
	printf("batch size %u\n", batch_size);
	if (batch_size > 0) {
		status = reldb->sort_multiple(reldb, &store_buf,
				NULL, DB_MULTIPLE_KEY);
		if (status != 0) {
			logprintf(obj, "DB bulk sort failed: %s\n",
					db_strerror(status));
			exit(-1);
		}

		status = reldb->put(reldb, NULL, &store_buf,
					NULL, DB_MULTIPLE_KEY);
		if (status != 0) {
			logprintf(obj, "DB bulk write failed: %s\n",
							db_strerror(status));
			exit(-1);
		}
	}
	free(store_buf.data);

	reldb->stat_print(reldb, DB_STAT_ALL);

	status = reldb->close(reldb, 0);
	if (status != 0) {
		logprintf(obj, "DB close failed: %s\n",
						db_strerror(status));
		exit(-1);
	}

	sprintf(buf, "%s.db", savefile->name);
	printf("file size %" PRIu64 "\n", get_file_size(buf));

	reldb = init_relation_db(obj, buf, 0);
	if (reldb == NULL) {
		printf("error opening relation DB\n");
		exit(-1);
	}

	printf("compacting DB...\n");

	reldb->compact(reldb, NULL, NULL, NULL, 
			NULL, DB_FREE_SPACE, NULL);
	if (status != 0) {
		logprintf(obj, "DB compact failed: %s\n",
						db_strerror(status));
		exit(-1);
	}

	reldb->stat_print(reldb, DB_STAT_ALL);

	status = reldb->close(reldb, 0);
	if (status != 0) {
		logprintf(obj, "DB close failed: %s\n",
						db_strerror(status));
		exit(-1);
	}

	printf("file size %" PRIu64 "\n", get_file_size(buf));

	mpz_clear(scratch);
	savefile_close(savefile);

	if (num_skipped_b > 0)
		logprintf(obj, "skipped %" PRIu64 
				" relations with b > 2^32\n",
				num_skipped_b);
	logprintf(obj, "found %" PRIu64 " relations\n", num_relations);

	/* the large prime cutoff for the rest of the filtering
	   process should be chosen here. We don't want the bound
	   to depend on an arbitrarily chosen factor base, since
	   that bound may be too large or much too small. The former
	   would make filtering take too long, and the latter 
	   could make filtering impossible.

	   Conceptually, we want the bound to be the point below
	   which large primes appear too often in the dataset. */

	i = 1 << (32 - LOG2_BIN_SIZE);
	bin_max = (double)BIN_SIZE * i /
			(log((double)BIN_SIZE * i) - 1);
	for (i--; i > 2; i--) {
		double bin_min = (double)BIN_SIZE * i /
				(log((double)BIN_SIZE * i) - 1);
		double hits_per_prime = (double)prime_bins[i] /
						(bin_max - bin_min);
		if (hits_per_prime > TARGET_HITS_PER_PRIME)
			break;
		bin_max = bin_min;
	}

	free(prime_bins);
	*num_relations_out = num_relations;
	return BIN_SIZE * (i + 0.5);
}
