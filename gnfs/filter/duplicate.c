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

#define LOG2_BIN_SIZE 17
#define BIN_SIZE (1 << (LOG2_BIN_SIZE))
#define TARGET_HITS_PER_PRIME 40.0

#define KEY_SIZE (sizeof(int64) + sizeof(uint32))

typedef struct {
	DB *read_db;
	DBC *read_curs;
	DBT curr_key;
	DBT curr_data;
} tmp_db_t;

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
		b2 += (uint32)decompress_p((uint8 *)compressed->data, &count);

	data_size = (uint32)decompress_p((uint8 *)compressed->data, &count);

	key->size = KEY_SIZE;
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
static DB_ENV * init_filter_env(msieve_obj *obj, char *dirname,
				uint64 cache_size)
{
	DB_ENV * env = NULL;
	int32 status;

	status = db_env_create(&env, 0);
	if (status != 0) {
		logprintf(obj, "DB_ENV creation failed: %s\n",
				db_strerror(status));
		goto error_finished;
	}

	status = env->set_cachesize(env, 
				(uint32)(cache_size >> 32), 
				(uint32)cache_size, 0);
	if (status != 0) {
		logprintf(obj, "DB_ENV set cache failed: %s\n",
				db_strerror(status));
		goto error_finished;
	}

	status = env->open(env, dirname, 
			DB_CREATE | DB_INIT_CDB | 
			DB_INIT_MPOOL | DB_PRIVATE, 0);
	if (status != 0) {
		logprintf(obj, "DB_ENV open failed: %s\n",
				db_strerror(status));
		goto error_finished;
	}

	return env;

error_finished:
	if (env)
		env->close(env, 0);
	return NULL;
}

/*--------------------------------------------------------------------*/
static void free_filter_env(msieve_obj *obj, DB_ENV *env)
{
	int32 status;

	status = env->memp_sync(env, NULL); 
	if (status != 0) {
		logprintf(obj, "DB_ENV sync failed: %s\n",
				db_strerror(status));
	}

	status = env->close(env, 0);
	if (status != 0) {
		logprintf(obj, "DB_ENV close failed: %s\n",
				db_strerror(status));
	}
}

/*--------------------------------------------------------------------*/
static DB * init_relation_db(msieve_obj *obj, 
				DB_ENV *filter_env, char *name,
				uint32 open_flags)
{
	DB * reldb = NULL;
	int32 status;

	status = db_create(&reldb, filter_env, 0);
	if (status != 0) {
		logprintf(obj, "DB creation failed: %s\n",
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
/* boilerplate code for managing a heap of relations */

#define HEAP_SWAP(a,b) { tmp = a; a = b; b = tmp; }
#define HEAP_PARENT(i)  (((i)-1) >> 1)
#define HEAP_LEFT(i)    (2 * (i) + 1)
#define HEAP_RIGHT(i)   (2 * (i) + 2)

static void heapify(tmp_db_t **h, uint32 index, uint32 size) {

	uint32 c;
	tmp_db_t * tmp;
	for (c = HEAP_LEFT(index); c < (size-1); 
			index = c, c = HEAP_LEFT(index)) {

		if (compare_relations(NULL, &h[c]->curr_key,
				&h[c+1]->curr_key) >= 0)
			c++;

		if (compare_relations(NULL, &h[index]->curr_key,
				&h[c]->curr_key) >= 0) {
			HEAP_SWAP(h[index], h[c]);
		}
		else
			return;
	}
	if (c == (size-1) && 
	    compare_relations(NULL, &h[index]->curr_key,
				&h[c]->curr_key) >= 0) {
		HEAP_SWAP(h[index], h[c]);
	}
}

static void make_heap(tmp_db_t **h, uint32 size) {

	uint32 i;
	for (i = HEAP_PARENT(size); i; i--)
		heapify(h, i-1, size);
}

/*--------------------------------------------------------------------*/
static uint32
merge_relation_files(msieve_obj *obj, DB_ENV *filter_env,
			DBT *store_buf, uint32 num_db, 
			uint64 *num_relations_out)
{
	uint32 i, j;
	int32 status;
	char buf[LINE_BUF_SIZE];
	uint64 num_relations;
	uint64 num_duplicates;
	DB *write_db;
	tmp_db_t *tmp_db;
	tmp_db_t **tmp_db_ptr;
	void *store_ptr;

	DBT start_key, end_key;
	uint8 prev_key[KEY_SIZE] = {0};
	uint8 start_buf[KEY_SIZE] = {0};
	uint8 end_buf[KEY_SIZE];

	uint32 *prime_bins;
	double bin_max;

	logprintf(obj, "commencing duplicate removal, pass 2\n");
	if (num_db > 1)
		logprintf(obj, "merging %u temporary DB files\n", num_db);

	tmp_db = (tmp_db_t *)xcalloc(num_db, sizeof(tmp_db_t));
	tmp_db_ptr = (tmp_db_t **)xcalloc(num_db, sizeof(tmp_db_t *));
	prime_bins = (uint32 *)xcalloc((size_t)1 << (32 - LOG2_BIN_SIZE),
					sizeof(uint32));

	/* open the destination DB */
	write_db = init_relation_db(obj, filter_env, "relations.db", DB_CREATE);
	if (write_db == NULL) {
		printf("error opening output relation DB\n");
		exit(-1);
	}

	/* open all the input DBs */

	for (i = 0; i < num_db; i++) {

		tmp_db_t *t = tmp_db + i;

		sprintf(buf, "tmpdb.%u", i);
		t->read_db = init_relation_db(obj, filter_env, 
					buf, DB_RDONLY);
		if (t->read_db == NULL) {
			printf("error opening input DB\n");
			exit(-1);
		}

		status = t->read_db->cursor(t->read_db, NULL, 
					&t->read_curs, DB_CURSOR_BULK);
		if (status != 0) {
			logprintf(obj, "read cursor open failed: %s\n",
						db_strerror(status));
			exit(-1);
		}
	}

	/* do the initial read from each DB. Squeeze out any
	   DBs that cannot deliver, heapify all the ones that can */

	for (i = j = 0; i < num_db; i++) {

		tmp_db_t *t = tmp_db + i;

		status = t->read_curs->get(t->read_curs, 
					&t->curr_key, 
					&t->curr_data, DB_NEXT);
		if (status == 0)
			tmp_db_ptr[j++] = t;
	}

	if (j > 1)
		make_heap(tmp_db_ptr, j);

	/* we have the relation from each DB that has the 
	   'smallest' key value; progressively read the smallest key 
	   from the list, make sure it isn't a duplicate, write it 
	   out and read the next relation from its DB, heapifying 
	   the new collection */

	memset(&start_key, 0, sizeof(DBT));
	start_key.data = start_buf;
	start_key.size = KEY_SIZE;
	start_key.ulen = KEY_SIZE;
	start_key.flags = DB_DBT_USERMEM;

	memset(&end_key, 0, sizeof(DBT));
	end_key.data = end_buf;
	end_key.size = KEY_SIZE;
	end_key.ulen = KEY_SIZE;
	end_key.flags = DB_DBT_USERMEM;

	i = 0;
	num_relations = 0;
	num_duplicates = 0;
	DB_MULTIPLE_WRITE_INIT(store_ptr, store_buf);
	while (i < j) {

		tmp_db_t *t = tmp_db_ptr[i];

		if (memcmp(t->curr_key.data, prev_key, KEY_SIZE) == 0) {
			num_duplicates++;
		}
		else {
			num_relations++;
			DB_MULTIPLE_KEY_WRITE_NEXT(store_ptr, store_buf,
					t->curr_key.data, t->curr_key.size,
					t->curr_data.data, t->curr_data.size);

			if (store_ptr == NULL) {

				/* relation doesn't fit, commit the batch and
				   start the next */

				status = write_db->put(write_db, NULL, store_buf,
						NULL, DB_MULTIPLE_KEY);
				if (status != 0) {
					logprintf(obj, "DB bulk write failed: %s\n",
							db_strerror(status));
					exit(-1);
				}

				/* all the data in write_db is written in sorted
				   order, so compact the range we just wrote */

				memcpy(end_buf, t->curr_key.data, KEY_SIZE);
				status = write_db->compact(write_db, NULL,
						&start_key, &end_key, NULL, 0, NULL);
				if (status != 0) {
					logprintf(obj, "DB compact failed: %s\n",
							db_strerror(status));
					exit(-1);
				}
				memcpy(start_buf, t->curr_key.data, KEY_SIZE);

				DB_MULTIPLE_WRITE_INIT(store_ptr, store_buf);

				DB_MULTIPLE_KEY_WRITE_NEXT(store_ptr, store_buf,
						t->curr_key.data, t->curr_key.size,
						t->curr_data.data, t->curr_data.size);
			}
			memcpy(prev_key, t->curr_key.data, KEY_SIZE);
		}

		status = t->read_curs->get(t->read_curs, 
					&t->curr_key, 
					&t->curr_data, DB_NEXT);
		if (status != 0)
			i++;
		else
			heapify(tmp_db_ptr + i, 0, j-i);
	}

	/* commit and compact the last batch; we make
	   an ending key of 'all bits set' */

	status = write_db->put(write_db, NULL, store_buf,
				NULL, DB_MULTIPLE_KEY);
	if (status != 0) {
		logprintf(obj, "DB bulk write failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	memset(end_buf, 0xff, KEY_SIZE);
	status = write_db->compact(write_db, NULL,
			&start_key, &end_key, NULL, 0, NULL);
	if (status != 0) {
		logprintf(obj, "DB compact failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	/* close DBs; delete the input ones */

	write_db->close(write_db, 0);

	for (i = 0; i < num_db; i++) {

		tmp_db_t *t = tmp_db + i;

		t->read_curs->close(t->read_curs);
		t->read_db->close(t->read_db, 0);

		sprintf(buf, "tmpdb.%u", i);
		filter_env->dbremove(filter_env, NULL, buf, NULL, 0);
	}
	free(tmp_db);
	free(tmp_db_ptr);

	logprintf(obj, "found %" PRIu64 " duplicates and %"
			PRIu64 " unique relations\n", 
			num_duplicates, num_relations);
#if 0
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
#endif
	free(prime_bins);

	*num_relations_out = num_relations;
	return BIN_SIZE * (i + 0.5);
}

/*--------------------------------------------------------------------*/
uint32 nfs_purge_duplicates(msieve_obj *obj, factor_base_t *fb,
				uint64 max_relations, 
				uint64 *num_relations_out) {

	int32 status;
	savefile_t *savefile = &obj->savefile;
	char buf[LINE_BUF_SIZE];
	uint64 curr_relation;
	uint64 num_relations;
	uint64 num_skipped_b;
	mpz_t scratch;

	uint32 bound;
	uint32 num_db;
	uint32 batch_size;
	uint64 curr_batch_size;
	uint64 cache_size;

	uint8 tmp_factors[COMPRESSED_P_MAX_SIZE];
	relation_t tmp_rel;

	DB * reldb;
	DB_ENV * filter_env;
	DBT store_buf;
	void *store_ptr;

	logprintf(obj, "commencing duplicate removal, pass 1\n");

	tmp_rel.factors = tmp_factors + 2;

	savefile_open(savefile, SAVEFILE_READ);

	/* pass 1: break up the input dataset into temporary
	   databases, each about the size of the DB cache and
	   each in sorted order */

	num_db = 0;
	cache_size = (uint64)100*1024*1024;
	batch_size = 30*1024*1024;

	sprintf(buf, "%s.filter", savefile->name);
	filter_env = init_filter_env(obj, buf, cache_size);
	if (filter_env == NULL) {
		printf("error opening filtering DB env\n");
		exit(-1);
	}

	sprintf(buf, "tmpdb.%u", num_db++);
	reldb = init_relation_db(obj, filter_env, buf, DB_CREATE);
	if (reldb == NULL) {
		printf("error opening relation DB\n");
		exit(-1);
	}

	memset(&store_buf, 0, sizeof(DBT));
	store_buf.data = xmalloc(batch_size);
	store_buf.ulen = batch_size;
	store_buf.flags = DB_DBT_USERMEM | DB_DBT_BULK;
	DB_MULTIPLE_WRITE_INIT(store_ptr, &store_buf);

	num_relations = 0;
	curr_relation = 0;
	num_skipped_b = 0;
	curr_batch_size = 0;
	mpz_init(scratch);
	savefile_read_line(buf, sizeof(buf), savefile);

	while (!savefile_eof(savefile)) {

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

		DB_MULTIPLE_KEY_RESERVE_NEXT(store_ptr, &store_buf,
				rel_key, KEY_SIZE,
				rel_data, array_size + 2);

		if (rel_key == NULL || rel_data == NULL) {

			/* relation doesn't fit, commit the batch and
			   start the next. Note that bulk writes
			   to a compressed DB will not work unless
			   the batch is pre-sorted */

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
					rel_key, KEY_SIZE,
					rel_data, array_size + 2);

			/* if the number of committed relations threatens
			   to overflow the cache, flush the DB to disk
			   and start another */

			curr_batch_size += batch_size;
			if (curr_batch_size > cache_size) {

				/* the DB is small enough to fit in cache,
				   so compacting it is a cheap operation. 
				   Berkeley DB only fills up about 75% of
				   disk pages, so we save space and improve
				   access times by compacting */ 
				   
				status = reldb->compact(reldb, NULL, NULL, 
							NULL, NULL, 
							DB_FREE_SPACE, NULL);
				if (status != 0) {
					logprintf(obj, "DB compact failed: %s\n",
							db_strerror(status));
					exit(-1);
				}

				status = reldb->close(reldb, 0);
				if (status != 0) {
					logprintf(obj, "DB close failed: %s\n",
							db_strerror(status));
					exit(-1);
				}

				sprintf(buf, "tmpdb.%u", num_db++);

				reldb = init_relation_db(obj, filter_env, 
							buf, DB_CREATE);
				if (reldb == NULL) {
					printf("error opening relation DB\n");
					exit(-1);
				}
				curr_batch_size = 0;
			}
		}

		memcpy(rel_key, &tmp_rel.a, sizeof(int64));
		memcpy(rel_key + sizeof(int64), &tmp_rel.b, sizeof(uint32));
		rel_data[0] = tmp_rel.num_factors_r;
		rel_data[1] = tmp_rel.num_factors_a;
		memcpy(rel_data + 2, tmp_rel.factors, array_size);

		/* get the next line */

		savefile_read_line(buf, sizeof(buf), savefile);
	}

	/* commit the final batch */

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

	status = reldb->compact(reldb, NULL, NULL, 
				NULL, NULL, 
				DB_FREE_SPACE, NULL);
	if (status != 0) {
		logprintf(obj, "DB compact failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	status = reldb->close(reldb, 0);
	if (status != 0) {
		logprintf(obj, "DB close failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	/* free pass 1 stuff */

	mpz_clear(scratch);
	savefile_close(savefile);

	if (num_skipped_b > 0)
		logprintf(obj, "skipped %" PRIu64 
				" relations with b > 2^32\n",
				num_skipped_b);
	logprintf(obj, "found %" PRIu64 " relations\n", num_relations);

	/* pass 2: merge the temporary DBs into a single output one */

	bound = merge_relation_files(obj, filter_env, 
				&store_buf, num_db, 
				num_relations_out);
	free(store_buf.data);
	free_filter_env(obj, filter_env);
	return bound;
}

#if 0
	reldb->stat_print(reldb, DB_STAT_ALL);
#endif

#if 0
		/* add the factors of tmp_rel to the counts of (32-bit) primes */
	   
		for (i = array_size = 0; i < num_r + num_a; i++) {
			uint64 p = decompress_p(tmp_rel.factors, 
						&array_size);

			if (p >= ((uint64)1 << 32))
				continue;

			prime_bins[p / BIN_SIZE]++;
		}
#endif

