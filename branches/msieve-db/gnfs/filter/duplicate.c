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

typedef struct {
	DB *read_db;
	DBC *read_curs;
	void *read_ptr;
	DBT read_key;
	DBT read_data;
	DBT curr_key;
	DBT curr_data;
} tmp_db_t;

#define BULK_BUF_ROUND(x) ((x) & ~(1024 - 1))

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
			size_t batch_size, uint32 num_db, 
			uint64 num_relations_in,
			uint64 *num_relations_out)
{
	uint32 i, j;
	int32 status;
	char buf[LINE_BUF_SIZE];
	uint64 num_relations;
	DB *write_db;
	tmp_db_t *tmp_db;
	tmp_db_t **tmp_db_ptr;
	void *store_ptr;
	DBT store_buf;
	size_t read_batch_size;
	size_t write_batch_size;

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

	/* allocate 60% of the available memory for write
	   buffers, with the remainder split evenly across all the
	   temporary DBs for read buffers */

	write_batch_size = 0.6 * batch_size;
	read_batch_size = (batch_size - write_batch_size) / num_db;
	read_batch_size = MAX(read_batch_size, 1048576);

	/* open the destination DB */

	memset(&store_buf, 0, sizeof(DBT));
	store_buf.ulen = BULK_BUF_ROUND(write_batch_size);
	store_buf.flags = DB_DBT_USERMEM | DB_DBT_BULK;
	store_buf.data = xmalloc(store_buf.ulen);
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

		t->read_data.ulen = BULK_BUF_ROUND(read_batch_size);
		t->read_data.flags = DB_DBT_USERMEM | DB_DBT_BULK;
		t->read_data.data = xmalloc(t->read_data.ulen);
	}

	/* do the initial read from each DB. Squeeze out any
	   DBs that cannot deliver, heapify all the ones that can */

	for (i = j = 0; i < num_db; i++) {

		tmp_db_t *t = tmp_db + i;

		status = t->read_curs->get(t->read_curs, 
					&t->read_key, 
					&t->read_data, 
					DB_NEXT | DB_MULTIPLE_KEY);
		if (status == 0) {
			DB_MULTIPLE_INIT(t->read_ptr, &t->read_data);
			DB_MULTIPLE_KEY_NEXT(t->read_ptr, &t->read_data,
					t->curr_key.data,
					t->curr_key.size,
					t->curr_data.data,
					t->curr_data.size);
			tmp_db_ptr[j++] = t;
		}
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
	DB_MULTIPLE_WRITE_INIT(store_ptr, &store_buf);
	while (i < j) {

		tmp_db_t *t = tmp_db_ptr[i];

		if (memcmp(t->curr_key.data, prev_key, KEY_SIZE) != 0) {
			num_relations++;
			DB_MULTIPLE_KEY_WRITE_NEXT(store_ptr, &store_buf,
					t->curr_key.data, t->curr_key.size,
					t->curr_data.data, t->curr_data.size);

			if (store_ptr == NULL) {

				/* relation doesn't fit, commit the batch and
				   start the next */

				status = write_db->put(write_db, NULL, 
							&store_buf, NULL, 
							DB_MULTIPLE_KEY);
				if (status != 0) {
					logprintf(obj, "DB bulk write 1 "
							"failed: %s\n",
							db_strerror(status));
					exit(-1);
				}

				/* all the data in write_db is written in sorted
				   order, so compact the range we just wrote */

				memcpy(end_buf, t->curr_key.data, KEY_SIZE);
				status = write_db->compact(write_db, NULL,
						&start_key, &end_key, NULL, 
						0, NULL);
				if (status != 0) {
					logprintf(obj, "compact 1 failed: %s\n",
							db_strerror(status));
					exit(-1);
				}
				memcpy(start_buf, t->curr_key.data, KEY_SIZE);

				DB_MULTIPLE_WRITE_INIT(store_ptr, &store_buf);

				DB_MULTIPLE_KEY_WRITE_NEXT(store_ptr, &store_buf,
							t->curr_key.data, 
							t->curr_key.size,
							t->curr_data.data, 
							t->curr_data.size);
			}
			memcpy(prev_key, t->curr_key.data, KEY_SIZE);
		}

		DB_MULTIPLE_KEY_NEXT(t->read_ptr, &t->read_data,
					t->curr_key.data,
					t->curr_key.size,
					t->curr_data.data,
					t->curr_data.size);

		if (t->read_ptr == NULL) {

			status = t->read_curs->get(t->read_curs, 
					&t->read_key, 
					&t->read_data, 
					DB_NEXT | DB_MULTIPLE_KEY);
			if (status != 0) {
				i++;
			}
			else {
				DB_MULTIPLE_INIT(t->read_ptr, 
						&t->read_data);
				DB_MULTIPLE_KEY_NEXT(t->read_ptr, 
						&t->read_data,
						t->curr_key.data,
						t->curr_key.size,
						t->curr_data.data,
						t->curr_data.size);
			}
		}

		if (j - i > 1)
			heapify(tmp_db_ptr + i, 0, j-i);
	}

	/* commit and compact the last batch */

	status = write_db->put(write_db, NULL, &store_buf,
				NULL, DB_MULTIPLE_KEY);
	if (status != 0) {
		logprintf(obj, "DB bulk write 2 failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	status = write_db->compact(write_db, NULL,
			&start_key, NULL, NULL, 
			DB_FREE_SPACE, NULL);
	if (status != 0) {
		logprintf(obj, "DB compact 2 failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	/* close DBs; delete the input ones */

	write_db->close(write_db, 0);

	for (i = 0; i < num_db; i++) {

		tmp_db_t *t = tmp_db + i;

		free(t->read_data.data);

		t->read_curs->close(t->read_curs);
		t->read_db->close(t->read_db, 0);

		sprintf(buf, "tmpdb.%u", i);
		filter_env->dbremove(filter_env, NULL, buf, NULL, 0);
	}
	free(tmp_db);
	free(tmp_db_ptr);
	free(store_buf.data);

	logprintf(obj, "wrote %" PRIu64 " duplicates and %"
			PRIu64 " unique relations\n", 
			num_relations_in - num_relations, 
			num_relations);
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
static uint32 db_fills_cache(DB_ENV *env, char *db_name)
{
	uint32 i = 0;
	uint32 evicted = 0;
	DB_MPOOL_FSTAT **fstats = NULL;

	/* a DB is considered to have filled up the cache when
	   some of its pages have been evicted. This is actually
	   somewhat convoluted to determine, because the cache
	   can still be full of dirty pages from other databases,
	   which are being evicted continuously */

	if (env->memp_stat(env, NULL, &fstats, 0) != 0)
		return 0;

	if (fstats != NULL) {

		while (1) {
			DB_MPOOL_FSTAT *f = fstats[i++];

			if (f == NULL)
				break;

			if (strcmp(f->file_name, db_name) == 0) {
				evicted = f->st_page_out;
				break;
			}
		}

		free(fstats);
	}

	return (evicted > 2000);
}

/*--------------------------------------------------------------------*/
uint32 nfs_purge_duplicates(msieve_obj *obj, factor_base_t *fb,
				uint64 mem_size,
				uint64 max_relations, 
				uint64 *num_relations_out) 
{

	int32 status;
	savefile_t *savefile = &obj->savefile;
	char buf[LINE_BUF_SIZE];
	char db_name[24];
	uint64 curr_relation;
	uint64 num_relations;
	uint64 num_skipped_b;
	mpz_t scratch;

	uint32 bound;
	uint32 num_db = 0;
	size_t batch_size;
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

	/* we divide the input memory bound into an allocation for
	   Berkeley DB, and an allocation for the bulk buffers that
	   read from and write to it. We split it 90-10, with an 
	   upper limit of 3GB for the bulk buffers */

	if (0.1 * mem_size > 0xc0000000)
		batch_size = 0xc0000000;
	else
		batch_size = 0.1 * mem_size;

	cache_size = mem_size - batch_size; 

	/* pass 1: break up the input dataset into temporary
	   databases, each about the size of the DB cache and
	   each in sorted order */

	sprintf(buf, "%s.filter", savefile->name);
	filter_env = init_filter_env(obj, buf, cache_size);
	if (filter_env == NULL) {
		printf("error opening filtering DB env\n");
		exit(-1);
	}

	sprintf(db_name, "tmpdb.%u", num_db++);
	reldb = init_relation_db(obj, filter_env, db_name, DB_CREATE);
	if (reldb == NULL) {
		printf("error opening relation DB\n");
		exit(-1);
	}

	memset(&store_buf, 0, sizeof(DBT));
	store_buf.ulen = BULK_BUF_ROUND(batch_size);
	store_buf.flags = DB_DBT_USERMEM | DB_DBT_BULK;
	store_buf.data = xmalloc(store_buf.ulen);
	DB_MULTIPLE_WRITE_INIT(store_ptr, &store_buf);

	num_relations = 0;
	curr_relation = 0;
	num_skipped_b = 0;
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
				logprintf(obj, "DB bulk sort 3 failed: %s\n",
						db_strerror(status));
				exit(-1);
			}

			status = reldb->put(reldb, NULL, &store_buf,
					NULL, DB_MULTIPLE_KEY);
			if (status != 0) {
				logprintf(obj, "DB bulk write 3 failed: %s\n",
						db_strerror(status));
				exit(-1);
			}

			DB_MULTIPLE_WRITE_INIT(store_ptr, &store_buf);

			DB_MULTIPLE_KEY_RESERVE_NEXT(store_ptr, &store_buf,
					rel_key, KEY_SIZE,
					rel_data, array_size + 2);

			/* if the number of committed relations overflows
			    the cache, flush the DB to disk and start another */

			if (db_fills_cache(filter_env, db_name)) {

				/* the DB is small enough to fit in cache,
				   so compacting it is a cheap operation. 
				   Berkeley DB only fills up about 75% of
				   disk pages, so we save space and improve
				   access times by compacting */ 

				status = reldb->compact(reldb, NULL, NULL, 
							NULL, NULL, 
							DB_FREE_SPACE, NULL);
				if (status != 0) {
					logprintf(obj, "compact 3 failed: %s\n",
							db_strerror(status));
					exit(-1);
				}

				status = reldb->close(reldb, 0);
				if (status != 0) {
					logprintf(obj, "DB close failed: %s\n",
							db_strerror(status));
					exit(-1);
				}

				sprintf(db_name, "tmpdb.%u", num_db++);

				reldb = init_relation_db(obj, filter_env, 
							db_name, DB_CREATE);
				if (reldb == NULL) {
					printf("error opening relation DB\n");
					exit(-1);
				}
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
		logprintf(obj, "DB bulk sort 4 failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	status = reldb->put(reldb, NULL, &store_buf,
				NULL, DB_MULTIPLE_KEY);
	if (status != 0) {
		logprintf(obj, "DB bulk write 4 failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	status = reldb->compact(reldb, NULL, NULL, 
				NULL, NULL, 
				DB_FREE_SPACE, NULL);
	if (status != 0) {
		logprintf(obj, "DB compact 4 failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	status = reldb->close(reldb, 0);
	if (status != 0) {
		logprintf(obj, "DB close 4 failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	/* free pass 1 stuff */

	mpz_clear(scratch);
	savefile_close(savefile);
	free(store_buf.data);

	if (num_skipped_b > 0)
		logprintf(obj, "skipped %" PRIu64 
				" relations with b > 2^32\n",
				num_skipped_b);
	logprintf(obj, "found %" PRIu64 " relations\n", num_relations);

	/* pass 2: merge the temporary DBs into a single output one */

	bound = merge_relation_files(obj, filter_env, 
				batch_size, num_db, 
				num_relations,
				num_relations_out);
	free_filter_env(obj, filter_env);
	return bound;
}

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

