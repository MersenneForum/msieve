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
/* boilerplate code for managing a heap of ideals */

#define HEAP_SWAP(a,b) { tmp = a; a = b; b = tmp; }
#define HEAP_PARENT(i)  (((i)-1) >> 1)
#define HEAP_LEFT(i)    (2 * (i) + 1)
#define HEAP_RIGHT(i)   (2 * (i) + 2)

static void heapify(tmp_db_t **h, uint32 index, uint32 size) {

	uint32 c;
	tmp_db_t * tmp;
	for (c = HEAP_LEFT(index); c < (size-1); 
			index = c, c = HEAP_LEFT(index)) {

		if (compare_ideals(NULL, &h[c]->curr_key,
				&h[c+1]->curr_key) >= 0)
			c++;

		if (compare_ideals(NULL, &h[index]->curr_key,
				&h[c]->curr_key) >= 0) {
			HEAP_SWAP(h[index], h[c]);
		}
		else
			return;
	}
	if (c == (size-1) && 
	    compare_ideals(NULL, &h[index]->curr_key,
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
static void
merge_ideal_files(msieve_obj *obj, size_t mem_size, 
		uint32 num_db, uint64 *num_ideals_out)
{
	uint32 i, j;
	int32 status;
	savefile_t *savefile = &obj->savefile;
	char buf[LINE_BUF_SIZE];
	uint64 num_ideals;
	DB *write_db;
	DB_ENV * filter_env;
	tmp_db_t *tmp_db;
	tmp_db_t **tmp_db_ptr;
	void *store_ptr;
	DBT store_buf;
	size_t cache_size;
	size_t batch_size;
	size_t read_batch_size;
	size_t write_batch_size;

	DBT start_key, end_key;
	ideal_t prev_ideal;
	ideal_t start_ideal;
	ideal_t end_ideal;

	logprintf(obj, "commencing ideal list build, pass 2\n");
	if (num_db > 1)
		logprintf(obj, "merging %u temporary DB files\n", num_db);

	tmp_db = (tmp_db_t *)xcalloc(num_db, sizeof(tmp_db_t));
	tmp_db_ptr = (tmp_db_t **)xcalloc(num_db, sizeof(tmp_db_t *));

	/* we divide the input memory bound into an allocation for
	   Berkeley DB, and an allocation for the bulk buffers that
	   read from and write to it. We split it 50-50, with an 
	   upper limit of 3GB for the bulk buffers */

	if (0.5 * mem_size > 0xc0000000)
		batch_size = 0xc0000000;
	else
		batch_size = 0.5 * mem_size;

	cache_size = mem_size - batch_size; 

	sprintf(buf, "%s.filter", savefile->name);
	filter_env = init_filter_env(obj, buf, cache_size);
	if (filter_env == NULL) {
		printf("error opening filtering DB env\n");
		exit(-1);
	}

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
	write_db = init_ideal_db(obj, filter_env, "ideals.db", DB_CREATE);

	if (write_db == NULL) {
		printf("error opening output ideal DB\n");
		exit(-1);
	}

	/* open the input DBs for reading; do *not* make them part
	   of the DB environment, as we will be streaming ideals
	   and don't want them evicting cache pages that are being
	   written with ideals */

	for (i = 0; i < num_db; i++) {

		tmp_db_t *t = tmp_db + i;

		sprintf(buf, "%s.filter/tmpdb.%u", savefile->name, i);
		t->read_db = init_ideal_db(obj, NULL, buf, DB_RDONLY);
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

	/* we have the ideal from each DB that has the 
	   'smallest' key value; progressively read the smallest key 
	   from the list, make sure it isn't a duplicate, write it 
	   out and read the next relation from its DB, heapifying 
	   the new collection */

	memset(&prev_ideal, 0, sizeof(ideal_t));

	memset(&start_ideal, 0, sizeof(ideal_t));
	memset(&start_key, 0, sizeof(DBT));
	start_key.data = &start_ideal;
	start_key.size = sizeof(ideal_t);
	start_key.ulen = sizeof(ideal_t);
	start_key.flags = DB_DBT_USERMEM;

	memset(&end_ideal, 0, sizeof(ideal_t));
	memset(&end_key, 0, sizeof(DBT));
	end_key.data = &end_ideal;
	end_key.size = sizeof(ideal_t);
	end_key.ulen = sizeof(ideal_t);
	end_key.flags = DB_DBT_USERMEM;

	i = 0;
	num_ideals = 0;
	DB_MULTIPLE_WRITE_INIT(store_ptr, &store_buf);
	while (i < j) {

		tmp_db_t *t = tmp_db_ptr[i];

		if (memcmp(t->curr_key.data, &prev_ideal, sizeof(ideal_t)) != 0) {

			num_ideals++;
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

				memcpy(&end_ideal, t->curr_key.data, sizeof(ideal_t));
				status = write_db->compact(write_db, NULL,
						&start_key, &end_key, NULL, 
						0, NULL);
				if (status != 0) {
					logprintf(obj, "compact 1 failed: %s\n",
							db_strerror(status));
					exit(-1);
				}
				start_ideal = end_ideal;

				DB_MULTIPLE_WRITE_INIT(store_ptr, &store_buf);

				DB_MULTIPLE_KEY_WRITE_NEXT(store_ptr, &store_buf,
							t->curr_key.data, 
							t->curr_key.size,
							t->curr_data.data, 
							t->curr_data.size);
			}
			memcpy(&prev_ideal, t->curr_key.data, sizeof(ideal_t));
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

		sprintf(buf, "%s.filter/tmpdb.%u", savefile->name, i);
		remove(buf);
	}
	free(tmp_db);
	free(tmp_db_ptr);
	free(store_buf.data);
	free_filter_env(obj, filter_env);

	logprintf(obj, "found %" PRIu64 " unique ideals\n", num_ideals);
	*num_ideals_out = num_ideals;
}

/*--------------------------------------------------------------------*/
void nfs_write_lp_file(msieve_obj *obj, 
			uint64 mem_size,
			filter_t *filter)
{

	int32 status;
	savefile_t *savefile = &obj->savefile;
	char db_name[24];

	uint32 num_db = 0;
	size_t batch_size;
	uint64 cache_size;

	uint8 buf[COMPRESSED_P_MAX_SIZE];
	relation_t tmp_rel;

	DB * reldb;
	DB * ideal_db;
	DBC * read_curs;
	DB_ENV * filter_env;
	DBT tmp_key;
	DBT store_buf;
	void *store_ptr;
	DBT read_buf;
	void *read_ptr;

	logprintf(obj, "commencing singleton removal\n");

	tmp_rel.factors = buf;

	/* we divide the input memory bound into an allocation for
	   Berkeley DB, and an allocation for the bulk buffers that
	   read from and write to it. We split it 50-50, with an 
	   upper limit of 3GB for the bulk buffers */

	if (0.5 * mem_size > 0xc0000000)
		batch_size = 0xc0000000;
	else
		batch_size = 0.5 * mem_size;

	cache_size = mem_size - batch_size; 

	/* pass 1: break up the collection of ideals into temporary
	   databases, each about the size of the DB cache and
	   each in sorted order */

	sprintf(buf, "%s.filter", savefile->name);
	filter_env = init_filter_env(obj, buf, cache_size);
	if (filter_env == NULL) {
		printf("error opening filtering DB env\n");
		exit(-1);
	}

	/* open relation DB for reading; do *not* make it part
	   of the DB environment, as we will be streaming relations
	   and don't want them evicting cache pages that are being
	   written with ideals */

	sprintf(buf, "%s.filter/relations.db", savefile->name);
	reldb = init_relation_db(obj, NULL, buf, DB_RDONLY);
	if (reldb == NULL) {
		printf("error opening relation DB\n");
		exit(-1);
	}

	status = reldb->cursor(reldb, NULL, 
				&read_curs, DB_CURSOR_BULK);
	if (status != 0) {
		logprintf(obj, "read cursor open failed: %s\n",
					db_strerror(status));
		exit(-1);
	}

	memset(&read_buf, 0, sizeof(DBT));
	read_buf.ulen = BULK_BUF_ROUND(batch_size / 2);
	read_buf.flags = DB_DBT_USERMEM | DB_DBT_BULK;
	read_buf.data = xmalloc(read_buf.ulen);
	DB_MULTIPLE_INIT(read_ptr, &read_buf);
	memset(&tmp_key, 0, sizeof(tmp_key));
	status = read_curs->get(read_curs, 
				&tmp_key, &read_buf, 
				DB_NEXT | DB_MULTIPLE_KEY);

	/* open first ideal DB to write */

	sprintf(db_name, "tmpdb.%u", num_db++);
	ideal_db = init_ideal_db(obj, filter_env, db_name, DB_CREATE);
	if (ideal_db == NULL) {
		printf("error opening ideal DB\n");
		exit(-1);
	}

	memset(&store_buf, 0, sizeof(DBT));
	store_buf.ulen = BULK_BUF_ROUND(batch_size / 2);
	store_buf.flags = DB_DBT_USERMEM | DB_DBT_BULK;
	store_buf.data = xmalloc(store_buf.ulen);
	DB_MULTIPLE_WRITE_INIT(store_ptr, &store_buf);

	while (1) {

		uint32 i;
		relation_lp_t tmp_ideal;
		uint8 * curr_key;
		size_t key_size;
		uint8 * curr_data;
		size_t data_size;

		/* read next relation; if batch is empty, read the
		   next batch */

		DB_MULTIPLE_KEY_NEXT(read_ptr, &read_buf,
					curr_key, key_size,
					curr_data, data_size);

		if (read_ptr == NULL) {
			status = read_curs->get(read_curs, 
					&tmp_key, &read_buf, 
					DB_NEXT | DB_MULTIPLE_KEY);
			if (status != 0) {
				break;
			}
			else {
				DB_MULTIPLE_INIT(read_ptr, &read_buf);
				DB_MULTIPLE_KEY_NEXT(read_ptr, &read_buf,
						curr_key, key_size,
						curr_data, data_size);
			}
		}

		/* pack the retrieved results into a relation_t and
		   get its ideal decomposition */

		memcpy(&tmp_rel.a, curr_key, sizeof(int64));
		memcpy(&tmp_rel.b, curr_key + sizeof(int64), sizeof(uint32));
		tmp_rel.num_factors_r = curr_data[0];
		tmp_rel.num_factors_a = curr_data[1];
		memcpy(tmp_rel.factors, curr_data + 2, data_size - 2);

		find_large_ideals(&tmp_rel, &tmp_ideal,
				filter->filtmin_r, 
				filter->filtmin_a);

		/* queue each large ideal for DB write */

		for (i = 0; i < tmp_ideal.ideal_count; i++) {

			DB_MULTIPLE_KEY_WRITE_NEXT(store_ptr, &store_buf,
					tmp_ideal.ideal_list + i,
					sizeof(ideal_t), NULL, 0);

			if (store_ptr == NULL) {

				/* ideal doesn't fit, commit the batch and
				   start the next. Note that bulk writes
				   to a compressed DB will not work unless
				   the batch is pre-sorted */

				status = ideal_db->sort_multiple(ideal_db, &store_buf,
						NULL, DB_MULTIPLE_KEY);
				if (status != 0) {
					logprintf(obj, "DB bulk sort 3 failed: %s\n",
							db_strerror(status));
					exit(-1);
				}

				status = ideal_db->put(ideal_db, NULL, &store_buf,
						NULL, DB_MULTIPLE_KEY);
				if (status != 0) {
					logprintf(obj, "DB bulk write 3 failed: %s\n",
							db_strerror(status));
					exit(-1);
				}

				DB_MULTIPLE_WRITE_INIT(store_ptr, &store_buf);

				DB_MULTIPLE_KEY_WRITE_NEXT(store_ptr, &store_buf,
						tmp_ideal.ideal_list + i,
						sizeof(ideal_t), NULL, 0);

				/* if the number of committed relations overflows
				    the cache, flush the DB and start another */

				if (db_fills_cache(filter_env, db_name, 2000)) {

					/* the DB is small enough to fit in cache,
					   so compacting it is a cheap operation. 
					   Berkeley DB only fills up about 75% of
					   disk pages, so we save space and improve
					   access times by compacting */ 

					status = ideal_db->compact(ideal_db, 
								NULL, NULL, 
								NULL, NULL, 
								DB_FREE_SPACE, NULL);
					if (status != 0) {
						logprintf(obj, "compact 3 failed:%s\n",
								db_strerror(status));
						exit(-1);
					}

					status = ideal_db->close(ideal_db, 0);
					if (status != 0) {
						logprintf(obj, "DB close failed: %s\n",
								db_strerror(status));
						exit(-1);
					}

					sprintf(db_name, "tmpdb.%u", num_db++);

					ideal_db = init_relation_db(obj, filter_env, 
								db_name, DB_CREATE);
					if (ideal_db == NULL) {
						printf("error opening ideal DB\n");
						exit(-1);
					}
				}
			}
		}
	}

	/* commit the final batch */

	status = ideal_db->sort_multiple(ideal_db, &store_buf,
			NULL, DB_MULTIPLE_KEY);
	if (status != 0) {
		logprintf(obj, "DB bulk sort 4 failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	status = ideal_db->put(ideal_db, NULL, &store_buf,
				NULL, DB_MULTIPLE_KEY);
	if (status != 0) {
		logprintf(obj, "DB bulk write 4 failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	status = ideal_db->compact(ideal_db, NULL, NULL, 
				NULL, NULL, 
				DB_FREE_SPACE, NULL);
	if (status != 0) {
		logprintf(obj, "DB compact 4 failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	status = ideal_db->close(ideal_db, 0);
	if (status != 0) {
		logprintf(obj, "DB close 4 failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	status = reldb->close(reldb, 0);
	if (status != 0) {
		logprintf(obj, "DB close 5 failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	/* free pass 1 stuff */

	free(read_buf.data);
	free(store_buf.data);
	free_filter_env(obj, filter_env);

	/* pass 2: merge the temporary DBs into a single output one */

	merge_ideal_files(obj, mem_size,
			num_db, &filter->num_ideals);
}
