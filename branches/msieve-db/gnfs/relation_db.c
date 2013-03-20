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

#include <errno.h>
#include "gnfs.h"

typedef int (*compare_cb)(DB *db, const DBT *rel1, const DBT *rel2);

typedef int (*compress_cb)(DB *db, 
			const DBT *prev_key, const DBT *prev_data, 
			const DBT *key, const DBT *data, 
			DBT *dest);

typedef int (*decompress_cb)(DB *db, 
			const DBT *prev_key, const DBT *prev_data, 
			DBT *compressed, 
			DBT *key, DBT *data);

/*--------------------------------------------------------------------*/
uint32 db_fills_cache(DB_ENV *env, char *db_name, uint32 min_evict)
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

	return (evicted > min_evict);
}

/*--------------------------------------------------------------------*/
int compare_relations(DB *db, const DBT *rel1,
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

	key->size = AB_KEY_SIZE;
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
int compare_ideals(DB *db, const DBT *i1,
			const DBT *i2)
{
	/* Ordering is by prime, then by root of prime,
	   then by rational or algebraic type. This ordering
	   is designed to put the smallest ideals first */

	ideal_t a, b;

	memcpy(&a, i1->data, sizeof(ideal_t));
	memcpy(&b, i2->data, sizeof(ideal_t));

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

	memcpy(&i1, prev_key->data, sizeof(ideal_t));
	memcpy(&i2, key->data, sizeof(ideal_t));

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
	uint32 data_size;

	r2 = decompress_p((uint8 *)compressed->data, &count);
	p_diff_zero = r2 & 1;
	i2.rat_or_alg = (r2 >> 1) & 1;
	i2.r = r2 >> 2;

	memcpy(&i1, prev_key->data, sizeof(ideal_t));
	i2.p = i1.p;
	if (!p_diff_zero)
		i2.p += decompress_p((uint8 *)compressed->data, &count);

	data_size = (uint32)decompress_p((uint8 *)compressed->data, &count);

	key->size = sizeof(ideal_t);
	data->size = data_size;
	if (key->ulen < key->size || data->ulen < data->size)
		return DB_BUFFER_SMALL;

	memcpy(data->data, (uint8 *)compressed->data + count, data_size);
	memcpy(key->data, &i2, sizeof(ideal_t));
	compressed->size = count + data_size;
	return 0;
}

/*--------------------------------------------------------------------*/
DB_ENV * init_filter_env(msieve_obj *obj, char *dirname,
				uint64 cache_size)
{
	DB_ENV * env = NULL;
	int32 status;
	uint32 flags = DB_CREATE | DB_INIT_MPOOL | DB_PRIVATE;

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

#if defined(WIN32) || defined(_WIN64)
	status = env->set_flags(env, DB_DIRECT_DB, 1);
	if (status != 0) {
		logprintf(obj, "can't do direct IO, try building BDB with "
			"--enable-direct_io\n");
	}
#endif

	/* pre-emptively make the environment directory; the open 
	   call will not fail if it doesn't exist */

	mkdir(dirname, 0755);
	status = env->open(env, dirname, flags, 0);
	if (status != 0) {
		logprintf(obj, "DB_ENV open failed: %s\n",
				db_strerror(status));
		goto error_finished;
	}

	/* do not cache databases in disk cache and also in
	   Berkeley DB's memory cache. Continue if the call
	   fails, not all OSes support it */

	return env;

error_finished:
	if (env)
		env->close(env, 0);
	return NULL;
}

/*--------------------------------------------------------------------*/
void free_filter_env(msieve_obj *obj, DB_ENV *env)
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
static DB * init_db_core(msieve_obj *obj, 
		DB_ENV *filter_env, char *name,
		uint32 open_flags,
		compare_cb compare_fcn,
		compress_cb compress_fcn,
		decompress_cb decompress_fcn)
{
	DB * db = NULL;
	int32 status;

	status = db_create(&db, filter_env, 0);
	if (status != 0) {
		logprintf(obj, "DB creation failed: %s\n",
				db_strerror(status));
		goto error_finished;
	}

	status = db->set_bt_compare(db, compare_fcn);
	if (status != 0) {
		logprintf(obj, "DB set compare fcn failed: %s\n",
				db_strerror(status));
		goto error_finished;
	}

	status = db->set_bt_compress(db, compress_fcn,
				decompress_fcn);
	if (status != 0) {
		logprintf(obj, "DB set compress fcn failed: %s\n",
				db_strerror(status));
		goto error_finished;
	}

	/* strip out DB_TRUNCATE, environments do not allow it */
	if (filter_env && (open_flags & DB_TRUNCATE)) {
		open_flags &= ~DB_TRUNCATE;
		status = filter_env->dbremove(filter_env, 
					NULL, name, NULL, 0);
		if (status != 0 && status != ENOENT) {
			logprintf(obj, "DB delete failed: %s\n",
					db_strerror(status));
			goto error_finished;
		}
	}

	status = db->open(db, NULL, name, NULL, DB_BTREE, 
				open_flags, 
				(open_flags & DB_RDONLY) ? 0444 : 0664);
	if (status != 0) {
		logprintf(obj, "DB open failed: %s\n",
				db_strerror(status));
		goto error_finished;
	}

	return db;

error_finished:
	if (db)
		db->close(db, 0);
	return NULL;
}

/*--------------------------------------------------------------------*/
DB * init_relation_db(msieve_obj *obj, 
			DB_ENV *filter_env, char *name,
			uint32 open_flags)
{
	return init_db_core(obj, filter_env, name,
				open_flags,
				compare_relations,
				compress_relation,
				decompress_relation);
}

/*--------------------------------------------------------------------*/
DB * init_ideal_db(msieve_obj *obj, 
			DB_ENV *filter_env, char *name,
			uint32 open_flags)
{
	return init_db_core(obj, filter_env, name,
				open_flags,
				compare_ideals,
				compress_ideal,
				decompress_ideal);
}
