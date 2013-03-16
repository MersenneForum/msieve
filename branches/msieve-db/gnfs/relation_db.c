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
	if (status == 0) {
		logprintf(obj, "can't do direct IO, try building BDB with "
			"--enable-direct_io\n");
	}
#elif defined(__linux__)
	/* Berkeley DB manages disk data in its own cache, and reads
	   and writes also go through the filesystem cache. On linux this 
	   means that by default, grinding through reading or writing a
	   large dataset will kick everything on your computer out of 
	   memory. The best way to stop this is not to store everything
	   in memory twice, which means using direct IO.
	   
	   Berkeley DB will not use direct IO in linux, because memory
	   buffers have to be aligned to a large boundary, which BDB
	   does not want to figure out on the fly. Making direct
	   IO work better in linux will never happen because it is
	   insufficiently aesthetically pleasing to Linus Torvalds. In
	   the meantime, the following is the least bad solution */

	logprintf(obj, "Make sure to do 'echo 0 > /proc/sys/vm/swappiness'\n");
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
DB * init_relation_db(msieve_obj *obj, 
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

	/* assume DB_CREATE means the DB should be truncated */
	if (open_flags & DB_CREATE) {
		status = filter_env->dbremove(filter_env, 
					NULL, name, NULL, 0);
		if (status != 0 && status != ENOENT) {
			logprintf(obj, "DB delete failed: %s\n",
					db_strerror(status));
			goto error_finished;
		}
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
