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

#ifndef _RELATION_DB_H
#define _RELATION_DB_H

#include <db.h>
#include <common.h>

#ifdef __cplusplus
extern "C"
{
#endif

DB_ENV * init_filter_env(msieve_obj *obj, char *dirname,
				uint64 cache_size);

void free_filter_env(msieve_obj *obj, DB_ENV *env);

/* placeholders for callbacks */

typedef int (*compare_cb)(DB *db, const DBT *rel1, const DBT *rel2);

typedef int (*compress_cb)(DB *db, 
			const DBT *prev_key, const DBT *prev_data, 
			const DBT *key, const DBT *data, 
			DBT *dest);

typedef int (*decompress_cb)(DB *db, 
			const DBT *prev_key, const DBT *prev_data, 
			DBT *compressed, 
			DBT *key, DBT *data);

/* batch stream interface to databases */

void * stream_db_init(msieve_obj *obj, 
			DB_ENV *filter_env, 
			uint32 num_existing_files,
			compare_cb compare,
			compress_cb compress,
			decompress_cb decompress,
			char *name_prefix);

void stream_db_free(void *s);

uint32 stream_db_num_files(void *s);

/* batch write interface */

void stream_db_write_init(void *s,
		uint32 db_flags,
		size_t buffer_size);

void stream_db_write_next(void *s,
		void *key, uint32 key_size,
		void *data, uint32 data_size);

void stream_db_write_close(void *s);

/* batch read interface */

void stream_db_read_init(void *s,
		uint32 db_flags,
		size_t buffer_size,
		DBT **first_key, 
		DBT **first_data);

void stream_db_read_next(void *s,
		DBT **next_key, 
		DBT **next_data);

void stream_db_read_close(void *s, uint32 delete_files);

#ifdef __cplusplus
}
#endif

#endif /* _RELATION_DB_H */
