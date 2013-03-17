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

#define AB_KEY_SIZE (sizeof(int64) + sizeof(uint32))

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

int compare_relations(DB *db, const DBT *rel1,
			const DBT *rel2);

int compare_ideals(DB *db, const DBT *rel1,
			const DBT *rel2);

DB_ENV * init_filter_env(msieve_obj *obj, char *dirname,
				uint64 cache_size);

void free_filter_env(msieve_obj *obj, DB_ENV *env);

DB * init_relation_db(msieve_obj *obj, 
			DB_ENV *filter_env, char *name,
			uint32 open_flags);

DB * init_ideal_db(msieve_obj *obj, 
			DB_ENV *filter_env, char *name,
			uint32 open_flags);

uint32 db_fills_cache(DB_ENV *env, char *db_name, uint32 min_evict);

#ifdef __cplusplus
}
#endif

#endif /* _RELATION_DB_H */
