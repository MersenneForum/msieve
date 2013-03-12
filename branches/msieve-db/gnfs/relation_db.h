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

#define KEY_SIZE (sizeof(int64) + sizeof(uint32))

int compare_relations(DB *db, const DBT *rel1,
			const DBT *rel2);

DB_ENV * init_filter_env(msieve_obj *obj, char *dirname,
				uint64 cache_size);

void free_filter_env(msieve_obj *obj, DB_ENV *env);

DB * init_relation_db(msieve_obj *obj, 
			DB_ENV *filter_env, char *name,
			uint32 open_flags);

#ifdef __cplusplus
}
#endif

#endif /* _RELATION_DB_H */
