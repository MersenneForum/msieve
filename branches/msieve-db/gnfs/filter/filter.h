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

/* implementation of number field sieve filtering */

#ifndef _GNFS_FILTER_FILTER_H_
#define _GNFS_FILTER_FILTER_H_

#include <common/filter/filter.h>
#include "gnfs.h"

#ifdef __cplusplus
extern "C" {
#endif

int compare_relnums(DB *db, const DBT *rel1,
			const DBT *rel2);

int compress_relnum(DB *db, 
		const DBT *prev_key, const DBT *prev_data, 
		const DBT *key, const DBT *data, 
		DBT *dest);

int decompress_relnum(DB *db, 
		const DBT *prev_key, const DBT *prev_data, 
		DBT *compressed, 
		DBT *key, DBT *data);

/* create a set of databases containing the line numbers of unique 
   relations. Duplicate removal only applies to the first max_relations 
   relations found (or all relations if zero). The number of files
   containing line numbers is returned. Relations only have cursory
   checks applied, the singleton removal does full verification.
   The return value is the automatically computed filtering bound */

uint32 nfs_purge_duplicates(msieve_obj *obj,
				uint64 mem_size,
				uint64 max_relations,
				uint32 *num_line_num_files);

/* read '<savefile_name>.d' and create '<savefile_name>.lp', a 
   binary file containing the relations surviving the singleton
   removal pass. If pass = 0, the .d file is assumed to contain
   relation numbers to skip; otherwise it contains relation numbers
   to keep */
   
void nfs_write_lp_file(msieve_obj *obj, factor_base_t *fb,
			uint64 mem_size,
			uint32 num_line_num_files,
			filter_t *filter);

#ifdef __cplusplus
}
#endif

#endif /* _GNFS_FILTER_FILTER_H_ */
