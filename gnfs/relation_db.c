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
#include "relation_db.h"

#define BULK_BUF_ROUND(x) ((x) & ~(1024 - 1))

typedef struct {
	DB *read_db;
	DBC *read_curs;
	void *read_ptr;
	DBT read_key;
	DBT read_data;
	DBT curr_key;
	DBT curr_data;
} tmp_db_t;

typedef struct {
	msieve_obj *obj;
	DB_ENV *env;
	char *name_prefix;
	compare_cb compare_fcn;
	compress_cb compress_fcn;
	decompress_cb decompress_fcn;

	uint32 num_db;
	uint32 curr_num_db;
	uint32 curr_db;
	uint32 db_flags;

	DB * write_db;
	DBT store_buf;
	void *store_ptr;

	tmp_db_t *tmp_db;
	tmp_db_t **tmp_db_ptr;
} stream_db_t;

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
		uint32 db_flags,
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

	if (compress_fcn != NULL && decompress_fcn != NULL) {
		status = db->set_bt_compress(db, compress_fcn,
					decompress_fcn);
		if (status != 0) {
			logprintf(obj, "DB set compress fcn failed: %s\n",
					db_strerror(status));
			goto error_finished;
		}
	}
	
	if (filter_env && (open_flags & DB_TRUNCATE)) {
		/* strip out DB_TRUNCATE, environments do not allow it */
		open_flags &= ~DB_TRUNCATE;
		status = filter_env->dbremove(filter_env, 
					NULL, name, NULL, 0);
		if (status != 0 && status != ENOENT) {
			logprintf(obj, "DB delete failed: %s\n",
					db_strerror(status));
			goto error_finished;
		}
	}

	if (db_flags != 0) {
		status = db->set_flags(db, db_flags);
		if (status != 0) {
			logprintf(obj, "DB set flags failed: %s\n",
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
void * stream_db_init(msieve_obj *obj, 
			DB_ENV *filter_env, 
			uint32 num_existing_files,
			compare_cb compare,
			compress_cb compress,
			decompress_cb decompress,
			char *name_prefix)
{
	stream_db_t *s = (stream_db_t *)xcalloc(1, sizeof(stream_db_t));

	if (s != NULL) {
		s->obj = obj;
		s->name_prefix = strdup(name_prefix);
		s->env = filter_env;
		s->compare_fcn = compare;
		s->compress_fcn = compress;
		s->decompress_fcn = decompress;
		s->num_db = num_existing_files;
	}
	return s;
}

/*--------------------------------------------------------------------*/
void stream_db_free(void *s_in)
{
	stream_db_t *s = (stream_db_t *)s_in;

	free(s->name_prefix);
	free(s);
}

/*--------------------------------------------------------------------*/
uint32 stream_db_num_files(void *s_in)
{
	stream_db_t *s = (stream_db_t *)s_in;

	return s->num_db;
}

/*--------------------------------------------------------------------*/
/* boilerplate code for managing a heap of DB items */

#define HEAP_SWAP(a,b) { tmp = a; a = b; b = tmp; }
#define HEAP_PARENT(i)  (((i)-1) >> 1)
#define HEAP_LEFT(i)    (2 * (i) + 1)
#define HEAP_RIGHT(i)   (2 * (i) + 2)

static void heapify(tmp_db_t **h, uint32 index, uint32 size,
			compare_cb compare_fcn) 
{
	uint32 c;
	tmp_db_t * tmp;
	for (c = HEAP_LEFT(index); c < (size-1); 
			index = c, c = HEAP_LEFT(index)) {

		if (compare_fcn(NULL, &h[c]->curr_key,
				&h[c+1]->curr_key) >= 0)
			c++;

		if (compare_fcn(NULL, &h[index]->curr_key,
				&h[c]->curr_key) >= 0) {
			HEAP_SWAP(h[index], h[c]);
		}
		else
			return;
	}
	if (c == (size-1) && 
	    compare_fcn(NULL, &h[index]->curr_key,
				&h[c]->curr_key) >= 0) {
		HEAP_SWAP(h[index], h[c]);
	}
}

static void make_heap(tmp_db_t **h, uint32 size,
			compare_cb compare_fcn) 
{
	uint32 i;
	for (i = HEAP_PARENT(size); i; i--)
		heapify(h, i-1, size, compare_fcn);
}

/*--------------------------------------------------------------------*/
void stream_db_read_init(void *s_in,
		uint32 db_flags,
		uint64 buffer_size,
		DBT **first_key, 
		DBT **first_data)
{
	stream_db_t *s = (stream_db_t *)s_in;
	uint32 num_db = s->num_db;
	int32 status;
	uint32 i, j;

	if (num_db > 1)
		logprintf(s->obj, "merging %u temporary databases\n", num_db);

	s->db_flags = db_flags;
	s->tmp_db = (tmp_db_t *)xcalloc(num_db, sizeof(tmp_db_t));
	s->tmp_db_ptr = (tmp_db_t **)xcalloc(num_db, sizeof(tmp_db_t *));

	buffer_size = buffer_size / num_db;
	buffer_size = MAX(buffer_size, 1 << 20);
	buffer_size = MIN(buffer_size, 0xe0000000);

	for (i = 0; i < num_db; i++) {

		char db_name[LINE_BUF_SIZE];
		tmp_db_t *t = s->tmp_db + i;

		sprintf(db_name, "%s.%u", s->name_prefix, i);
		t->read_db = init_db_core(s->obj, s->env, db_name, 
					DB_RDONLY, s->db_flags,
					s->compare_fcn,
					s->compress_fcn,
					s->decompress_fcn);
		if (t->read_db == NULL) {
			logprintf(s->obj, "error opening input DB\n");
			exit(-1);
		}

		status = t->read_db->cursor(t->read_db, NULL, 
					&t->read_curs, 
					DB_CURSOR_BULK);
		if (status != 0) {
			logprintf(s->obj, "read cursor open failed: %s\n",
						db_strerror(status));
			exit(-1);
		}

		t->read_data.ulen = BULK_BUF_ROUND(buffer_size);
		t->read_data.flags = DB_DBT_USERMEM | DB_DBT_BULK;
		t->read_data.data = xmalloc(t->read_data.ulen);
	}

	/* do the initial read from each DB. Squeeze out any
	   DBs that cannot deliver, heapify all the ones that can */

	for (i = j = 0; i < num_db; i++) {

		tmp_db_t *t = s->tmp_db + i;

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
			s->tmp_db_ptr[j++] = t;
		}
	}

	if (j != 0) {
		if (j > 1)
			make_heap(s->tmp_db_ptr, j, s->compare_fcn);

		*first_key = &s->tmp_db_ptr[0]->curr_key;
		*first_data = &s->tmp_db_ptr[0]->curr_data;
		s->curr_num_db = j;
		s->curr_db = 0;
	}
	else {
		*first_key = NULL;
		*first_data = NULL;
	}
}

/*--------------------------------------------------------------------*/
void stream_db_read_next(void *s_in,
		DBT **next_key, DBT **next_data)
{
	stream_db_t *s = (stream_db_t *)s_in;
	int32 status;
	uint32 i = s->curr_db;
	uint32 j = s->curr_num_db;
	tmp_db_t *t;

	if (j == i) {
		*next_key = NULL;
		*next_data = NULL;
		return;
	}

	t = s->tmp_db_ptr[i];
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
			s->curr_db = ++i;
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

	if (j == i) {
		*next_key = NULL;
		*next_data = NULL;
	}
	else {
		if (j - i > 1) {
			heapify(s->tmp_db_ptr + i, 0, j - i, 
					s->compare_fcn);
		}

		*next_key = &s->tmp_db_ptr[i]->curr_key;
		*next_data = &s->tmp_db_ptr[i]->curr_data;
	}
}

/*--------------------------------------------------------------------*/
void stream_db_read_close(void *s_in, uint32 delete_files)
{
	stream_db_t *s = (stream_db_t *)s_in;
	uint32 i;

	for (i = 0; i < s->num_db; i++) {

		tmp_db_t *t = s->tmp_db + i;

		free(t->read_data.data);

		t->read_curs->close(t->read_curs);
		t->read_db->close(t->read_db, 0);

		if (delete_files) {
			char db_name[LINE_BUF_SIZE];

			sprintf(db_name, "%s.%u", s->name_prefix, i);
			remove(db_name);
		}
	}
	free(s->tmp_db);
	free(s->tmp_db_ptr);
}

/*--------------------------------------------------------------------*/
void stream_db_write_init(void *s_in,
		uint32 db_flags,
		uint64 buffer_size)
{
	stream_db_t *s = (stream_db_t *)s_in;
	char db_name[LINE_BUF_SIZE];

	buffer_size = MIN(buffer_size, 0xe0000000);

	s->num_db = 0;
	s->db_flags = db_flags;
	sprintf(db_name, "%s.%u", s->name_prefix, s->num_db++);

	s->write_db = init_db_core(s->obj, s->env, db_name, 
				DB_CREATE | DB_TRUNCATE, s->db_flags,
				s->compare_fcn,
				s->compress_fcn,
				s->decompress_fcn);
	if (s->write_db == NULL) {
		printf("error opening DB\n");
		exit(-1);
	}

	memset(&s->store_buf, 0, sizeof(DBT));
	s->store_buf.ulen = BULK_BUF_ROUND(buffer_size);
	s->store_buf.flags = DB_DBT_USERMEM | DB_DBT_BULK;
	s->store_buf.data = xmalloc(s->store_buf.ulen);
	DB_MULTIPLE_WRITE_INIT(s->store_ptr, &s->store_buf);
}

/*--------------------------------------------------------------------*/
void stream_db_write_next(void *s_in,
		void *key, uint32 key_size,
		void *data, uint32 data_size)
{
	stream_db_t *s = (stream_db_t *)s_in;
	int32 status;

	DB_MULTIPLE_KEY_WRITE_NEXT(s->store_ptr, 
				&s->store_buf,
				key, key_size,
				data, data_size);

	if (s->store_ptr == NULL) {

		char db_name[LINE_BUF_SIZE];

		/* entry doesn't fit, commit the batch and open a new DB */

		status = s->write_db->sort_multiple(s->write_db, 
					&s->store_buf,
					NULL, DB_MULTIPLE_KEY);
		if (status != 0) {
			logprintf(s->obj, "DB bulk sort failed: %s\n",
					db_strerror(status));
			exit(-1);
		}

		status = s->write_db->put(s->write_db, 
					NULL, &s->store_buf,
					NULL, DB_MULTIPLE_KEY);
		if (status != 0) {
			logprintf(s->obj, "DB bulk write failed: %s\n",
					db_strerror(status));
			exit(-1);
		}

		DB_MULTIPLE_WRITE_INIT(s->store_ptr, &s->store_buf);

		DB_MULTIPLE_KEY_WRITE_NEXT(s->store_ptr, 
					&s->store_buf,
					key, key_size,
					data, data_size);

		status = s->write_db->close(s->write_db, 0);
		if (status != 0) {
			logprintf(s->obj, "DB close failed: %s\n",
					db_strerror(status));
			exit(-1);
		}

		sprintf(db_name, "%s.%u", 
				s->name_prefix, s->num_db++);

		s->write_db = init_db_core(s->obj, s->env, db_name, 
				DB_CREATE | DB_TRUNCATE, s->db_flags,
				s->compare_fcn,
				s->compress_fcn,
				s->decompress_fcn);
		if (s->write_db == NULL) {
			printf("error opening new DB\n");
			exit(-1);
		}
	}
}

/*--------------------------------------------------------------------*/
void stream_db_write_close(void *s_in)
{
	stream_db_t *s = (stream_db_t *)s_in;
	int32 status;

	status = s->write_db->sort_multiple(s->write_db, 
				&s->store_buf,
				NULL, DB_MULTIPLE_KEY);
	if (status != 0) {
		logprintf(s->obj, "DB bulk sort failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	status = s->write_db->put(s->write_db, 
				NULL, &s->store_buf,
				NULL, DB_MULTIPLE_KEY);
	if (status != 0) {
		logprintf(s->obj, "DB bulk write failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	status = s->write_db->close(s->write_db, 0);
	if (status != 0) {
		logprintf(s->obj, "DB close failed: %s\n",
				db_strerror(status));
		exit(-1);
	}

	s->write_db = NULL;
	free(s->store_buf.data);
}

