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

#define AB_KEY_SIZE (sizeof(int64) + sizeof(uint32))

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
int compare_relnums(DB *db, const DBT *rel1,
			const DBT *rel2)
{
	uint64 i1, i2;

	memcpy(&i1, rel1->data, sizeof(int64));
	memcpy(&i2, rel2->data, sizeof(int64));

	if (i1 > i2)
		return 1;
	if (i1 < i2)
		return -1;

	return 0;
}

/*--------------------------------------------------------------------*/
static int compress_relnum(DB *db, 
			const DBT *prev_key, const DBT *prev_data, 
			const DBT *key, const DBT *data, 
			DBT *dest)
{
	uint64 i1, i2;
	uint32 count = 0;
	uint8 tmp[24];

	memcpy(&i1, (uint8 *)prev_key->data, sizeof(uint64));
	memcpy(&i2, (uint8 *)key->data, sizeof(uint64));

	count = compress_p(tmp, i2 - i1, count);

	dest->size = count;
	if (dest->ulen < dest->size)
		return DB_BUFFER_SMALL;

	memcpy(dest->data, tmp, count);
	return 0;
}

/*--------------------------------------------------------------------*/
static int decompress_relnum(DB *db, 
		const DBT *prev_key, const DBT *prev_data, 
		DBT *compressed, 
		DBT *key, DBT *data)
{
	uint64 i2;
	uint32 count = 0;

	memcpy(&i2, (uint8 *)prev_key->data, sizeof(uint64));
	i2 += decompress_p((uint8 *)compressed->data, &count);

	key->size = sizeof(uint64);
	data->size = 0;
	if (key->ulen < key->size || data->ulen < data->size)
		return DB_BUFFER_SMALL;

	memcpy(key->data, &i2, sizeof(uint64));
	compressed->size = count;
	return 0;
}

/*--------------------------------------------------------------------*/
static int compare_ideals(DB *db, const DBT *i1,
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
static void
purge_pass2(msieve_obj *obj,
		void *pass1_stream,
		size_t mem_size,
		uint64 num_relations_in)
{
	char db_name[LINE_BUF_SIZE];
	uint64 num_relations;
	size_t read_batch_size;
	size_t write_batch_size;
	void *write_stream;
	DBT *read_key;
	DBT *read_data;
	uint8 prev_key[AB_KEY_SIZE] = {0};

	logprintf(obj, "commencing duplicate removal, pass 2\n");

	sprintf(db_name, "%s.filter/tmpi", obj->savefile.name);

	if (.5 * mem_size > 0xc0000000)
		read_batch_size = 0xc0000000;
	else
		read_batch_size = 0.5 * mem_size;
	write_batch_size = read_batch_size;

	write_stream = stream_db_init(obj, NULL, 
				compare_relnums,
				compress_relnum,
				decompress_relnum,
				db_name);
	if (write_stream == NULL) {
		logprintf(obj, "error opening write stream 2\n");
		exit(-1);
	}

	stream_db_write_init(write_stream, 0, write_batch_size);

	stream_db_read_init(pass1_stream, 0,
				read_batch_size,
				&read_key, &read_data);
	num_relations = 0;
	while (read_key != NULL && read_data != NULL) {

		if (memcmp(read_key->data, prev_key, AB_KEY_SIZE) != 0) {

			uint64 relnum;
			uint32 count = 0;

			relnum = decompress_p(read_data->data, &count);

			stream_db_write_next(write_stream,
					&relnum, sizeof(uint64),
					NULL, 0);

			memcpy(prev_key, read_key->data, AB_KEY_SIZE);
			num_relations++;
		}

		stream_db_read_next(pass1_stream, 
				&read_key, &read_data);
	}

	stream_db_write_close(write_stream);
	stream_db_read_close(pass1_stream);

	logprintf(obj, "wrote %" PRIu64 " duplicates and %"
			PRIu64 " unique relations\n", 
			num_relations_in - num_relations, 
			num_relations);

	stream_db_free(write_stream);
}

/*--------------------------------------------------------------------*/
void nfs_purge_duplicates(msieve_obj *obj, factor_base_t *fb,
				uint64 mem_size,
				uint64 max_relations)
{

	savefile_t *savefile = &obj->savefile;
	char buf[LINE_BUF_SIZE];
	char db_name[LINE_BUF_SIZE];
	uint64 curr_relation;
	uint64 num_relations;
	size_t batch_size;
	void *write_stream;

	logprintf(obj, "commencing duplicate removal, pass 1\n");

	savefile_open(savefile, SAVEFILE_READ);

	/* upper limit of 3GB for bulk writes */

	if (mem_size > 0xc0000000)
		batch_size = 0xc0000000;
	else
		batch_size = mem_size;

	/* pass 1: break up the input dataset into temporary
	   databases, each about the size of the DB cache and
	   each in sorted order */

	sprintf(db_name, "%s.filter/tmprel", savefile->name);

	write_stream = stream_db_init(obj, NULL, 
				compare_relations,
				compress_relation,
				decompress_relation,
				db_name);
	if (write_stream == NULL) {
		logprintf(obj, "error opening write stream\n");
		exit(-1);
	}

	num_relations = 0;
	curr_relation = 0;
	savefile_read_line(buf, sizeof(buf), savefile);
	stream_db_write_init(write_stream, 0, batch_size);

	while (!savefile_eof(savefile)) {

		uint8 db_key[AB_KEY_SIZE];
		uint8 db_data[12];
		uint32 relnum_count = 0;
		char *ptr;
		int64 a;
		uint32 b;

		if (buf[0] != '-' && !isdigit(buf[0])) {

			/* no relation on this line */

			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}

		/* read the relation coordinates (skip free relations) */

		curr_relation++;
		if (max_relations && curr_relation >= max_relations)
			break;

		ptr = buf;
		a = strtoll(ptr, &ptr, 10);

		if (ptr[0] != ',') {
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}

		b = strtoul(ptr + 1, &ptr, 10);

		if (ptr[0] != ':' || b == 0) {
			savefile_read_line(buf, sizeof(buf), savefile);
			continue;
		}

		if (curr_relation % 10000000 == 0) {
			printf("read %" PRIu64 "M relations\n", 
					curr_relation / 1000000);
		}

		/* relation coords are good; queue for DB write */

		num_relations++;
		relnum_count = compress_p(db_data, curr_relation, relnum_count);
		memcpy(db_key, &a, sizeof(int64));
		memcpy(db_key + sizeof(int64), &b, sizeof(uint32));

		stream_db_write_next(write_stream,
				db_key, AB_KEY_SIZE,
				db_data, relnum_count);

		/* get the next line */

		savefile_read_line(buf, sizeof(buf), savefile);
	}

	/* free pass 1 stuff */

	savefile_close(savefile);
	stream_db_write_close(write_stream);

	/* pass 2: merge the temporary DBs into a single output one */

	logprintf(obj, "found %" PRIu64 " relations\n", num_relations);

	purge_pass2(obj, write_stream, mem_size, num_relations);

	stream_db_free(write_stream);
}
