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
static uint32 pack_relation(uint8 *array, int64 a, uint32 b)
{
	uint32 a_neg = 0;
	uint32 count = 0;

	if (a < 0) {
		a_neg = 1;
		a = -a;
	}
	count = compress_p(array, a << 1 | a_neg, count);
	count = compress_p(array, b, count);
	return count;
}

static void unpack_relation(uint8 *array, 
			int64 *a_out, uint32 *b_out)
{
	uint32 count = 0;
	int64 a = decompress_p(array, &count);
	uint32 b = decompress_p(array, &count);

	if (a & 1)
		a = -a;

	*a_out = a >> 1;
	*b_out = b;
}

/*--------------------------------------------------------------------*/
int compare_relations(DB *db, const DBT *rel1,
			const DBT *rel2)
{
	int64 a1, a2;
	uint32 b1, b2;

	unpack_relation((uint8 *)rel1->data, &a1, &b1);
	unpack_relation((uint8 *)rel2->data, &a2, &b2);

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
	int64 a1, a2;
	uint32 b1, b2;
	uint32 count = 0;
	uint8 a_neg = 0;
	uint8 b_diff_zero = 0;
	uint8 tmp[24];

	unpack_relation((uint8 *)prev_key->data, &a1, &b1);
	unpack_relation((uint8 *)key->data, &a2, &b2);

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
	uint64 a1, a2;
	uint32 b2;
	uint32 b_diff_zero;
	uint32 a_neg;
	uint32 count = 0;
	uint32 data_size;
	uint8 tmp[24];

	a2 = decompress_p((uint8 *)compressed->data, &count);
	a_neg = a2 & 2;
	b_diff_zero = a2 & 1;
	a2 >>= 2;
	if (a_neg)
		a2 = -a2;

	unpack_relation((uint8 *)prev_key->data, (uint64 *)&a1, &b2);
	if (!b_diff_zero)
		b2 += (uint32)decompress_p((uint8 *)compressed->data, &count);

	data_size = (uint32)decompress_p((uint8 *)compressed->data, &count);

	key->size = pack_relation(tmp, a2, b2);
	data->size = data_size;
	if (key->ulen < key->size || data->ulen < data->size)
		return DB_BUFFER_SMALL;

	memcpy(data->data, (uint8 *)compressed->data + count, data_size);
	memcpy(key->data, tmp, key->size);

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
	uint64 i1 = decompress_one_p((uint8 *)rel1->data);
	uint64 i2 = decompress_one_p((uint8 *)rel2->data);

	if (i1 > i2)
		return 1;
	if (i1 < i2)
		return -1;

	return 0;
}

/*--------------------------------------------------------------------*/
int compress_relnum(DB *db, 
		const DBT *prev_key, const DBT *prev_data, 
		const DBT *key, const DBT *data, 
		DBT *dest)
{
	uint8 tmp[24];
	uint64 i1 = decompress_one_p((uint8 *)prev_key->data);
	uint64 i2 = decompress_one_p((uint8 *)key->data);
	uint32 count = compress_one_p(tmp, i2 - i1);

	dest->size = count;
	if (dest->ulen < dest->size)
		return DB_BUFFER_SMALL;

	memcpy(dest->data, tmp, count);
	return 0;
}

/*--------------------------------------------------------------------*/
int decompress_relnum(DB *db, 
		const DBT *prev_key, const DBT *prev_data, 
		DBT *compressed, 
		DBT *key, DBT *data)
{
	uint32 count;
	uint32 c_count = 0;
	uint8 tmp[24];
	uint64 i2 = decompress_one_p((uint8 *)prev_key->data);

	i2 += decompress_p((uint8 *)compressed->data, &c_count);
	count = compress_one_p(tmp, i2);

	key->size = count;
	data->size = 0;
	if (key->ulen < key->size || data->ulen < data->size)
		return DB_BUFFER_SMALL;

	memcpy(key->data, tmp, count);
	compressed->size = c_count;
	return 0;
}

/*--------------------------------------------------------------------*/
static void
purge_pass2(msieve_obj *obj,
		void *pass1_stream,
		size_t mem_size,
		uint64 num_relations_in,
		uint32 * num_line_num_files)
{
	char db_name[LINE_BUF_SIZE];
	uint64 num_relations;
	size_t read_batch_size;
	size_t write_batch_size;
	void *write_stream;
	DBT *read_key;
	DBT *read_data;
	uint8 prev_key[24] = {0};
	uint32 prev_key_size = 0;

	logprintf(obj, "commencing duplicate removal, pass 2\n");

	sprintf(db_name, "%s.filter/tmpi", obj->savefile.name);

	if (0.5 * mem_size > 0xc0000000)
		read_batch_size = 0xc0000000;
	else
		read_batch_size = 0.5 * mem_size;
	write_batch_size = read_batch_size;

	write_stream = stream_db_init(obj, NULL, 0,
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

		if (prev_key_size != read_key->size ||
		    memcmp(prev_key, read_key->data, read_key->size) != 0) {

			stream_db_write_next(write_stream,
					read_data->data, read_data->size,
					NULL, 0);

			memcpy(prev_key, read_key->data, read_key->size);
			prev_key_size = read_key->size;
			num_relations++;
		}

		stream_db_read_next(pass1_stream, 
				&read_key, &read_data);
	}

	stream_db_write_close(write_stream);
	stream_db_read_close(pass1_stream, 1);

	logprintf(obj, "wrote %" PRIu64 " duplicates and %"
			PRIu64 " unique relations\n", 
			num_relations_in - num_relations, 
			num_relations);

	*num_line_num_files = stream_db_num_files(write_stream);
	stream_db_free(write_stream);
}

/*--------------------------------------------------------------------*/
#define LOG2_BIN_SIZE 17
#define BIN_SIZE (1 << (LOG2_BIN_SIZE))
#define TARGET_HITS_PER_PRIME 40.0

uint32 nfs_purge_duplicates(msieve_obj *obj, uint64 mem_size,
				uint64 max_relations,
				uint32 *num_line_num_files)
{

	uint32 i;
	savefile_t *savefile = &obj->savefile;
	char buf[LINE_BUF_SIZE];
	char db_name[LINE_BUF_SIZE];
	uint64 curr_relation;
	uint64 num_relations;
	size_t batch_size;
	void *write_stream;

	uint32 *prime_bins;
	double bin_max;

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

	write_stream = stream_db_init(obj, NULL, 0,
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

	prime_bins = (uint32 *)xcalloc((size_t)1 << (32 - LOG2_BIN_SIZE),
					sizeof(uint32));

	while (!savefile_eof(savefile)) {

		uint8 db_key[24];
		uint8 db_data[24];
		uint32 ab_count;
		uint32 relnum_count;
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
		ab_count = pack_relation(db_key, a, b);
		relnum_count = compress_one_p(db_data, curr_relation);

		stream_db_write_next(write_stream,
				db_key, ab_count,
				db_data, relnum_count);

		/* parse all the factors < 2^32 that the relation 
		   contains and keep a histogram of their counts. 
		   Note that we have no idea whether the relation
		   is a duplicate or even whether the factors are
		   valid, but it doesn't really matter for our
		   purposes */

		while (ptr[0] == ':' || ptr[0] == ',') {
			uint64 p = strtoull(ptr + 1, &ptr, 16);

			if (p >= ((uint64)1 << 32))
				continue;

			prime_bins[p / BIN_SIZE]++;
		}

		/* get the next line */

		savefile_read_line(buf, sizeof(buf), savefile);
	}

	/* free pass 1 stuff */

	savefile_close(savefile);
	stream_db_write_close(write_stream);

	/* pass 2: merge the temporary DBs into a single output one */

	logprintf(obj, "found %" PRIu64 " relations\n", num_relations);

	purge_pass2(obj, write_stream, mem_size, 
			num_relations, num_line_num_files);

	stream_db_free(write_stream);

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
			log((double)BIN_SIZE * i);
	for (i--; i > 2; i--) {
		double bin_min = (double)BIN_SIZE * i /
				log((double)BIN_SIZE * i);
		double hits_per_prime = (double)prime_bins[i] /
						(bin_max - bin_min);
		if (hits_per_prime > TARGET_HITS_PER_PRIME)
			break;
		bin_max = bin_min;
	}

	free(prime_bins);
	return BIN_SIZE * (i + 0.5);
}
