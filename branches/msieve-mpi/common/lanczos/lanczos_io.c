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

#include "lanczos.h"

/*--------------------------------------------------------------------*/
#ifdef HAVE_MPI

typedef struct {
	uint32 col_start;
	uint64 mat_file_offset;
} mat_block_t;

typedef struct {
	uint32 sparse_per_proc;
	uint32 curr_mpi;
	uint32 curr_sparse;
	uint32 target_sparse;
	uint32 curr_col;
	mat_block_t idx_entries[MAX_MPI_PROCS + 1];
} mat_idx_t;

static mat_idx_t * mat_idx_init(uint32 num_sparse) {

	uint32 i;
	mat_idx_t *m = (mat_idx_t *)xcalloc(MAX_MPI_PROCS,
					sizeof(mat_idx_t));

	for (i = 1; i <= MAX_MPI_PROCS; i++)
		m[i-1].sparse_per_proc = num_sparse / i + 100;

	return m;
}

static void mat_idx_update(mat_idx_t *m, FILE *mat_fp,
			uint32 curr_sparse) {

	uint32 i;

	for (i = 0; i < MAX_MPI_PROCS; i++) {
		mat_idx_t *curr_m = m + i;

		if (curr_m->curr_sparse >= curr_m->target_sparse) {
			mat_block_t *curr_block = curr_m->idx_entries +
							curr_m->curr_mpi++;
			curr_block->col_start = curr_m->curr_col;
			curr_block->mat_file_offset = ftello(mat_fp);

			curr_m->target_sparse = curr_m->curr_sparse +
						curr_m->sparse_per_proc;
		}

		curr_m->curr_col++;
		curr_m->curr_sparse += curr_sparse;
	}
}

static void mat_idx_final(msieve_obj *obj, mat_idx_t *m,
			uint32 ncols, uint64 mat_file_size) {

	uint32 i;
	char buf[256];
	FILE *idx_fp;

	sprintf(buf, "%s.mat.idx", obj->savefile.name);
	idx_fp = fopen(buf, "wb");
	if (idx_fp == NULL) {
		logprintf(obj, "error: can't open matrix index file\n");
		exit(-1);
	}

	i = MAX_MPI_PROCS;
	fwrite(&i, sizeof(uint32), (size_t)1, idx_fp);

	for (i = 1; i <= MAX_MPI_PROCS; i++) {
		mat_idx_t *curr_m = m + (i-1);

		curr_m->idx_entries[i].col_start = ncols;
		curr_m->idx_entries[i].mat_file_offset = mat_file_size;

		fwrite(curr_m->idx_entries, sizeof(mat_block_t), 
					(size_t)(i+1), idx_fp);

	}

	fclose(idx_fp);
	free(m);
}

static void find_submatrix_bounds(msieve_obj *obj, uint32 *ncols,
			uint32 *start_col, uint64 *mat_file_offset) {

	mat_block_t mat_block;
	mat_block_t next_mat_block;
	char buf[256];
	FILE *matrix_idx_fp;
	uint32 num_mpi_procs;

	sprintf(buf, "%s.mat.idx", obj->savefile.name);
	matrix_idx_fp = fopen(buf, "rb");
	if (matrix_idx_fp == NULL) {
		logprintf(obj, "error: can't open matrix index file\n");
		exit(-1);
	}

	fread(&num_mpi_procs, sizeof(uint32), (size_t)1, matrix_idx_fp);
	if (num_mpi_procs < obj->mpi_size) {
		logprintf(obj, "error: matrix expects MPI procs <= %u\n",
				num_mpi_procs);
		exit(-1);
	}

	/* get the matrix file offset for MPI process obj->mpi_rank 
	   Also fast-forward to the list specific to obj->mpi_size
	   MPI processes */

	fseek(matrix_idx_fp, 
		(long)((obj->mpi_size *
		        (obj->mpi_size + 1) / 2 - 1 +
			obj->mpi_rank) * sizeof(mat_block_t)), 
		SEEK_CUR);

	fread(&mat_block, sizeof(mat_block_t), 
				(size_t)1, matrix_idx_fp);
	fread(&next_mat_block, sizeof(mat_block_t), 
				(size_t)1, matrix_idx_fp);
	fclose(matrix_idx_fp);

	*start_col = mat_block.col_start;
	*ncols = next_mat_block.col_start - mat_block.col_start;
	*mat_file_offset = mat_block.mat_file_offset;
}

#endif

/*--------------------------------------------------------------------*/
void dump_cycles(msieve_obj *obj, la_col_t *cols, uint32 ncols) {

	uint32 i;
	char buf[256];
	FILE *cycle_fp;

	sprintf(buf, "%s.cyc", obj->savefile.name);
	cycle_fp = fopen(buf, "wb");
	if (cycle_fp == NULL) {
		logprintf(obj, "error: can't open cycle file\n");
		exit(-1);
	}

	fwrite(&ncols, sizeof(uint32), (size_t)1, cycle_fp);

	for (i = 0; i < ncols; i++) {
		la_col_t *c = cols + i;
		uint32 num = c->cycle.num_relations;
		
		fwrite(&num, sizeof(uint32), (size_t)1, cycle_fp);
		fwrite(c->cycle.list, sizeof(uint32), (size_t)num, cycle_fp);
	}
	fclose(cycle_fp);
}

/*--------------------------------------------------------------------*/
void dump_matrix(msieve_obj *obj, 
		uint32 nrows, uint32 num_dense_rows,
		uint32 ncols, la_col_t *cols,
		uint32 sparse_weight) {

	uint32 i;
	uint32 dense_row_words;
	char buf[256];
	FILE *matrix_fp;
#ifdef HAVE_MPI
	mat_idx_t *mpi_idx_data = mat_idx_init(sparse_weight);
#endif

	dump_cycles(obj, cols, ncols);

	sprintf(buf, "%s.mat", obj->savefile.name);
	matrix_fp = fopen(buf, "wb");
	if (matrix_fp == NULL) {
		logprintf(obj, "error: can't open matrix file\n");
		exit(-1);
	}

	fwrite(&nrows, sizeof(uint32), (size_t)1, matrix_fp);
	fwrite(&num_dense_rows, sizeof(uint32), (size_t)1, matrix_fp);
	fwrite(&ncols, sizeof(uint32), (size_t)1, matrix_fp);
	dense_row_words = (num_dense_rows + 31) / 32;

	for (i = 0; i < ncols; i++) {
		la_col_t *c = cols + i;
		uint32 num = c->weight + dense_row_words;

#ifdef HAVE_MPI
		mat_idx_update(mpi_idx_data, matrix_fp, c->weight);
#endif
		fwrite(&c->weight, sizeof(uint32), (size_t)1, matrix_fp);
		fwrite(c->data, sizeof(uint32), (size_t)num, matrix_fp);
	}

#ifdef HAVE_MPI
	mat_idx_final(obj, mpi_idx_data, ncols, ftello(matrix_fp));
#endif
	fclose(matrix_fp);
}

/*--------------------------------------------------------------------*/
#define MAX_CYCLE_SIZE 500

void read_cycles(msieve_obj *obj, 
		uint32 *num_cycles_out, 
		la_col_t **cycle_list_out, 
		uint32 dependency,
		uint32 *colperm) {

	uint32 i;
	uint32 num_cycles;
	uint32 curr_cycle;
	uint32 rel_index[MAX_CYCLE_SIZE];
	char buf[256];
	FILE *cycle_fp;
	FILE *dep_fp = NULL;
	la_col_t *cycle_list = *cycle_list_out;
	uint64 mask = 0;

	if (dependency > 0 && colperm != NULL) {
		logprintf(obj, "error: cannot read dependency with permute\n");
		exit(-1);
	}

	sprintf(buf, "%s.cyc", obj->savefile.name);
	cycle_fp = fopen(buf, "rb");
	if (cycle_fp == NULL) {
		logprintf(obj, "error: read_cycles can't open cycle file\n");
		exit(-1);
	}

	if (dependency) {
		sprintf(buf, "%s.dep", obj->savefile.name);
		dep_fp = fopen(buf, "rb");
		if (dep_fp == NULL) {
			logprintf(obj, "error: read_cycles can't "
					"open dependency file\n");
			exit(-1);
		}
		mask = (uint64)1 << (dependency - 1);
	}

	/* read the number of cycles to expect. If necessary,
	   allocate space for them */

	fread(&num_cycles, sizeof(uint32), (size_t)1, cycle_fp);
	if (cycle_list == NULL) {
		cycle_list = (la_col_t *)xcalloc((size_t)num_cycles, 
						sizeof(la_col_t));
	}

	/* read the relation numbers for each cycle */

	for (i = curr_cycle = 0; i < num_cycles; i++) {

		la_col_t *c;
		uint32 num_relations;

		if (fread(&num_relations, sizeof(uint32), 
					(size_t)1, cycle_fp) != 1)
			break;

		if (num_relations > MAX_CYCLE_SIZE) {
			printf("error: cycle too large; corrupt file?\n");
			exit(-1);
		}

		if (fread(rel_index, sizeof(uint32), (size_t)num_relations, 
					cycle_fp) != num_relations)
			break;

		/* all the relation numbers for this cycle
		   have been read; save them and start the
		   count for the next cycle. If reading in 
		   relations to produce a particular dependency
		   from the linear algebra phase, skip any
		   cycles that will not appear in the dependency */

		if (dependency) {
			uint64 curr_dep;

			if (fread(&curr_dep, sizeof(uint64), 
						(size_t)1, dep_fp) == 0) {
				printf("dependency file corrupt\n");
				exit(-1);
			}
			if (!(curr_dep & mask))
				continue;
		}

		if (colperm != NULL)
			c = cycle_list + colperm[i];
		else
			c = cycle_list + curr_cycle;

		curr_cycle++;
		c->cycle.num_relations = num_relations;
		c->cycle.list = (uint32 *)xmalloc(num_relations * 
						sizeof(uint32));
		memcpy(c->cycle.list, rel_index, 
				num_relations * sizeof(uint32));
	}
	logprintf(obj, "read %u cycles\n", curr_cycle);
	num_cycles = curr_cycle;

	/* check that all cycles have a nonzero number of relations */
	for (i = 0; i < num_cycles; i++) {
		if (cycle_list[i].cycle.num_relations == 0) {
			logprintf(obj, "error: empty cycle encountered\n");
			exit(-1);
		}
	}

	fclose(cycle_fp);
	if (dep_fp) {
		fclose(dep_fp);
	}
	if (num_cycles == 0) {
		free(cycle_list);
		*num_cycles_out = 0;
		*cycle_list_out = NULL;
		return;
	}

	*num_cycles_out = num_cycles;
	*cycle_list_out = (la_col_t *)xrealloc(cycle_list, 
				num_cycles * sizeof(la_col_t));
}
/*--------------------------------------------------------------------*/
static int compare_uint32(const void *x, const void *y) {
	uint32 *xx = (uint32 *)x;
	uint32 *yy = (uint32 *)y;
	if (*xx > *yy)
		return 1;
	if (*xx < *yy)
		return -1;
	return 0;
}

/*--------------------------------------------------------------------*/
uint32 read_matrix(msieve_obj *obj, 
		uint32 *nrows_out, uint32 *num_dense_rows_out,
		uint32 *ncols_out, uint32 *start_col_out, 
		la_col_t **cols_out, uint32 *rowperm, uint32 *colperm) {

	uint32 i, j;
	uint32 dense_row_words;
	uint32 ncols;
	uint32 max_ncols;
	la_col_t *cols;
	char buf[256];
	FILE *matrix_fp;

	if (start_col_out != NULL && colperm != NULL) {
		logprintf(obj, "error: cannot read submatrix with permute\n");
		exit(-1);
	}

	sprintf(buf, "%s.mat", obj->savefile.name);
	matrix_fp = fopen(buf, "rb");
	if (matrix_fp == NULL) {
		logprintf(obj, "error: cannot open matrix file\n");
		exit(-1);
	}

	fread(nrows_out, sizeof(uint32), (size_t)1, matrix_fp);
	fread(num_dense_rows_out, sizeof(uint32), (size_t)1, matrix_fp);
	dense_row_words = (*num_dense_rows_out + 31) / 32;
	fread(&max_ncols, sizeof(uint32), (size_t)1, matrix_fp);
	ncols = max_ncols;

	if (start_col_out != NULL) {
#ifdef HAVE_MPI
		/* read in only a subset of the matrix */

		uint64 mat_file_offset;

		find_submatrix_bounds(obj, &ncols, start_col_out,
					&mat_file_offset);
		fseeko(matrix_fp, mat_file_offset, SEEK_SET);
#else
		*start_col_out = 0;
#endif
	}
	cols = (la_col_t *)xcalloc((size_t)ncols, sizeof(la_col_t));

	for (i = 0; i < ncols; i++) {
		la_col_t *c;
		uint32 num;
		
		if (colperm != NULL)
			c = cols + colperm[i];
		else
			c = cols + i;

		fread(&num, sizeof(uint32), (size_t)1, matrix_fp);
		c->weight = num;
		c->data = (uint32 *)xmalloc((num + dense_row_words) * 
					sizeof(uint32));
		fread(c->data, sizeof(uint32), (size_t)(num + 
				dense_row_words), matrix_fp);

		if (rowperm != NULL) {
			for (j = 0; j < num; j++)
				c->data[j] = rowperm[c->data[j]];
	
			if (num > 1) {
				qsort(c->data, (size_t)num, 
					sizeof(uint32), compare_uint32);
			}
		}
	}
	fclose(matrix_fp);
	*ncols_out = ncols;
	*cols_out = cols;
	return max_ncols;
}

/*--------------------------------------------------------------------*/
void dump_dependencies(msieve_obj *obj, 
			uint64 *deps, uint32 ncols) {

	char buf[256];
	FILE *deps_fp;

	/* we allow up to 64 dependencies, even though the
	   average case will have (64 - POST_LANCZOS_ROWS) */

	sprintf(buf, "%s.dep", obj->savefile.name);
	deps_fp = fopen(buf, "wb");
	if (deps_fp == NULL) {
		logprintf(obj, "error: can't open deps file\n");
		exit(-1);
	}

	fwrite(deps, sizeof(uint64), (size_t)ncols, deps_fp);
	fclose(deps_fp);
}

