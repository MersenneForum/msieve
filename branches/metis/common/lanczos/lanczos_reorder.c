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

#ifdef HAVE_METIS

#include "lanczos.h"
#include "metis/metis.h"

typedef struct {
	uint16 num_edges;
	uint16 curr_off;
	uint32 edge_offset;
	uint32 orig_index;
	uint32 perm;
} vertex_t;

typedef struct {
	uint32 num_rows;
	uint32 num_cols;
	uint32 num_rows_orig;
	uint32 num_cols_orig;

	vertex_t *vertex;
	uint32 *edges;
} graph_t;


#define MIN_ROW_IDX 50000
#define MAX_EDGES 20
#define INVALID_INDEX ((uint32)(-1))

/*--------------------------------------------------------------------*/
#define MAX_ROW_ENTRIES 1000

static void graph_init(msieve_obj *obj, graph_t *graph) {

	uint32 i, j, k;
	uint32 nrows;
	uint32 dense_row_words;
	uint32 ncols;
	char buf[256];
	FILE *matrix_fp;

	uint32 num_edges;
	uint32 num_edges_alloc;
	uint32 *edges;
	vertex_t *vertex;
	uint32 num_vertex;

	sprintf(buf, "%s.mat", obj->savefile.name);
	matrix_fp = fopen(buf, "rb");
	if (matrix_fp == NULL) {
		logprintf(obj, "error: can't open matrix file\n");
		exit(-1);
	}

	fread(&nrows, sizeof(uint32), (size_t)1, matrix_fp);
	fread(&dense_row_words, sizeof(uint32), (size_t)1, matrix_fp);
	fread(&ncols, sizeof(uint32), (size_t)1, matrix_fp);
	dense_row_words = (dense_row_words + 31) / 32;

	vertex = (vertex_t *)xcalloc((size_t)ncols + nrows - MIN_ROW_IDX, 
					sizeof(vertex_t));

	num_vertex = 0;
	num_edges = 0;
	num_edges_alloc = 100000;
	edges = (uint32 *)xmalloc(num_edges_alloc * sizeof(uint32));

	for (i = 0; i < ncols; i++) {

		uint32 num_row_idx;
		uint32 row_entries[MAX_ROW_ENTRIES];
		vertex_t *v;

		fread(&num_row_idx, sizeof(uint32), 
				(size_t)1, matrix_fp);
		fread(row_entries, sizeof(uint32), 
				(size_t)(num_row_idx + dense_row_words), 
				matrix_fp);

		if (num_edges + num_row_idx >= num_edges_alloc) {
			num_edges_alloc *= 1.4;
			edges = (uint32 *)xrealloc(edges, num_edges_alloc *
						sizeof(uint32));
		}

		for (j = num_row_idx, k = 0; j; j--) {
			uint32 row = row_entries[j-1];

			if (row < MIN_ROW_IDX)
				break;

			v = vertex + ncols + row - MIN_ROW_IDX;
			if (v->num_edges < MAX_EDGES) {
				v->num_edges++;
				edges[num_edges + k] = row - MIN_ROW_IDX;
				k++;
			}
		}
		if (k == 0 || k >= MAX_EDGES)
			continue;

		v = vertex + num_vertex++;
		v->num_edges = k;
		v->orig_index = i;
		v->edge_offset = num_edges;
		num_edges += k;
	}
	fclose(matrix_fp);

	for (i = num_edges = 0; i < num_vertex; i++) {
		vertex_t *v = vertex + i;
		uint32 curr_num_edges = v->num_edges;
		uint32 *row_list = edges + v->edge_offset;

		for (j = k = 0; j < curr_num_edges; j++) {
			uint32 row = row_list[j];
			vertex_t *row_v = vertex + ncols + row;

			if (row_v->num_edges < MAX_EDGES) {
				edges[num_edges + k] = row;
				k++;
			}
		}
		v->num_edges = k;
		v->edge_offset = num_edges;
		num_edges += k;
	}

	for (i = 0; i < nrows - MIN_ROW_IDX; i++) {
		vertex[ncols + i].num_edges = 0;
	}

	for (i = j = 0; i < num_vertex; i++) {
		vertex_t *v = vertex + i;
		if (v->num_edges > 0) {
			uint32 curr_num_edges = v->num_edges;
			uint32 *row_list = edges + v->edge_offset;

			for (k = 0; k < curr_num_edges; k++) {
				vertex_t *row_v = vertex + ncols + row_list[k];
				row_v->num_edges++;
			}
			vertex[j++] = *v;
		}
	}
	num_vertex = j;

	for (i = 0, j = num_edges, k = num_vertex; 
				i < nrows - MIN_ROW_IDX; i++) {
		vertex_t *v = vertex + ncols + i;
		if (v->num_edges > 0 && v->num_edges < MAX_EDGES) {
			v->orig_index = i + MIN_ROW_IDX;
			v->edge_offset = j;
			v->perm = k++;
			j += v->num_edges;
		}
	}

	for (i = 0; i < num_vertex; i++) {
		vertex_t *v = vertex + i;
		uint32 curr_num_edges = v->num_edges;
		uint32 *row_list = edges + v->edge_offset;

		for (j = 0; j < curr_num_edges; j++) {
			vertex_t *row_v = vertex + ncols + row_list[j];
			row_list[j] = row_v->perm;
		}
	}

	for (i = j = 0; i < nrows - MIN_ROW_IDX; i++) {
		vertex_t *v = vertex + ncols + i;
		if (v->num_edges > 0 && v->num_edges < MAX_EDGES) {
			vertex[num_vertex + j] = *v;
			j++;
		}
	}

	graph->num_cols_orig = ncols;
	graph->num_rows_orig = nrows;
	logprintf(obj, "matrix is %u x %u\n", nrows, ncols);
	graph->num_cols = ncols = num_vertex;
	graph->num_rows = nrows = j;
	logprintf(obj, "sparse core is %u x %u\n", nrows, ncols);
	logprintf(obj, "graph has %u edges\n", num_edges);
	logprintf(obj, "memory use: %.2lf MB\n", 
			((nrows + ncols) * sizeof(vertex_t) +
			 2 * num_edges * sizeof(uint32)) / 1048576.0);

	vertex = (vertex_t *)xrealloc(vertex, (nrows + ncols) *
					sizeof(vertex_t));
	edges = (uint32 *)xrealloc(edges, 2 * num_edges * sizeof(uint32));

	for (i = 0; i < ncols; i++) {
		vertex_t *v = vertex + i;
		uint32 curr_num_edges = v->num_edges;
		uint32 *row_list = edges + v->edge_offset;

		for (j = 0; j < curr_num_edges; j++) {
			vertex_t *row_v = vertex + row_list[j];
			uint32 *row_edges = edges + row_v->edge_offset;
			row_edges[row_v->curr_off++] = i;
		}
	}

	graph->edges = edges;
	graph->vertex = vertex;
}

/*--------------------------------------------------------------------*/
static void graph_partition(msieve_obj *obj, graph_t *graph) {

	uint32 i;
	vertex_t *vertex = graph->vertex;
	uint32 num_vertex = graph->num_rows + graph->num_cols;
	uint32 *edge_offsets;
	uint32 *partition;
	uint32 num_partitions = graph->num_cols / 100000;
	uint32 num_constraints = 2;
	float balance[2] = {1.05, 1.05};
	uint32 constraint_flag = 0;
	uint32 *constraints;
	uint32 numbering = 0;
	uint32 edge_cut;
	uint32 use_kway = (num_partitions > 8);

	edge_offsets = (uint32 *)xmalloc((num_vertex + 1) *
					sizeof(uint32));
	partition = (uint32 *)xmalloc(num_vertex *
					sizeof(uint32));
	constraints = (uint32 *)xmalloc(num_constraints * 
					num_vertex * sizeof(uint32));

	for (i = 0; i < num_vertex; i++) {
		edge_offsets[i] = vertex[i].edge_offset;

		if (i < graph->num_cols) {
			constraints[2*i+0] = 0;
			constraints[2*i+1] = 1;
		}
		else {
			constraints[2*i+0] = 1;
			constraints[2*i+1] = 0;
		}
	}
	edge_offsets[i] = edge_offsets[i-1] + vertex[i-1].num_edges;

	logprintf(obj, "computing a %u-way partitioning\n", num_partitions);

	if (use_kway) {
		uint32 options[5] = {1, 
				MTYPE_SBHEM_ONENORM, 
				ITYPE_RANDOM, 
				RTYPE_KWAYRANDOM, 
				0};

		METIS_mCPartGraphKway(
				(int *)&num_vertex,
				(int *)&num_constraints,
				(idxtype *)edge_offsets,
				(idxtype *)graph->edges,
				(idxtype *)constraints,
				NULL,
				(int *)&constraint_flag,
				(int *)&numbering,
				(int *)&num_partitions,
				balance,
				(int *)options,
				(int *)&edge_cut,
				(idxtype *)partition
				);
	}
	else {
		uint32 options[5] = {1, 
				MTYPE_SBHEM_ONENORM, 
				ITYPE_RANDOM, 
				RTYPE_FM, 
				0};

		METIS_mCPartGraphRecursive(
				(int *)&num_vertex,
				(int *)&num_constraints,
				(idxtype *)edge_offsets,
				(idxtype *)graph->edges,
				(idxtype *)constraints,
				NULL,
				(int *)&constraint_flag,
				(int *)&numbering,
				(int *)&num_partitions,
				(int *)options,
				(int *)&edge_cut,
				(idxtype *)partition
				);
	}

	logprintf(obj, "partition complete, edge cut = %u\n", edge_cut);

	for (i = 0; i < num_vertex; i++)
		vertex[i].perm = partition[i];

	free(partition);
	free(edge_offsets);
	free(constraints);
}

/*--------------------------------------------------------------------*/
static int compare_vertex(const void *x, const void *y) {
	vertex_t *xx = (vertex_t *)x;
	vertex_t *yy = (vertex_t *)y;

	return (int)xx->perm - (int)yy->perm;
}

static void graph_free(graph_t *graph, 
			uint32 **rowperm_out,
			uint32 **colperm_out) {

	uint32 i, j;
	uint32 nrows = graph->num_rows;
	uint32 ncols = graph->num_cols;
	uint32 nrows_orig = graph->num_rows_orig;
	uint32 ncols_orig = graph->num_cols_orig;
	uint32 rowgap = nrows_orig - nrows;
	uint32 colgap = ncols_orig - ncols;
	uint32 *rowperm;
	uint32 *colperm;

	free(graph->edges);
	rowperm = (uint32 *)xmalloc(nrows_orig * sizeof(uint32));
	colperm = (uint32 *)xmalloc(ncols_orig * sizeof(uint32));

	for (i = 0; i < nrows_orig; i++)
		rowperm[i] = INVALID_INDEX;

	qsort(graph->vertex + ncols, nrows, 
			sizeof(vertex_t), compare_vertex);

	for (i = 0; i < nrows; i++) {
		vertex_t *v = graph->vertex + ncols + i;
		rowperm[v->orig_index] = i + rowgap;
	}

	for (i = j = 0; i < nrows_orig; i++) {
		if (rowperm[i] == INVALID_INDEX)
			rowperm[i] = j++;
	}


	for (i = 0; i < ncols_orig; i++)
		colperm[i] = INVALID_INDEX;

	qsort(graph->vertex, ncols, 
			sizeof(vertex_t), compare_vertex);

	for (i = 0; i < ncols; i++) {
		vertex_t *v = graph->vertex + i;
		colperm[v->orig_index] = i + colgap;
	}

	for (i = j = 0; i < ncols_orig; i++) {
		if (colperm[i] == INVALID_INDEX)
			colperm[i] = j++;
	}

	free(graph->vertex);
	*rowperm_out = rowperm;
	*colperm_out = colperm;
}

/*--------------------------------------------------------------------*/
void reorder_matrix(msieve_obj *obj, 
		    uint32 **rowperm, 
		    uint32 **colperm) {

	graph_t graph;

	logprintf(obj, "permuting matrix for faster multiplies\n");

	graph_init(obj, &graph);

	graph_partition(obj, &graph);

	graph_free(&graph, rowperm, colperm);
}

#endif /* HAVE_METIS */
