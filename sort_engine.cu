#include <stdio.h>
#include <b40c/util/error_utils.cuh>
#include <b40c/util/multi_buffer.cuh>
#include <b40c/radix_sort/enactor.cuh>

#include "sort_engine.h"

typedef unsigned int uint32;

#if defined(_WIN32) || defined (_WIN64)
	#define SORT_ENGINE_DECL __declspec(dllexport)
	typedef unsigned __int64 uint64;
#else
	#define SORT_ENGINE_DECL __attribute__((visibility("default")))
	typedef unsigned long long uint64;
#endif

typedef uint64 KEY_TYPE;
typedef uint32 DATA_TYPE;

using namespace b40c;

typedef struct
{
	radix_sort::Enactor enactor;
} sort_engine;

extern "C"
{

SORT_ENGINE_DECL void * 
sort_engine_init(void)
{
	return new sort_engine;
}

SORT_ENGINE_DECL void 
sort_engine_free(void * e)
{
	delete (sort_engine *)e;
}

SORT_ENGINE_DECL void 
sort_engine_run(void * e, sort_data_t * data)
{
	sort_engine *engine = (sort_engine *)e;
	bool need_swap;

	// arrays are assumed packed together; check
	// they would all start on a power-of-two boundary

	if (data->num_arrays > 1 && data->num_elements % 16) {
		printf("sort_engine: invalid array size\n");
		exit(-1);
	}

	if (data->sort_keys_only) {
		for (size_t i = 0; i < data->num_arrays; i++) {

			cudaError_t status;
			util::DoubleBuffer<KEY_TYPE> ptrs;

			ptrs.d_keys[0] = (KEY_TYPE *)data->keys_in +
						i * data->num_elements;
			ptrs.d_keys[1] = (KEY_TYPE *)data->keys_in_scratch +
						i * data->num_elements;

			if (data->num_elements < 32000) {
				status = engine->enactor.Sort<radix_sort::SMALL_PROBLEM, 
				       		8 * sizeof(KEY_TYPE), 0>(ptrs, data->num_elements);
			}
			else {
				status = engine->enactor.Sort<radix_sort::LARGE_PROBLEM, 
				       		8 * sizeof(KEY_TYPE), 0>(ptrs, data->num_elements);
			}

			need_swap = (ptrs.selector > 0);
			if (status != CUDA_SUCCESS) {
				util::B40CPerror(status, "sort engine: ", __FILE__, __LINE__);
				exit(-1);
			}
		}
	}
	else {
		for (size_t i = 0; i < data->num_arrays; i++) {

			cudaError_t status;
			util::DoubleBuffer<KEY_TYPE, DATA_TYPE> ptrs;

			ptrs.d_keys[0] = (KEY_TYPE *)data->keys_in +
						i * data->num_elements;
			ptrs.d_keys[1] = (KEY_TYPE *)data->keys_in_scratch +
						i * data->num_elements;
			ptrs.d_values[0] = (DATA_TYPE *)data->data_in +
						i * data->num_elements;
			ptrs.d_values[1] = (DATA_TYPE *)data->data_in_scratch +
						i * data->num_elements;

			if (data->num_elements < 32000) {
				status = engine->enactor.Sort<radix_sort::SMALL_PROBLEM, 
				       		8 * sizeof(KEY_TYPE), 0>(ptrs, data->num_elements);
			}
			else {
				status = engine->enactor.Sort<radix_sort::LARGE_PROBLEM, 
				       		8 * sizeof(KEY_TYPE), 0>(ptrs, data->num_elements);
			}

			need_swap = (ptrs.selector > 0);
			if (status != CUDA_SUCCESS) {
				util::B40CPerror(status, "sort engine: ", __FILE__, __LINE__);
				exit(-1);
			}
		}
	}

	if (need_swap == true) {
		std::swap(data->keys_in, data->keys_in_scratch);
		std::swap(data->data_in, data->data_in_scratch);
	}
}

} // extern "C"
