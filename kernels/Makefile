
# Compiler options

PTX = .

NVCCFLAGS = -ptx -arch=compute_20 -code=sm_20 

T128 = -D NUM_THREADS=128
T256 = -D NUM_THREADS=256
INDEX = -D VALUE_TYPE_INDEX
KEY = -D VALUE_TYPE_NONE
SINGLE = -D VALUE_TYPE_SINGLE
MULTI = -D VALUE_TYPE_MULTI
LIST = -D SCATTER_TRANSACTION_LIST
SIMPLE = -D SCATTER_SIMPLE

SORTFLAGS = $(NVCCFLAGS) -D VALUES_PER_THREAD=8

# output files

COUNTS = \
	$(PTX)/count.ptx

HISTS = \
	$(PTX)/hist_simple.ptx
	
SORTS = \
	$(PTX)/sort_128_8_index_simple.ptx \
	$(PTX)/sort_128_8_key_simple.ptx \
	$(PTX)/sort_128_8_single_simple.ptx \
	$(PTX)/sort_128_8_multi_simple.ptx \
	$(PTX)/sort_256_8_index_simple.ptx \
	$(PTX)/sort_256_8_key_simple.ptx \
	$(PTX)/sort_256_8_single_simple.ptx \
	$(PTX)/sort_256_8_multi_simple.ptx

# Define dependencies for count, hist, and sort
COUNT_DEPS = count.cu countgen.cu common.cu

HIST_DEPS = histgen.cu histogram.cu common.cu params.cu \
	hist1.cu hist2.cu hist3.cu hist3_64.cu

SORT_DEPS = common.cu params.cu sort.cu sortcommon.cu sortgen.cu \
	sortgeneric.cu sortscan1.cu sortscan2.cu sortscan3.cu sortstore.cu
	 
# Generate all kernels.
all : $(SORTS) $(HISTS) $(COUNTS)

clean :
	rm -f $(PTX)/*.ptx	

	
# Build count ptx
$(PTX)/count.ptx : $(COUNT_DEPS)
	nvcc $(NVCCFLAGS) -o $@ countgen.cu
	
# Build hist ptx
$(PTX)/hist_simple.ptx : $(HIST_DEPS)
	nvcc $(NVCCFLAGS) -o $@ histgen.cu

# Build sort ptx
$(PTX)/sort_128_8_index_simple.ptx : $(SORT_DEPS)
	nvcc $(SORTFLAGS) $(T128) $(INDEX) $(SIMPLE) -o $@ sortgen.cu
	
$(PTX)/sort_128_8_key_simple.ptx : $(SORT_DEPS)
	nvcc $(SORTFLAGS) $(T128) $(KEY) $(SIMPLE) -o $@ sortgen.cu
	
$(PTX)/sort_128_8_single_simple.ptx : $(SORT_DEPS)
	nvcc $(SORTFLAGS) $(T128) $(SINGLE) $(SIMPLE) -o $@ sortgen.cu
	
$(PTX)/sort_128_8_multi_simple.ptx : $(SORT_DEPS)
	nvcc $(SORTFLAGS) $(T128) $(MULTI) $(SIMPLE) -o $@ sortgen.cu
	
$(PTX)/sort_256_8_index_simple.ptx : $(SORT_DEPS)
	nvcc $(SORTFLAGS) $(T256) $(INDEX) $(SIMPLE) -o $@ sortgen.cu
	
$(PTX)/sort_256_8_key_simple.ptx : $(SORT_DEPS)
	nvcc $(SORTFLAGS) $(T256) $(KEY) $(SIMPLE) -o $@ sortgen.cu
	
$(PTX)/sort_256_8_single_simple.ptx : $(SORT_DEPS)
	nvcc $(SORTFLAGS) $(T256) $(SINGLE) $(SIMPLE) -o $@ sortgen.cu
	
$(PTX)/sort_256_8_multi_simple.ptx : $(SORT_DEPS)
	nvcc $(SORTFLAGS) $(T256) $(MULTI) $(SIMPLE) -o $@ sortgen.cu
	
