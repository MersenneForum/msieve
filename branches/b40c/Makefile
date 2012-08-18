
NVCC = "$(shell which nvcc)"

GEN_SM20 = -gencode=arch=compute_20,code=\"sm_20,compute_20\" 
GEN_SM13 = -gencode=arch=compute_13,code=\"sm_13,compute_13\" 
GEN_SM10 = -gencode=arch=compute_10,code=\"sm_10,compute_10\" 

CUDA_INC = "$(shell dirname $(NVCC))/../include"
INC = -I$(CUDA_INC) -I. 

NVCCFLAGS = -Xptxas -v -Xcudafe -\# -shared -Xptxas -abi=no

ifeq ($(verbose), 1)
    NVCCFLAGS += -v
endif
ifeq ($(keep), 1)
    NVCCFLAGS += -keep
endif
ifdef maxregisters
    NVCCFLAGS += -maxrregcount $(maxregisters)
endif

DEPS = ./Makefile \
	sort_engine.cu \
	sort_engine.h \
	$(wildcard b40c/util/*.cuh) \
	$(wildcard b40c/util/**/*.cuh) \
	$(wildcard b40c/radix_sort/*.cuh) \
	$(wildcard b40c/radix_sort/**/*.cuh) 

LIBNAME = sort_engine

all: $(LIBNAME)_sm10.dll $(LIBNAME)_sm13.dll $(LIBNAME)_sm20.dll

clean :
	rm -f  *.dll *.lib *.exp

$(LIBNAME)_sm10.dll : $(DEPS)
	$(NVCC) $(GEN_SM10) -o $@ sort_engine.cu $(NVCCFLAGS) $(INC) -O3  

$(LIBNAME)_sm13.dll : $(DEPS)
	$(NVCC) $(GEN_SM13) -o $@ sort_engine.cu $(NVCCFLAGS) $(INC) -O3  

$(LIBNAME)_sm20.dll : $(DEPS)
	$(NVCC) $(GEN_SM20) -o $@ sort_engine.cu $(NVCCFLAGS) $(INC) -O3  

