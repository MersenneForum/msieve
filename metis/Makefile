
CC = gcc
CFLAGS = -O3 -march=native -fomit-frame-pointer -D_FILE_OFFSET_BITS=64

INC_DIRS = -I. -Ilibmetis -IGKlib

LIB = libmetis.a

METIS_SRCS = \
	GKlib/b64.c \
	GKlib/blas.c \
	GKlib/dfkvkselect.c \
	GKlib/dlmalloc.c \
	GKlib/error.c \
	GKlib/htable.c \
	GKlib/memory.c \
	GKlib/omp.c \
	GKlib/pqueue.c \
	GKlib/sort.c \
	GKlib/string.c \
	GKlib/timers.c \
	GKlib/util.c \
	libmetis/balance.c \
	libmetis/bucketsort.c \
	libmetis/ccgraph.c \
	libmetis/coarsen.c \
	libmetis/compress.c \
	libmetis/debug.c \
	libmetis/estmem.c \
	libmetis/fm.c \
	libmetis/fortran.c \
	libmetis/frename.c \
	libmetis/graph.c \
	libmetis/initpart.c \
	libmetis/kmetis.c \
	libmetis/kvmetis.c \
	libmetis/kwayfm.c \
	libmetis/kwayrefine.c \
	libmetis/kwayvolfm.c \
	libmetis/match.c \
	libmetis/mbalance.c \
	libmetis/mbalance2.c \
	libmetis/mcoarsen.c \
	libmetis/memory.c \
	libmetis/mfm.c \
	libmetis/mfm2.c \
	libmetis/mincover.c \
	libmetis/minitpart.c \
	libmetis/minitpart2.c \
	libmetis/mkmetis.c \
	libmetis/mkwayfmh.c \
	libmetis/mkwayrefine.c \
	libmetis/mmatch.c \
	libmetis/mpmetis.c \
	libmetis/mrefine.c \
	libmetis/mrefine2.c \
	libmetis/mrkmetis.c \
	libmetis/mutil.c \
	libmetis/myqsort.c \
	libmetis/pmetis.c \
	libmetis/pqueue.c \
	libmetis/refine.c \
	libmetis/rkmetis.c \
	libmetis/separator.c \
	libmetis/sfm.c \
	libmetis/srefine.c \
	libmetis/stat.c \
	libmetis/subdomains.c \
	libmetis/timing.c \
	libmetis/util.c

METIS_OBJS = $(METIS_SRCS:.c=.o)

METIS_HDRS = \
	GKlib/GKlib.h \
	GKlib/gk_arch.h \
	GKlib/gk_blas.h \
	GKlib/gk_defs.h \
	GKlib/gk_dlmalloc.h \
	GKlib/gk_externs.h \
	GKlib/gk_getopt.h \
	GKlib/gk_macros.h \
	GKlib/gk_memory.h \
	GKlib/gk_proto.h \
	GKlib/gk_sort.h \
	GKlib/gk_struct.h \
	GKlib/gk_types.h \
	libmetis/defs.h \
	libmetis/macros.h \
	libmetis/metislib.h \
	libmetis/proto.h \
	libmetis/rename.h \
	libmetis/stdheaders.h \
	libmetis/struct.h \
	metis.h

$(LIB): $(METIS_OBJS)
	ar r $(LIB) $(METIS_OBJS)
	ranlib $(LIB)

clean:
	rm -rf $(LIB) $(METIS_OBJS)

%.o: %.c $(METIS_HDRS)
	$(CC) $(CFLAGS) $(INC_DIRS) -c -o $@ $<
