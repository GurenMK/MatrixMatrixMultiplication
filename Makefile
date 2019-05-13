# On Bridges we will check versus your performance versus Intel MKL library's BLAS. 

CC = cc 
OPT = -O3
CFLAGS = -Wall -fno-signed-zeros -march=native -ffast-math -funroll-loops -std=gnu99 $(OPT)
# -ftree-vectorizer-verbose=2
# Attempted compiler flagfs (no change)
# -fopt-info-vec-missed
# -Ofast
# -fopenmp
# -D_GLIBCXX_PARALLEL
# -frename-registers
# -fno-trapping-math
# -falign-functions=32 (16/24) 
# -falign-loops=32 
# -felide-constructors
# -fstrict-aliasing
# -Wdisabled-optimization

#MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
#LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm
LDLIBS = -lrt  -I$(MKLROOT)/include -Wl,-L$(MKLROOT)/lib/intel64/ -lmkl_intel_lp64 -lmkl_core -lmkl_sequential -lpthread -lm -ldl

targets = benchmark-naive benchmark-blocked benchmark-blas
objects = benchmark.o dgemm-naive.o dgemm-blocked.o dgemm-blas.o

.PHONY : default
default : all

.PHONY : all
all : clean $(targets)

benchmark-naive : benchmark.o dgemm-naive.o 
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-blocked : benchmark.o dgemm-blocked.o
	$(CC) -o $@ $^ $(LDLIBS)
benchmark-blas : benchmark.o dgemm-blas.o
	$(CC) -o $@ $^ $(LDLIBS)

%.o : %.c
	$(CC) -c $(CFLAGS) $<

.PHONY : clean
clean:
	rm -f $(targets) $(objects) *.stdout
