UNAME_S := $(shell uname -s)

CC = gcc -std=c99 -pedantic -Wall -O3

ifeq ($(UNAME_S), Linux)
    LINKER = gcc -I/opt/intel/mkl/include -L/opt/intel/mkl/lib/intel64 -lmkl_lapack95_lp64 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread
endif
ifeq ($(UNAME_S), Darwin)
    LINKER = gfortran /usr/local/lib/liblapacke.a /usr/local/lib/liblapack.a /usr/local/lib/libblas.a
endif

all: qed

fields.o: fields.c fields.h lattice.h linalg.h complex/complex.h rand/ranlxd.h Makefile
	$(CC) -c $< -o $@

mc.o: mc.c mc.h lattice.h fields.h linalg.h complex/complex.h rand/ranlxd.h rand/gauss.h Makefile
	$(CC) -c $< -o $@

lattice.o: lattice.c lattice.h
	$(CC) -c $< -o $@

linalg.o: linalg.c linalg.h lattice.h complex/complex.h Makefile
	$(CC) -c $< -o $@

ranlxd.o: rand/ranlxd.c rand/ranlxd.h Makefile
	$(CC) -c $< -o $@ -I rand

gauss.o: rand/gauss.c rand/gauss.h Makefile
	$(CC) -c $< -o $@ -I rand

fermion.o: fermion.c fermion.h lattice.h linalg.h complex/complex.h rand/ranlxd.h Makefile
	$(CC) -c $< -o $@

measurement.o: measurement.c measurement.h lattice.h fermion.h fields.h complex/complex.h Makefile
	$(CC) -c $< -o $@

qed.o: qed.c fields.h lattice.h linalg.h mc.h complex/complex.h Makefile rand/ranlxd.h
	$(CC) -c $< -o $@

qed: fields.o qed.o mc.o lattice.o linalg.o  ranlxd.o gauss.o fermion.o measurement.o Makefile
	$(LINKER) qed.o fields.o mc.o lattice.o linalg.o ranlxd.o gauss.o fermion.o measurement.o -o qed -lm

clean:
	rm -f *.o qed test1 test2

