#
# Computers with Red Hat Enterprise Linux 5 in the computer room 648, KTH Forum, Kista
#
CC = g++
MPCC =  mpicc -cc=gcc

TARGETS = serial pthreads openmp mpi

all:	$(TARGETS)

serial:
	$(CC) -O3 -o $@ -lm serial.cpp common.cpp
pthreads:
	$(CC) -O3 -o $@ pthreads.cpp common.cpp -lm -lpthread
openmp:
	$(CC) -O3 -o $@ openmp.cpp common.cpp -lm -fopenmp
mpi:
	$(CC) -O3 -o $@ mpi.cpp common.cpp -lmpicxx -lmpi

clean:
	rm -f *.o $(TARGETS)
