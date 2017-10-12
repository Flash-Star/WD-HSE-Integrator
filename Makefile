CC = gcc
F77 = gfortran -O

GLS_PATH = 
BLAS_PATH = /usr/lib64/atlas
#LFLAGS = -lgsl -L$(BLAS_PATH) -llapack -lcblas -lf77blas -latlas
LFLAGS = -lgsl -L$(BLAS_PATH) -llapack -lgslcblas -lblas -lsatlas
FFLAGS = 
CFLAGS = -Wall -O2

all: SNproj

SNproj: SNproj.o Teos_f.o helm_eos.o
	$(F77) $(LFLAGS) SNproj.o Teos_f.o helm_eos.o -o SNproj

SNproj.o: SNproj.c
	$(CC) $(CFLAGS) -c SNproj.c

Teos_f.o: Teos_f.c composition.h Teos.h
	$(CC) $(CFLAGS) -c Teos_f.c

helm_eos.o: helm_eos.f
	$(F77) $(FFLAGS) -c helm_eos.f

clean:
	rm *.o SNproj
