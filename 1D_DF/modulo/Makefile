FC=gfortran
CFLAG=-c

all: modulo.o

modulo.o: modulo.f90
	$(FC) $(CFLAG) modulo.f90

clean:	
	find . -type f -not -name 'Makefile' \
		 -not -name 'modulo.f90' -delete
