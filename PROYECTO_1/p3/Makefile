FC = gfortran

MODDIR =$(shell pwd)/modulo
SUBDIR=$(shell pwd)/subrutinas

all: P3.f90
	$(FC) -o P3 \
	$(SUBDIR)/rk4.f90 \
	P3.f90
