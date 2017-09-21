FC = gfortran

MODDIR =$(shell pwd)/modulo

all: modulo.o run

modulo.o:
	cd modulo/ && $(MAKE)

run:
	@echo $(MODDIR)
	$(FC) -o run -I$(MODDIR) cfd.f90

clean:
	cd modulo/ && $(MAKE) clean
	rm run
