FC = gfortran

MODDIR =$(shell pwd)/modulo
OBJDIR =$(shell pwd)/subrutinas

all: modulo.o algoritmo_thomas.o run

modulo.o:
	cd modulo/ && $(MAKE)

algoritmo_thomas.o:
	cd subrutinas/ && $(MAKE)

run:
	@echo $(MODDIR)
	$(FC) -o run -I$(MODDIR) $(OBJDIR)/algoritmo_thomas.o cfd_2.f90

clean:
	cd modulo/ && $(MAKE) clean
	rm run
