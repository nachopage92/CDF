FC = gfortran

MODDIR =$(shell pwd)/modulo
SUBDIR=$(shell pwd)/subrutinas

all: modulo_pregunta2.mod P2

S1: P2_S1.f90
	$(FC) -o S1 -I$(MODDIR) \
	$(SUBDIR)/euler_implicito.f90 \
	$(SUBDIR)/crank_nicolson.f90 \
	P2_S1.f90
	
S2: P2_S2.f90
	$(FC) -o S2 -I$(MODDIR) \
	$(SUBDIR)/euler_implicito.f90 \
	$(SUBDIR)/crank_nicolson.f90 \
	P2_S2.f90
	
S3: P2_S3.f90
	$(FC) -o S3 -I$(MODDIR) \
	$(SUBDIR)/discretizacion.f90 \
	$(SUBDIR)/leap_frog.f90 \
	P2_S3.f90
	
modulo_pregunta2.mod:
	cd $(MODDIR) && $(MAKE)

clean:
	rm P2
	cd $(MODDIR) && $(MAKE) clean
