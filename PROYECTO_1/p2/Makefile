FC = gfortran

MODDIR =$(shell pwd)/modulo
SUBDIR=$(shell pwd)/subrutinas

P2_1: P2_1.f90
	$(FC) -o P2_1 -I$(MODDIR) \
	$(SUBDIR)/eingval.f90 \
	$(SUBDIR)/euler_implicito.f90 \
	$(SUBDIR)/crank_nicolson.f90 \
	$(SUBDIR)/exportar_datos_S1.f90 \
	P2_1.f90
	
P2_2: P2_2.f90
	$(FC) -o P2_2 -I$(MODDIR) \
	$(SUBDIR)/discretizacion.f90 \
	$(SUBDIR)/leap_frog.f90 \
	$(SUBDIR)/newmark.f90 \
	P2_2.f90
	
modulo_pregunta2.mod:
	cd $(MODDIR) && $(MAKE)

clean:
	cd datos/S1 && rm ./*
	cd datos/S2 && rm ./*
	cd graficos/S1 && rm ./*
	cd graficos/S2 && rm ./*
