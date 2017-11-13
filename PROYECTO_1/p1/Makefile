FC = gfortran
SUBDIR=./subrutinas

all: P1_1 P1_2 P1_3
	$(FC) -o P1_1 P1_1.f90
	$(FC) -o P1_2 P1_2.f90
	$(FC) -o P1_3 $(SUBDIR)/matrix_mult.f90 P1_3.f90
	./P1_1 && ./P1_2 && \
	cd gnuplot/ && gnuplot plot_p1.txt
	./P1_3
	
P1_1: P1_1.f90
	$(FC) -o P1_1 P1_1.f90 

P1_2: P1_2.f90
	$(FC) -o P1_2 P1_2.f90
	
P1_3: P1_3.f90
	$(FC) -o P1_3 $(SUBDIR)/matrix_mult.f90 P1_3.f90 
	
clean:
	cd datos/ && rm ./*
	cd graficos/ && rm ./*
		

