!----------------------------------------------------------
! 
! 	U.T.F.S.M. - IPM-468 - CFD_2_2017 
!	Ignacio Apablaza
! 
!	Aproximacion numerica de la ecuacion difusion para 1D
!					dU/dt = D d2u/dx2
! 
!	utilizando un esquema implicito:
!		DF hacia adelante (tiempo)
!		DF centrado (espacio)
! 
!	tipo de discretizacion
!		Malla regular
! 
! 	solver
!		Algoritmo de Thomas
! 
!----------------------------------------------------------

program CDF_1D

	use cfd_modulo 
	
	implicit none
	
	!parametros del problema
	integer,parameter			:: n = 11	!n de nodos
	double precision,parameter	:: L = 10 , D = 0.25	
	
	!variables del programa
	integer							:: i,k
	double precision				:: d_x,d_t,n_fourier
	double precision,dimension(n)	:: x,d2,b,T_1,T_2
	double precision,dimension(2:n)	:: d1
	double precision,dimension(1:n-1):: d3
	integer,dimension(n)			:: indx
	
	!discretizacion del dominio (malla regular)
	d_x = L / (n-1)
	x(:) = (/ (d_x*(i-1) , i=1,n) /)
	indx(:) = (/ (i,i=1,n) /)
	
	!se define el paso de tiempo
	d_t = 1d0
	
	!mostrar numero de fourier
	n_fourier = D*d_t/d_x
	if (n_fourier .ge. 0.5) then
		write(*,*) 'revisar condicion de estabilidad'
		write(*,*) 'numero de fourier :  ' , n_fourier 
		return
	end if
	
	!condicion inicial
	T_1(:) = (/ (0 , i=1,n) /)
	
	!se escriben los vectores que definen a la 
	!matriz A del sistema (Ax=b)
		
	d1(2:n-1) = (/ (-D*d_t/d_x , i=1,n-2) /)
	d1(n) = 0d0
	
	d3(1) = 0d0
	d3(2:n-1) = (/ (-D*d_t/d_x , i=1,n-2) /)
	
	d2(1) = 1d0
	d2(2:n-1) = (/ (1d0 + 2d0*D*d_t/d_x , i=1,n-2) /)
	d2(n) = 1d0
		
	do k = 1,10	!(pasos de tiempo)
	
		!condiciones de contorno
		! nodo_1: dirichlet (esencial)
		!AQUI ALGO
		! nodo_n: neumann (natural)
		!AQUI ALGO
		
		b(1) = 10d0
		b(2:n-1) = T_1(2:n-1)
		b(n) = 0d0
	
		call thomas(d1,d2,d3,T_2,b,n)
	
		!exportar datos
		write(1,*) 'paso de tiempo : ', k
		write(1,*) T_2
		write(1,*) 
		write(1,*)
		
		T_1=T_2
	
	end do
	 
	

end program

	
