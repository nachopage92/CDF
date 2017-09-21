!----------------------------------------------------------
! 
! 	U.T.F.S.M. - IPM-468 - CFD_2_2017 
!	Ignacio Apablaza
! 
!	Aproximacion numerica de la ecuacion difusion para 1D
!					dU/dt = D d2u/dx2
! 
!	utilizando un esquema explicito:
!		DF hacia adelante (tiempo)
!		DF centrado (espacio)
! 
!	tipo de discretizacion
!		Malla regular
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
	double precision,dimension(n)	:: x,T_1,T_2
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
		
	do k = 1,10	!(pasos de tiempo)
	
		!condiciones de contorno
		! nodo_1: dirichlet (esencial)
		!AQUI ALGO
		! nodo_n: neumann (natural)
		!AQUI ALGO
	
		!nodo_1:
		T_2(1) = 10d0
		
		!nodos interiores:
		do i = 2 , n-1
			T_2(i) = (d_t*D/d_x**2d0) * &
			& (T_1(i-1)-2d0*T_1(i)+T_1(i+1)) + T_1(i)
		end do
		
		!nodo_n:
		T_2(n) = 0d0
	
		!exportar datos
		write(1,*) 'paso de tiempo : ', k
		write(1,*) T_2
		write(1,*) 
		write(1,*)
		
		T_1=T_2
	
	end do
	 
	

end program

