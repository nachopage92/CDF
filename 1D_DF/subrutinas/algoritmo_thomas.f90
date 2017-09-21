!----------------------------------------------------------
! 
! 	U.T.F.S.M. - IPM-468 - CFD_2_2017 
!	Ignacio Apablaza
! 
!	Implementacion del algoritmo de Thomas
!	
! 	descripcion:
!		breve descripcion
! 
!	input:
!		nombrar input
! 
! 	output:
!		nombrar output
! 
!----------------------------------------------------------

subroutine thomas(d,a,c,x,b,n)

	implicit none
	
	!input
	integer										:: n
	double precision,dimension(n),intent(in)	:: a,b
	double precision,dimension(1:n-1),intent(in):: c
	double precision,dimension(2:n),intent(in)	:: d
	
	!output
	double precision,dimension(n),intent(out)	:: x
	
	!variables locales
	integer										:: i
	double precision,dimension(n)				:: y,alfa
	double precision,dimension(2:n)				:: beta
	
	!calculo de los coeficientes de la matriz L y U 
	!(alfa, beta y gama) 
	alfa(1) = a(1)
	beta(2) = d(2)/alfa(1)	
	do i=2,n
		beta(i)=d(i)/alfa(i-1)
		alfa(i)=a(i)-b(i)*c(i-1)
	end do

	!resolver sistema Ly=b
	y(1)=b(1)
	do i=2,n
		y(i) = b(i) - beta(i) * y(i-1)
	end do
	
	!resolver sistema Ux=y
	x(n)=y(n)/alfa(n)
	do i=n-1,1,-1
		x(i) = ( y(i)-c(i)*x(i+1) ) / alfa(i)
	end do
	
end subroutine
