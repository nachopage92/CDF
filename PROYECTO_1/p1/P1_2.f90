program EJERCICIOS_EN_FORTRAN_2
	!-------------------------------------------
	!PREAMBULO
	
	! PREGUNTA 2
	integer::n,u0,u1,u2,j
	character(len=20)::datos_p2='./datos/datos_p2.dat'
	
	!-------------------------------------------

!	!PREGUNTA 2
	n=100
	u0=0 ; u1=1
	open(unit=10,file=datos_p2,action='write')
	write(10,*)n
	write(10,*)0,u0
	write(10,*)1,u1
	do j=2,n
		u2=u0+u1
		write(10,*)j,u2
		u0=u1
		u1=u2
	end do
	close(unit=10)
	
end program
