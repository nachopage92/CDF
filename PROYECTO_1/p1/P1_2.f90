program EJERCICIOS_EN_FORTRAN_2
	!-------------------------------------------
	!PREAMBULO
	integer::n,u0,u1,u2,j
	real(kind=8)::u0_r,u1_r,u2_r
	character(len=22)::datos_p2_1='./datos/datos_p2_1.dat'
	character(len=22)::datos_p2_2='./datos/datos_p2_2.dat'
	!-------------------------------------------

	!se calcula la serie para número enteros
	n=100
	u0=0 ; u1=1
	open(unit=10,file=datos_p2_1,action='write')
		write(10,*)0,u0
		write(10,*)1,u1
		do j=2,n
			u2=u0+u1
			write(10,*)j,u2
			u0=u1
			u1=u2
		end do
	close(unit=10)
	
	!se calcula la serie para número reales d. precision
	n=1500
	u0_r = 0._8 ; u1_r = 1._8
	open(unit=10,file=datos_p2_2,action='write')
		write(10,*)0,u0_r
		write(10,*)1,u1_r
		do j=2,n
			u2_r = u0_r + u1_r
			write(10,*) j , u2_r
			u0_r = u1_r
			u1_r = u2_r
		end do
	close(unit=10)
	
	!-------------------------------------------
	
	! se grafica la curva en el directorio
	! ./graficos/graficos_p1.dat
	call system ( ' cd ./gnuplot && gnuplot plot_p1_2.txt ' )
	
	!-------------------------------------------
	
end program
