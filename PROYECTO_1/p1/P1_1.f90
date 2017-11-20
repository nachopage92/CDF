program EJERCICIOS_EN_FORTRAN_1
	!-------------------------------------------
	!PREAMBULO
	integer::i,inf,contador
	real(kind=4)::A !simple precision
	real(kind=8)::B,error !doble precision
	character(len=22)::datos_p1='./datos/datos_p1.dat'
	character(len=22)::datos_pe='./datos/datos_p1_error.dat'
	!-------------------------------------------
 	
 	!simple y doble precisión
	A=1._4
	B=1._8
	inf=4000000
	contador=1
	
	! los datos generados se exportan al directorio
	! ./datos/datos_p1.dat
	open(unit=10,file=datos_p1,action='write')
	do i=2,inf
		!contador
		contador = contador + 1
		!desarrollo de la serie
		A = A + ( 1._4 / real(i,kind=4) )
		B = B + ( 1._8 / real(i,kind=8) )
		error = (B-real(A,kind=8))/B
		!cada 1000 operaciones exportar datos
		if ( contador .eq. 10000 ) then
			if ( i .gt. 250000 ) then
				write(10,*) i, A, B, error
			end if
			contador=0
		end if
	end do 
	close(unit=10)
	
	!-------------------------------------------
	
	! se grafica la curva en el directorio
	! ./graficos/graficos_p1.dat
	call system ( ' cd ./gnuplot && gnuplot plot_p1_1.txt ' )
	
	!-------------------------------------------
	
end program

