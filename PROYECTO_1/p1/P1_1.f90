program EJERCICIOS_EN_FORTRAN_1
	!-------------------------------------------
	!PREAMBULO
	
	! PREGUNTA 1
	integer*8::i,inf1,contador1
	integer*16::j,inf2,contador2
	real(kind=4)::A !simple precision
	real(kind=8)::B !doble precision
	real(kind=16)::C !cuadruple precision
	character(len=22)::datos_p1_1='./datos/datos_p1_1.dat'
	character(len=22)::datos_p1_2='./datos/datos_p1_2.dat'

	!-------------------------------------------
 	
 	!simple y doble precisión
	A=1._4
	B=1._8
	inf1=100000000
	contador1=0
	
	!los datos generados los guarda en el
	!directorio ./datos/datos_p1.dat
	open(unit=10,file=datos_p1_1,action='write')
	do i=2,inf1
		!contador
		contador1=contador1+1
		!desarrollo de la serie
		A=A+(1._4/i)
		B=B+(1._8/i)
		!cada 1000 operaciones exportar datos
		if (contador1.eq.1000) then
			if (i.gt.10000) then
				write(10,*)i,B,C
			end if
			contador2=0
		end if
	end do 
	close(unit=10)
	
	!doble y cuadruple precision
	B=1._8
	C=1._16
	inf2=1000000000
	contador2=0
	open(unit=10,file=datos_p1_2,action='write')
	do j=2,inf2
		!contador
		contador2=contador2+1
		!desarrollo de la serie
		B=B+(1._8/j)
		C=C+(1._16/j)
		!cada 1000 operaciones exportar datos
		if (contador2.eq.100000) then
			if (j.gt.10000000) then
				write(10,*)j,B,C
			end if
			contador2=0
		end if
	end do 
	close(unit=10)
	
	
end program

