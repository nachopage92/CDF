program Pregunta_3
	!--------------------------------
	implicit none
	integer::n,m,i,j
	real(kind=8) :: tiempo_total,dt,Pr,beta,t
	real(kind=8),dimension(3) :: z0,z1,Ra
	real(kind=8),dimension(11)::Ra_2
	character(len=8)::folder='./datos/'
	character(len=2),dimension(3) :: ra_n
	character(len=2),dimension(11) :: ra_n2
	!--------------------------------

	!discretizacion tiempo	
	tiempo_total = 6d1
	t	= 0._8
	dt	= 0.01_8
	m	= nint(tiempo_total/dt)

	!parametros
	Pr	= 10._8
	Ra	= (/ 0.5_8 , 10._8 , 28._8 /)
	beta = 8._8/3._8

	!resolucion atractor de lorenz
	ra_n(1) = '05' ; ra_n(2) = '10' ; ra_n(3) = '28' 
	do j=1,3
		open(unit=10, file=folder//'datos_P3_'//ra_n(j)//'.dat', action='write')
		!inicializacion
		n = 0
		z0(:) = (/ 0._8 , 1._8 , 0._8 /)
		write(10,*) n , t , 10._8*z0(1) , 10._8*z0(2) , 10._8*z0(3)
		!subrutina RK4 -> Ra = 0.5
		do i=1,m
			call rk4(Pr,Ra(j),beta,t,dt,z0,z1)
			z0 = z1
			n = n+1
			t = t+dt
			if ( mod(n,5) .eq. 0 ) then
				write(10,*) n , t , 10._8*z0(1) , 10._8*z0(2) , 10._8*z0(3)
			end if
		end do
		close(unit=10)
	end do

!----------------------------------------

	!variar lentamente Ra
	Ra_2(:) = (/ (3.0*(j-1) , j=1,11) /)

	ra_n2(1) = '00'
	ra_n2(2) = '03'
	ra_n2(3) = '06'
	ra_n2(4) = '09'
	ra_n2(5) = '12'
	ra_n2(6) = '15'
	ra_n2(7) = '18'
	ra_n2(8) = '21'
	ra_n2(9) = '24'
	ra_n2(10) = '27'
	ra_n2(11) = '30'

	do j=1,11
		open(unit=10, file=folder//'datos_P3_'//ra_n2(j)//'.dat', action='write')
		!inicializacion
		n = 0
		z0(:) = (/ 0._8 , 1._8 , 0._8 /)
		write(10,*) n , t , 10._8*z0(1) , 10._8*z0(2) , 10._8*z0(3)
		!subrutina RK4 -> Ra = 0.5
		do i=1,m
			call rk4(Pr,Ra_2(j),beta,t,dt,z0,z1)
			z0 = z1
			n = n+1
			t = t+dt
			if ( mod(n,5) .eq. 0 ) then
				write(10,*) n , t , 10._8*z0(1) , 10._8*z0(2) , 10._8*z0(3)
			end if
		end do
		close(unit=10)
	end do
	
	call system('cd gnuplot/ && gnuplot plot_P3.txt')
end program

