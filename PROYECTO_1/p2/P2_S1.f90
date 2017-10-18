program simulacion_1

	use modulo_pregunta2
	implicit none
	
	integer::i,j
	real(kind=np)::alfa,beta,gama,dt,dx
	real(kind=np),dimension(11)::x
	real(kind=np),dimension(26)::t
	real(kind=np),dimension(26,11,2)::y_euler,y_cn
	character(len=22)::data_file1='./datos/datos_S1_1.txt'
	character(len=22)::data_file2='./datos/datos_S1_2.txt'
	
	!--------------------------------------------------
	
	!parametros del problema
	alfa=E/(rho_w*R_0**2._8)
	beta=alfa
	gama=1/(rho_w*H)
	
	!discretizacion espacio-tiempo
	dt=0.0001
	dx=5d-3
	t(:)=(/ ((i-1)*dt,i=1,26) /)
	x(:)=(/ ((i-1)*dx,i=1,11) /)
	
	!condiciones iniciales
	do i=1,11
		y_euler(1,i,:) = (/ 0._8 , 0._8 /)
		y_cn(1,i,:) = (/ 0._8 , 0._8 /)
	end do

	!euler implicito
	do j=2,26
		do i=1,11
			call euler_implicito(alfa,beta,gama,t(j),dt,x(i),y_euler((j-1),i,:),y_euler(j,i,:))
		end do
	end do
	
	!crank nicholson
	do j=2,26
		do i=1,11
			call crank_nicolson(alfa,beta,gama,t(j-1:j),x(i),y_cn((j-1),i,:),y_cn(j,i,:))
		end do
	end do
	
	!exportar datos
	open(unit=10, status="unknown", action="write", file=data_file1)
	write(10,*) 'x',t(:)
	do i=1,11
		write(10,*) x(i),y_euler(:,i,1)
	end do
	close(10)
	open(unit=10, status="unknown", action="write", file=data_file2)
	write(10,*) 'x',t(:)
	do i=1,11
		write(10,*)	x(i),y_cn(:,i,1)
	end do
	close(10)
	
	!graficar
	call system('gnuplot plot_p2s1.txt')
	
end program
