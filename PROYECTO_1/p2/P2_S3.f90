program simulacion_3
!-------------------------------------------------
	use modulo_pregunta2
	implicit none
	interface
		subroutine leap_frog(gama,lambda,t,dt,x,y_0,y_1,y_2)
			implicit none
			real(kind=8),intent(in)	:: gama,lambda,t,dt
			real(kind=8),dimension(:),intent(in) :: x
			real(kind=8),dimension(:),intent(inout) :: y_0,y_1,y_2
		end subroutine
	end interface
	integer::nx,nt,k,i,j
	real(kind=8)::dx,dt,lambda,gama
	real(kind=8),dimension(:),allocatable::x,t
	real(kind=8),dimension(:,:),allocatable::y_lf,y_nw,y_exacta
	character(len=22)::data_file1='./datos/datos_S3_1.txt'
	character(len=22)::data_file2='./datos/datos_S3_2.txt'
!-------------------------------------------------
	!discretizacion
	k=	3 ! k=0,1,2,3
	call discretizacion(k,nx,nt,dx,dt)
	allocate(x(nx+1),t(nt+1),y_lf(nt+1,nx+1),y_nw(nt+1,nx+1),y_exacta(nt+1,nx+1))
	x(:) = (/ (dx*(i-1),i=1,nx+1) /)
	t(:) = (/ (dt*(i-1),i=1,nt+1) /)
	y_lf(1,:) = (/ (sin(pi*x(i)),i=1,nx+1) /)
	y_lf(2,:) = (/ (dx*pi*cos(pi*x(i))+y_lf(1,i),i=1,nx+1) /)
	lambda=dt/dx
	gama=1._8 !sqrt(1._8/(rho_w*H))
	
	!integracion temporal: esquema LEAP-FROG
	do i=3,nt+1
		call leap_frog(gama,lambda,t(i),dt,x,y_lf(i-2,:),y_lf(i-1,:),y_lf(i,:))
	end do
	
	!solucion analitica
	do i=1,nt+1
		y_exacta(i,:) = (/ (exp(-t(i))*sin(pi*x(j)),j=1,nx+1) /)
	end do
	
	!exportar datos
	open(unit=10, status="unknown", action="write", file=data_file1)
		write(10,*) 'x',t(:)
		do i=1,nx+1
			write(10,*) x(i),y_exacta(:,i)
		end do
	close(10)
	open(unit=10, status="unknown", action="write", file=data_file2)
		write(10,*) 'x',t(:)
		do i=1,nx+1
			write(10,*) x(i),y_lf(:,i)
		end do
	close(10)
	
	!generar graficos
	call system('cd gnuplot/ && gnuplot plot_p3s1.txt')
!-------------------------------------------------
end program
