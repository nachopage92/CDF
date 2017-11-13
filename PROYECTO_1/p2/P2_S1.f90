program simulacion_1

	use modulo_pregunta2
	implicit none
	
	integer::i,j,k
	integer,parameter::nx=11,nt=26
	real(kind=np)::alfa,gama,dt,dx
	real(kind=np),dimension(2)::beta
	complex(kind=np),dimension(2,2)::lambda
	real(kind=np),dimension(nx)::x
	real(kind=np),dimension(nt)::t
	real(kind=np),dimension(nt,nx,2,2)::y_euler,y_cn
	
	!--------------------------------------------------
	
	!parametros del problema
	alfa=E/(rho_w*R_0**2._8)
	gama=1._8/(rho_w*H)
	beta(1) = sqrt(alfa)
	beta(2) = alfa	
	
	!discretizacion espacio-tiempo
	dt=0.0001_8
	dx=0.005_8
	t(:)=(/ ((i-1)*dt,i=1,nt) /)
	x(:)=(/ ((i-1)*dx,i=1,nx) /)
	
	!condiciones iniciales
	do i=1,nx
		y_euler(1,i,:,1) = (/ 0._8 , 0._8 /)
		y_euler(1,i,:,2) = (/ 0._8 , 0._8 /)
		y_cn(1,i,:,1) = (/ 0._8 , 0._8 /)
		y_cn(1,i,:,2) = (/ 0._8 , 0._8 /)
	end do

	!calculo valores propios de A
	call eigenval(alfa,beta(1),lambda(1,1),lambda(1,2))
	call eigenval(alfa,beta(2),lambda(2,1),lambda(2,2))

	write(*,*) alfa
	write(*,*) sqrt(alfa)
	write(*,*)
	write(*,*)lambda(1,1)
	write(*,*)lambda(1,2)
	write(*,*)
	write(*,*)lambda(2,1)
	write(*,*)lambda(2,2)
	return	

	!euler implicito
	do k=1,2			! beta=sqrt(alfa) y beta=alfa
		do j=2,nt		! discretizacion tiempo
			do i=1,nx	! discretizacion espacio
				call euler_implicito(alfa,beta(k),gama,t(j),dt,x(i),y_euler((j-1),i,:,k),y_euler(j,i,:,k))
			end do
		end do
	end do
	
	
	!crank nicholson
	do k=1,2			! beta=sqrt(alfa) y beta=alfa
		do j=2,nt		! discretizacion tiempo
			do i=1,nx	! discretizacion espacio
				call crank_nicolson(alfa,beta(k),gama,t(j-1:j),x(i),y_cn((j-1),i,:,k),y_cn(j,i,:,k))
			end do
		end do
	end do
	
	!--------------------------- exportar datos--------------------------------------

	call exportar_datos_S1(x(:),t(:),nx,nt,y_euler(:,:,:,:),y_cn(:,:,:,:))
	
	!--------------------------- graficar ----------------------------------------
	
	call system('cd gnuplot/ && gnuplot plot_p2s1.txt')
	
end program
