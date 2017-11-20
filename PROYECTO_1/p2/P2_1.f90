program simulacion_1

	use modulo_pregunta2
	implicit none
	
	integer::i,j,k,nt
	integer,parameter::nx=11
	real(kind=np)::alfa,gama,dt,dx
	real(kind=np),dimension(2)::beta
	complex(kind=np),dimension(2,2)::lambda
	real(kind=np),dimension(nx)::x
	real(kind=np),dimension(:),allocatable::t
	real(kind=np),dimension(:,:,:,:),allocatable::y_euler,y_cn
	character(len=11)::data_file1='./datos/S1/'
	character(len=11)::data_file2='./datos/S2/'
	character(len=1),dimension(2)::data_beta
	
	!--------------------------------------------------
	!				SIMULACION 1
	!--------------------------------------------------
	
	nt = 26
	allocate(t(nt),y_euler(nt,nx,2,2),y_cn(nt,nx,2,2))
	
	!parametros del problema
	alfa=E/(rho_w*R_0**2._8)
	gama=1._8/(rho_w*H)
	beta(1) = sqrt(alfa)
	beta(2) = alfa	
	
	!discretizacion espacio-tiempo (t=0)
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

!	!calculo valores propios de A
	call eigenval(alfa,beta(1),lambda(1,1),lambda(1,2))
	call eigenval(alfa,beta(2),lambda(2,1),lambda(2,2))

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
	
	!--------------- exportar datos--------------------

	data_beta(1) = '1' ; data_beta(2) = '2' 
	do k=1,2
		!grafica para y_euler(x=0.005,t) vs t 
		open(unit=10, action="write", &
		&file=data_file1//'datos_euler_S1_y_b'//data_beta(k)//'.txt')
			do i=1,nt
				write(10,*) t(i) , y_euler(i,nx,1,k)
			end do
		close(10)
		!grafica para dy_euler(x=0.005,t) vs t
		open(unit=10, action="write", &
		&file=data_file1//'datos_euler_S1_dy_b'//data_beta(k)//'.txt')
			do i=1,nt
				write(10,*) t(i) , y_euler(i,nx,2,k)
			end do
		close(10)
		!grafica para y_cn(x=0.005,t) vs t 
		open(unit=10, action="write",&
		&file=data_file1//'datos_cn_S1_y_b'//data_beta(k)//'.txt')
			do i=1,nt
				write(10,*) t(i) , y_cn(i,nx,1,k)
			end do
		close(10)
		!grafica para dy_cn(x=0.005,t) vs t
		open(unit=10, action="write",&
		&file=data_file1//'datos_cn_S1_dy_b'//data_beta(k)//'.txt')
			do i=1,nt
				write(10,*) t(i) , y_cn(i,nx,2,k)
			end do
		close(10)
	end do
		
	deallocate(t,y_euler,y_cn)
		
	!---------------- graficar ------------------------
	
	call system('cd gnuplot/ && gnuplot plot_p2_1.txt')
	

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	
	!--------------------------------------------------
	!				SIMULACION 2
	!--------------------------------------------------
	
	!se conservan los parámetros, cambia el paso y el 
	!intervalo de tiempo
	nt = 100
	allocate(t(nt),y_euler(nt,nx,2,2),y_cn(nt,nx,2,2))
	dt = 0.1_8
	t(:)=(/ ((i-1)*dt,i=1,nt) /)
	
	!condiciones iniciales (t=0)
	do i=1,nx
		y_euler(1,i,:,1) = (/ 0._8 , 0._8 /)
		y_euler(1,i,:,2) = (/ 0._8 , 0._8 /)
		y_cn(1,i,:,1) = (/ 0._8 , 0._8 /)
		y_cn(1,i,:,2) = (/ 0._8 , 0._8 /)
	end do
	
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
	
	!--------------- exportar datos--------------------

	data_beta(1) = '1' ; data_beta(2) = '2' 
	do k=1,2
		!grafica para y_euler(x=0.005,t) vs t 
		open(unit=10, action="write", &
		&file=data_file2//'datos_euler_S2_y_b'//data_beta(k)//'.txt')
			do i=1,nt
				write(10,*) t(i) , y_euler(i,nx,1,k)
			end do
		close(10)
		!grafica para dy_euler(x=0.005,t) vs t
		open(unit=10, action="write", &
		&file=data_file2//'datos_euler_S2_dy_b'//data_beta(k)//'.txt')
			do i=1,nt
				write(10,*) t(i) , y_euler(i,nx,2,k)
			end do
		close(10)
		!grafica para y_cn(x=0.005,t) vs t 
		open(unit=10, action="write",&
		&file=data_file2//'datos_cn_S2_y_b'//data_beta(k)//'.txt')
			do i=1,nt
				write(10,*) t(i) , y_cn(i,nx,1,k)
			end do
		close(10)
		!grafica para dy_cn(x=0.005,t) vs t
		open(unit=10, action="write",&
		&file=data_file2//'datos_cn_S2_dy_b'//data_beta(k)//'.txt')
			do i=1,nt
				write(10,*) t(i) , y_cn(i,nx,2,k)
			end do
		close(10)
	end do
		
	deallocate(t,y_euler,y_cn)
	
	!---------------- graficar ------------------------
	
	call system('cd gnuplot/ && gnuplot plot_p2_1.txt')
		
end program
