subroutine exportar_datos_S1(x,t,nx,nt,y_euler,y_cn)

	implicit none
	integer,parameter::np=8
	integer,intent(in)::nx,nt
	real(kind=np),dimension(11),intent(in)::x
	real(kind=np),dimension(26),intent(in)::t
	real(kind=np),dimension(26,11,2,2),intent(in)::y_euler,y_cn
	integer::i
	character(len=11)::data_file='./datos/S1/'

!--------------------------- euler implicito

	!-------- beta = sqrt(alfa)

	!grafica para y(x,t) vs x (parametro t) , beta=sqrt(alfa)
!	open(unit=10, status="unknown", action="write", file=data_file//'datos_euler_S1_y1_b1.txt')
!	write(10,*) 'x',t(:)
!	do i=1,nx
!		write(10,*) x(i),y_euler(:,i,1,1)
!	end do
!	close(10)
!	
!	!grafica para dy(x,t) vs x (parametro t) , beta=sqrt(alfa)
!	open(unit=10, status="unknown", action="write", file=data_file//'datos_euler_S1_dy1_b1.txt')
!	write(10,*) 'x',t(:)
!	do i=1,nx
!		write(10,*) x(i),y_euler(:,i,2,1)
!	end do
!	close(10)
	
	!grafica para y(x=0.005,t) vs t , beta=sqrt(alfa)
	open(unit=10, status="unknown", action="write", file=data_file//'datos_euler_S1_y2_b1.txt')
	!write(10,*) 'x',t(:)
	do i=1,nt
		write(10,*) t(i),y_euler(i,nx,1,1)
	end do
	close(10)
	
	!grafica para dy(x=0.005,t) vs t , beta=sqrt(alfa)
	open(unit=10, status="unknown", action="write", file=data_file//'datos_euler_S1_dy2_b1.txt')
	!write(10,*) 'x',t(:)
	do i=1,nt
		write(10,*) t(i),y_euler(i,nx,2,1)
	end do
	close(10)
	
	!---------- beta = alfa
	
	!grafica para y(x,t) vs x (parametro t) , beta = alfa
!	open(unit=10, status="unknown", action="write", file=data_file//'datos_euler_S1_y1_b2.txt')
!	write(10,*) 'x',t(:)
!	do i=1,nx
!		write(10,*) x(i),y_euler(:,i,1,2)
!	end do
!	close(10)
!	
!	!grafica para dy(x,t) vs x (parametro t) , beta = alfa
!	open(unit=10, status="unknown", action="write", file=data_file//'datos_euler_S1_dy1_b2.txt')
!	write(10,*) 'x',t(:)
!	do i=1,nx
!		write(10,*) x(i),y_euler(:,i,2,2)
!	end do
!	close(10)
	
	!grafica para y(x=0.005,t) vs t , beta = alfa
	open(unit=10, status="unknown", action="write", file=data_file//'datos_euler_S1_y2_b2.txt')
	!write(10,*) 'x',t(:)
	do i=1,nt
		write(10,*) t(i),y_euler(i,nx,1,2)
	end do
	close(10)
	
	!grafica para dy(x=0.005,t) vs t , beta = alfa
	open(unit=10, status="unknown", action="write", file=data_file//'datos_euler_S1_dy2_b2.txt')
	!write(10,*) 'x',t(:)
	do i=1,nt
		write(10,*) t(i),y_euler(i,nx,2,2)
	end do
	close(10)

	!--------------------------- crank nicolson

	!-------- beta = sqrt(alfa)

	!grafica para y(x,t) vs x (parametro t) , beta=sqrt(alfa)
!	open(unit=10, status="unknown", action="write", file=data_file//'datos_cn_S1_y1_b1.txt')
!	write(10,*) 'x',t(:)
!	do i=1,nx
!		write(10,*) x(i),y_cn(:,i,1,1)
!	end do
!	close(10)
!	
!	!grafica para dy(x,t) vs x (parametro t) , beta=sqrt(alfa)
!	open(unit=10, status="unknown", action="write", file=data_file//'datos_cn_S1_dy1_b1.txt')
!	write(10,*) 'x',t(:)
!	do i=1,nx
!		write(10,*) x(i),y_cn(:,i,2,1)
!	end do
!	close(10)
	
	!grafica para y(x=0.005,t) vs t , beta=sqrt(alfa)
	open(unit=10, status="unknown", action="write", file=data_file//'datos_cn_S1_y2_b1.txt')
	!write(10,*) 'x',t(:)
	do i=1,nt
		write(10,*) t(i),y_cn(i,nx,1,1)
	end do
	close(10)
	
	!grafica para dy(x=0.005,t) vs t , beta=sqrt(alfa)
	open(unit=10, status="unknown", action="write", file=data_file//'datos_cn_S1_dy2_b1.txt')
	!write(10,*) 'x',t(:)
	do i=1,nt
		write(10,*) t(i),y_cn(i,nx,2,1)
	end do
	close(10)
	
	!---------- beta = alfa
	
	!grafica para y(x,t) vs x (parametro t) , beta = alfa
!	open(unit=10, status="unknown", action="write", file=data_file//'datos_cn_S1_y1_b2.txt')
!	write(10,*) 'x',t(:)
!	do i=1,nx
!		write(10,*) x(i),y_cn(:,i,1,2)
!	end do
!	close(10)
!	
!	!grafica para dy(x,t) vs x (parametro t) , beta = alfa
!	open(unit=10, status="unknown", action="write", file=data_file//'datos_cn_S1_dy1_b2.txt')
!	write(10,*) 'x',t(:)
!	do i=1,nx
!		write(10,*) x(i),y_cn(:,i,2,2)
!	end do
!	close(10)
	
	!grafica para y(x=0.005,t) vs t , beta = alfa
	open(unit=10, status="unknown", action="write", file=data_file//'datos_cn_S1_y2_b2.txt')
	!write(10,*) 'x',t(:)
	do i=1,26
		write(10,*) t(i),y_cn(i,nx,1,2)
	end do
	close(10)
	
	!grafica para dy(x=0.005,t) vs t , beta = alfa
	open(unit=10, status="unknown", action="write", file=data_file//'datos_cn_S1_dy2_b2.txt')
	!write(10,*) 'x',t(:)
	do i=1,nt
		write(10,*) t(i),y_cn(i,nx,2,2)
	end do
	close(10)
	
end subroutine
