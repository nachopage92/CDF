program simulacion_3_2
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
		subroutine newmark(gama,lambda,t,dt,x,theta,beta,u_0,u_1,v_0,v_1)
			implicit none
			real(kind=8),intent(in)	:: gama,lambda,t,dt,theta,beta
			real(kind=8),dimension(:),intent(in) :: x
			real(kind=8),dimension(:),intent(inout) :: u_0,u_1,v_0,v_1
		end subroutine
	end interface
	
	type y_data
		real(kind=8),dimension(:,:),allocatable::kk
	end type
	type xt_data
		real(kind=8),dimension(:),allocatable::kk_
	end type
	
	integer::nx,nt_lf,nt_nm,k,i,j,q1,q2,contador1,contador2,ii,jj
	integer,dimension(4)::datos_nx
	integer,dimension(2,4)::datos_nt
	real(kind=8),dimension(4)::datos_dx
	real(kind=8),dimension(2,4)::datos_dt
	real(kind=8)::dx,dt_lf,dt_nm,lambda_lf,lambda_nm,gama,theta,beta
	real(kind=8),dimension(4,10)::e_lf,e_nm,p_lf,p_nm
	real(kind=8),dimension(2,10)::tiempo
	type(y_data),dimension(4)::y_lf,y_nm,dy_nm,y_exacta_lf,y_exacta_nm
	type(xt_data),dimension(4)::x,t_lf,t_nm
	character(len=11)::data_folder='./datos/S4/'
	character(len=20)::data_folder_informe='./../INFORME/parte3/'
	
!-------------------------------------------------

	!SIMULACIÓN 1

	do k=1,4

		!discretizacion
			dx=0.05_8
			nx=11
			dt_lf=0.01_8
			nt_lf=101
			dt_nm=0.01_8
			nt_nm=101
		
!		!espacial
!		datos_dx(k)=dx
!		datos_nx(k)=nx
!		
!		!leapfrog
!		datos_dt(1,k)=dt_lf
!		datos_nt(1,k)=nt_lf
!		
!		!newman
!		datos_dt(2,k)=dt_nm
!		datos_nt(2,k)=nt_nm
		
		!parametros
		lambda_lf	= dt_lf/dx
		lambda_nm	= dt_nm/dx
		gama	= sigma/(rho_w*H)
		
		allocate(x(k)%kk_(nx+1),t_lf(k)%kk_(nt_lf+1),t_nm(k)%kk_(nt_nm+1),&
			&y_lf(k)%kk(nt_lf+1,nx+1),y_nm(k)%kk(nt_nm+1,nx+1),dy_nm(k)%kk(nt_nm+1,nx+1),&
			&y_exacta_lf(k)%kk(nt_lf+1,nx+1),y_exacta_nm(k)%kk(nt_nm+1,nx+1))
			
		!discretizacion espacio-tiempo
		x(k)%kk_(:) = (/ (dx*(i-1),i=1,nx+1) /)
		t_lf(k)%kk_(:) = (/ (dt_lf*(i-1),i=1,nt_lf+1) /)
		t_nm(k)%kk_(:) = (/ (dt_nm*(i-1),i=1,nt_nm+1) /)
		
		!condicion inicial esquema newmark
		y_nm(k)%kk(1,:) = (/ (sin(pi*x(k)%kk_(i)),i=1,nx+1) /)
		dy_nm(k)%kk(1,:)= (/ (-sin(pi*x(k)%kk_(i)),i=1,nx+1) /)
		
		!condiciones iniciales esquema leapfrog
		y_lf(k)%kk(1,:) = (/ (sin(pi*x(k)%kk_(i)),i=1,nx+1) /)
		y_lf(k)%kk(2,:) = y_lf(k)%kk(1,:) + dt_lf*(/ (-sin(pi*x(k)%kk_(i)),i=1,nx+1) /)
	
		!integracion temporal LEAP-FROG
		do i=3,nt_lf+1
			call leap_frog(gama,lambda_lf,t_lf(k)%kk_(i),dt_lf,x(k)%kk_(:),&
				&y_lf(k)%kk(i-2,:),y_lf(k)%kk(i-1,:),y_lf(k)%kk(i,:))
		end do

		!integracion temporal NEWMARK
		theta=0.5_8
		beta=0.25_8
		do i=2,nt_nm+1
			call newmark(gama,lambda_nm,t_nm(k)%kk_(i),dt_nm,x(k)%kk_(:),theta,beta,&
				&y_nm(k)%kk(i-1,:),y_nm(k)%kk(i,:),dy_nm(k)%kk(i-1,:),dy_nm(k)%kk(i,:))
		end do
		
		if (k.eq.1) then
			!exportar datos
			
			open(unit=10, status="unknown", action="write",&
				& file=data_folder//'datos_P3_S2_1.dat')
				do i=1,nx+1
					do j=1,nt_nm+1
						write(10,*) x(k)%kk_(i),t_lf(k)%kk_(j),y_lf(k)%kk(j,i)
					end do
				end do
			close(10)
			open(unit=10, status="unknown", action="write",&
				& file=data_folder//'datos_P3_S2_2.dat')
				do i=1,nx+1
					do j=1,nt_nm+1
						write(10,*) x(k)%kk_(i),t_nm(k)%kk_(j),y_nm(k)%kk(j,i)
					end do
				end do
			close(10)
		end if
		
		!call system('cd gnuplot/ && gnuplot plot_p3s1.txt')
		
		return
		
	end do

!-------------------------------------------------

end program
