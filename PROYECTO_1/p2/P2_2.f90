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
	character(len=11)::data_folder='./datos/S3/'
	character(len=20)::data_folder_informe='./../INFORME/parte3/'
	
	!para la simulación 2
	integer::nt
	real(kind=8)::y_lf1(101,11),y_lf2(401,11),y_nm1(101,11),dy_nm1(101,11),&
		&y_nm2(401,11),dy_nm2(401,11),dt,TT,xx(11),vt1(101),vt2(401),lambda
	
!-------------------------------------------------

	!SIMULACIÓN 1

	do k=1,4

		!discretizacion
		call discretizacion(k-1,nx,dx,nt_lf,dt_lf,nt_nm,dt_nm)
		
		!espacial
		datos_dx(k)=dx
		datos_nx(k)=nx
		
		!leapfrog
		datos_dt(1,k)=dt_lf
		datos_nt(1,k)=nt_lf
		
		!newman
		datos_dt(2,k)=dt_nm
		datos_nt(2,k)=nt_nm
		
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

		!solucion analitica
		do i=1,nt_lf+1
			y_exacta_lf(k)%kk(i,:) = (/ (exp(-t_lf(k)%kk_(i))*sin(pi*x(k)%kk_(j)),j=1,nx+1) /)
		end do
		do i=1,nt_nm+1
			y_exacta_nm(k)%kk(i,:) = (/ (exp(-t_nm(k)%kk_(i))*sin(pi*x(k)%kk_(j)),j=1,nx+1) /)
		end do
		
		if (k.eq.4) then
			!exportar datos
			open(unit=10, status="unknown", action="write",&
				& file=data_folder//'datos_S3_1.txt')
				write(10,*) 'x',t_nm(k)%kk_(:)
				do i=1,nx+1
					write(10,*) x(k)%kk_(i),y_exacta_nm(k)%kk(:,i)
				end do
			close(10)
			open(unit=10, status="unknown", action="write",&
				& file=data_folder//'datos_S3_2.txt')
				write(10,*) 'x',t_lf(k)%kk_(:)
				do i=1,nx+1
					write(10,*) x(k)%kk_(i),y_lf(k)%kk(:,i)
				end do
			close(10)
			open(unit=10, status="unknown", action="write",&
				& file=data_folder//'datos_S3_3.txt')
				write(10,*) 'x',t_nm(k)%kk_(:)
				do i=1,nx+1
					write(10,*) x(k)%kk_(i),y_nm(k)%kk(:,i)
				end do
			close(10)
		end if
	end do
	
	!------------------------------------
	
	!calculo del error
	!p(1,:) => error e_0
	do i=1,10
		!p_lf(1,i) = maxval( abs( y_lf(1)%kk(i+1,:) - y_exacta(1)%kk(i+1,:) ) )
		!p_nm(1,i) = maxval( abs( y_nm(1)%kk(i+1,:) - y_exacta(1)%kk(i+1,:) ) )
		p_lf(1,i) = sqrt(sum( ( y_lf(1)%kk(i+1,:) - y_exacta_lf(1)%kk(i+1,:) )**2._8 ))/datos_nx(1)
		p_nm(1,i) = sqrt(sum( ( y_nm(1)%kk(i+1,:) - y_exacta_nm(1)%kk(i+1,:) )**2._8 ))/datos_nx(1)
	end do
	
	!p(k,j) convergencia en cada paso de tiempo j para cada malla k
	do k=2,4
	
		q1=datos_nt(1,k)/10 
		q2=datos_nt(2,k)/10
		contador1 = 1
		contador2 = 1
		j=1
		jj=1

		!LEAP-FROG		
		do i=1,datos_nt(1,k)
			if (j.eq.q1) then
				!p_lf(k,contador1) = maxval( abs( y_lf(k)%kk(i,:) - y_exacta_lf(k)%kk(i,:) ) )
				p_lf(k,contador1) = sqrt(sum( ( y_lf(k)%kk(i,:) - y_exacta_lf(k)%kk(i,:) )**2._8 ))/datos_nx(k)
				p_lf(k,contador1) = log(p_lf(1,contador1)/p_lf(k,contador1)/log(2._8**dfloat(k-1)))
				tiempo(1,contador1) = t_lf(k)%kk_(i+1)
				contador1 = contador1 + 1
				j = 0
			end if
			j = j + 1
		end do
		
		!NEWMARK
		do i=1,datos_nt(2,k)
			if (jj.eq.q2) then
				!p_nm(k,contador2) = maxval( abs( y_nm(k)%kk(i,:) - y_exacta_nm(k)%kk(i,:) ) )
				p_nm(k,contador2) = sqrt(sum( ( y_nm(k)%kk(i,:) - y_exacta_nm(k)%kk(i,:) )**2._8 ))/datos_nx(k)
				p_nm(k,contador2) = log(p_nm(1,contador2)/p_nm(k,contador2)/log(2._8**dfloat(k-1)))
				tiempo(2,contador2) = t_nm(k)%kk_(i+1)
				contador2 = contador2 + 1
				jj = 0
			end if
			jj= jj + 1
		end do
		
	end do
	
	!genera tabla LF
	open(unit=10, status="unknown", action="write", file=data_folder//'TABLA_CONVERGENCIA_LEAPFROG.txt')
	!'t_j(0)', 'p_LF(1)' , 'p_LF(2)' , 'p_LF(3)'
	do i=1,10
		write(10,'(4F8.4)') tiempo(1,i) , p_lf(2,i) , p_lf(3,i) , p_lf(4,i)
	end do
	close(10)
	
	!genera tabla NW
	open(unit=10, status="unknown", action="write", file=data_folder//'TABLA_CONVERGENCIA_NEWMARK.txt')
	! 't_j(0)', 'p_NM(1)' , 'p_NM(2)' , 'p_NM(3)'
	do i=1,10
		write(10,'(4F8.4)') tiempo(2,i) , p_nm(2,i) , p_nm(3,i) , p_nm(4,i)
	end do
	close(10)
	
	!exportar tabla al informe
100 format(' ',2(F8.4,A3,F8.4,A3,F8.4,A3,F8.4,A4))
	open(unit=10, status="unknown", action="write", file=data_folder_informe//'TABLA_CONVERGENCIA.tex')
	!'t_j(0)', 'p_LF(1)' , 'p_LF(2)' , 'p_LF(3)'
	do i=1,10
		write(10,100) tiempo(1,i),'&',p_lf(2,i),'&',p_lf(3,i),'&',p_lf(4,i),&
			&'&', tiempo(2,i),'&',p_nm(2,i),'&',p_nm(3,i),'&',p_nm(4,i),'\\'
	end do
	close(10)
	
	!generar animaciones
	call system('cd gnuplot/ && gnuplot plot_p3s1.txt')
	
!-------------------------------------------------

end program
