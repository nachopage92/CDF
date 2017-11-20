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
	
	integer::nx,nt,k,i,j,q,contador
	integer,dimension(4)::datos_nx,datos_nt
	real(kind=8)::dx,dt,lambda,gama,theta,beta
	real(kind=8),dimension(4)::datos_dx,datos_dt
	real(kind=8),dimension(4,10)::e_lf,e_nm,p_lf,p_nm
	
	real(kind=8),dimension(10)::tiempo
	
	type(y_data),dimension(4)::y_lf,y_nm,dy_nm,y_exacta
	type(xt_data),dimension(4)::x,t
	character(len=11)::data_folder='./datos/S3/'
	
!-------------------------------------------------

	do k=1,4

		!discretizacion
		call discretizacion(k-1,nx,nt,dx,dt)
		
		datos_dx(k)=dx ; datos_dt(k)=dt
		datos_nx(k)=nx ; datos_nt(k)=nt
		
		lambda	= dt/dx
		gama	= sigma/(rho_w*H)
		
		allocate(x(k)%kk_(nx+1),t(k)%kk_(nt+1),y_lf(k)%kk(nt+1,nx+1),&
			&y_nm(k)%kk(nt+1,nx+1),dy_nm(k)%kk(nt+1,nx+1),y_exacta(k)%kk(nt+1,nx+1))
			
		!discretizacion espacio-tiempo
		x(k)%kk_(:) = (/ (dx*(i-1),i=1,nx+1) /)
		t(k)%kk_(:) = (/ (dt*(i-1),i=1,nt+1) /)
		
		!condiciones iniciales esquema leapfrog
		y_lf(k)%kk(1,:) = (/ (sin(pi*x(k)%kk_(i)),i=1,nx+1) /)
		y_lf(k)%kk(2,:) = (/ (exp(-dt)*sin(pi*x(k)%kk_(i)),i=1,nx+1) /)
		
		!condicion inicial esquema newmark
		y_nm(k)%kk(1,:) = (/ (sin(pi*x(k)%kk_(i)),i=1,nx+1) /)
		dy_nm(k)%kk(1,:)= (/ (-sin(pi*x(k)%kk_(i)),i=1,nx+1) /)
	
		!integracion temporal LEAP-FROG
		do i=3,nt+1
			call leap_frog(gama,lambda,t(k)%kk_(i),dt,x(k)%kk_(:),y_lf(k)%kk(i-2,:),y_lf(k)%kk(i-1,:),y_lf(k)%kk(i,:))
		end do

		!integracion temporal NEWMARK
		theta=0.5_8
		beta=0.25_8
		do i=2,nt+1
			call newmark(gama,lambda,t(k)%kk_(i),dt,x(k)%kk_(:),theta,beta,&
				&y_nm(k)%kk(i-1,:),y_nm(k)%kk(i,:),dy_nm(k)%kk(i-1,:),dy_nm(k)%kk(i,:))
		end do

		!solucion analitica
		do i=1,nt+1
			y_exacta(k)%kk(i,:) = (/ (exp(-t(k)%kk_(i))*sin(pi*x(k)%kk_(j)),j=1,nx+1) /)
		end do
		
		if (k.eq.4) then
			!exportar datos
			open(unit=10, status="unknown", action="write",&
				& file=data_folder//'datos_S3_1.txt')
				write(10,*) 'x',t(k)%kk_(:)
				do i=1,nx+1
					write(10,*) x(k)%kk_(i),y_exacta(k)%kk(:,i)
				end do
			close(10)
			open(unit=10, status="unknown", action="write",&
				& file=data_folder//'datos_S3_2.txt')
				write(10,*) 'x',t(k)%kk_(:)
				do i=1,nx+1
					write(10,*) x(k)%kk_(i),y_lf(k)%kk(:,i)
				end do
			close(10)
			open(unit=10, status="unknown", action="write",&
				& file=data_folder//'datos_S3_3.txt')
				write(10,*) 'x',t(k)%kk_(:)
				do i=1,nx+1
					write(10,*) x(k)%kk_(i),y_nm(k)%kk(:,i)
				end do
			close(10)
		end if
	end do

	!generar graficos
	call system('cd gnuplot/ && gnuplot plot_p3s1.txt')

	return
	
	!------------------------------------
	
	!calculo del error
	!p(1,:) => error e_0
	do i=1,10
		p_lf(1,i) = maxval( abs( y_lf(1)%kk(i+1,:) - y_exacta(1)%kk(i+1,:) ) )
	end do
	
	!p(k,j) convergencia en cada paso de tiempo j para cada malla k
	do k=2,4
		q=datos_nt(k)/10 
		contador = 1
		j=1
		do i=1,datos_nt(k)
			if (j.eq.q) then
				p_lf(k,contador) = maxval( abs( y_lf(k)%kk(i,:) - y_exacta(k)%kk(i,:) ) )
				p_lf(k,contador) = log(p_lf(1,contador)/p_lf(k,contador)/log(2._8**dfloat(k-1)))
				tiempo(contador) = t(k)%kk_(i+1)
				contador = contador + 1
				j = 0
			end if
			j = j + 1
		end do
	end do
	
	open(unit=10, status="unknown", action="write", file=data_folder//'TABLA_CONVERGENCIA')
	write(10,*) 't_j(0)', 'p_LF(1)' , 'p_LF(2)' , 'p_LF(3)'
	do i=1,10
		write(10,*) tiempo(i) , p_lf(2,i) , p_lf(3,i) , p_lf(4,i)
	end do
	close(10)
	
!-------------------------------------------------

end program
