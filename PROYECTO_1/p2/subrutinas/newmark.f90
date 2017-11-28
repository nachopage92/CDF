subroutine newmark(gama,lambda,t,dt,x,theta,beta,u_0,u_1,v_0,v_1)
!---------------------------
	use modulo_pregunta2
	implicit none
	!entrada
	real(kind=8),intent(in)	:: gama,lambda,t,dt,theta,beta
	real(kind=8),dimension(:),intent(in) :: x
	!entrada-salida
	real(kind=8),dimension(:),intent(inout) :: u_0,u_1,v_0,v_1
	!variables locales
	integer::n,nn,i,j,k
	integer,dimension(:),allocatable::indx
	real(kind=8)::alfa,tmp,d
	real(kind=8),dimension(:),allocatable::f_u,f_v
	real(kind=8),dimension(:,:),allocatable::M,M_inv,AA,tmp_m
!---------------------------
	n=size(x)
	nn=n-2
	alfa = lambda*gama**2._8
	
	allocate(f_u(n),f_v(n))
	f_u(:) = (/ (  beta*(1._8+pi**2._8*gama**2._8)*exp(-t-dt)*sin(pi*x(i)) &
		&+ (0.5_8-beta)*(1._8+pi**2._8*gama**2._8)*exp( -t  )*sin(pi*x(i)) , i=1,n	 ) /)
	f_v(:) = (/ (     theta*(1._8+pi**2._8*gama**2._8)*exp(-t-dt)*sin(pi*x(i)) &
			&+ (1._8-theta)*(1._8+pi**2._8*gama**2._8)*exp( -t  )*sin(pi*x(i)) , i=1,n	 ) /)
	
	!-----------------------------------

	!crear matriz M (matriz del lado izquierdo de la ecuacion)
	allocate(M(nn,nn))
	M(:,:) = 0d0
	M(1,1) = 1._8 + 2._8*alfa*beta ; M(1,2) = - beta*alfa
	do i=2,nn-1
		M(i,i-1:i+1) = (/ - beta*alfa , 1._8 + 2._8*alfa*beta , - beta*alfa /)
	end do
	M(nn,nn-1) = - beta*alfa ; M(nn,nn) = 1._8 + 2._8*alfa*beta
	
	!invertir matriz M
	allocate(indx(nn))
	allocate(M_inv(nn,nn))
	call invertir_matriz(M,nn,1d-10,M_inv)
		
!	tmp_m
	allocate(tmp_m(nn,nn))
	do i=1,nn
		tmp_m(i,:) = (/ (0._8,j=1,nn) /)
	end do
	tmp_m(1,1) = 1._8 - alfa + 2._8*alfa*beta ; tmp_m(1,2) = 0.5_8*alfa - alfa*beta
	do i=2,nn-1
		tmp_m(i,i-1:i+1) = (/ 0.5_8*alfa - alfa*beta , 1._8 - alfa + 2._8*alfa*beta , 0.5_8*alfa - alfa*beta /)
	end do
	tmp_m(nn,nn-1) = 0.5_8*alfa - alfa*beta ; tmp_m(nn,nn) = 1._8 - alfa + 2._8*alfa*beta

!	AA = M_inv * tmp_m
	allocate(AA(nn,nn))
	do i=1,nn
		do j=1,nn
			tmp=0._8
			do k=1,nn
				tmp = tmp + M_inv(i,k)*tmp_m(k,j)
			end do
			AA(i,j) = tmp
		end do
	end do
			
!	u_(n+1) = AA * u_(n) + M_inv * v(n) + dt * M_inv * f_u
	u_1(1) = 0._8
	do i=1,nn
		tmp = 0._8
		do j=1,nn
			tmp = tmp + AA(i,j)*u_0(j+1) + M_inv(i,j)*v_0(j+1)*dt + M_inv(i,j)*f_u(j+1)*dt**2._8
		end do
		u_1(i+1) = tmp
	end do
	u_1(n) = 0._8
	
!	v_(n+1) es explicito
!	v_1(1) = v_0(1) + dt*alfa*theta*(-2._8*u_1(1)+u_1(2)) &
!			&+ dt*alfa*(1._8-theta)*(-2._8*u_0(1)+u_0(2)) &
!			&+ dt*f_v(1)*dt

!	v_1(1) = (3._8*u_1(1)-4._8*u_1(2)+1._8*u_1(3))/(2._8*dt)

	v_1(1)=0d0
	do i=2,n-1
		v_1(i) = v_0(i) + alfa*theta*(1._8/dt)*(u_1(i-1)-2._8*u_1(i)+u_1(i+1)) &
				&+ (1._8/dt)*alfa*(1._8-theta)*(u_0(i-1)-2._8*u_0(i)+u_0(i+1)) &
				&+ dt*f_v(i)
	end do
!	v_1(n) = v_0(n) + dt*alfa*theta*(-2._8*u_1(n)+u_1(n-1)) &
!			&+ dt*alfa*(1._8-theta)*(-2._8*u_0(n)+u_0(n-1)) &
!			&+ f_v(n)*dt
	
!	v_1(n) = (3._8*u_1(n)-4._8*u_1(n-1)+1._8*u_1(n-2))/(2._8*dt)
	
	v_1(n) = 0d0
	
end subroutine

subroutine Ludecomp(a,b,n,tol,x)

	integer::n,er
	integer::o(n)
	double precision::a(n,n),b(n),x(n),tol,s(n)

	er=0
	call Decompose(a,n,tol,o,s,er)
	if(er/=-1)then
		call Substitute(a,o,n,b,x)
	endif

end subroutine

subroutine Decompose(a,n,tol,o,s,er)

	integer::n,er
	integer::o(n)
	double precision::a(n,n),tol,s(n)

	do i=1,n
		o(i)=i
		s(i)=abs(a(i,1))
		do j=2,n
			if(abs(a(i,j))>s(i))then
				s(i)=abs(a(i,j))
			endif
		enddo
	enddo
	do k=1,n-1
		call Pivot(a,o,s,n,k)
		if(abs(a(o(k),k)/s(o(k)))<tol)then
			er=-1
!			write(*,*)a(o(k),k)/s(o(k))
			exit
		endif
		do i=k+1,n
			factor=a(o(i),k)/a(o(k),k)
			a(o(i),k)=factor
			do j=k+1,n
				a(o(i),j)=a(o(i),j)-factor*a(o(k),j)
			enddo
		enddo
	enddo
	if(abs(a(o(k),k)/s(o(k)))<tol)then
		er=-1
!		write(*,*)a(o(k),k)/s(o(k))
	endif

end subroutine

subroutine Pivot(a,o,s,n,k)

	integer::n,k,p
	integer::o(n),dummy1
	double precision::a(n,n),s(n),big,dummy2

	p=k
	big=abs(a(o(k),k)/s(o(k)))
	do ii=k+1,n
		dummy2=abs(a(o(ii),k)/s(o(ii)))
		if(dummy>big)then
			big=dummy2
			p=ii
		endif
	enddo
	dummy1=o(p)
	o(p)=o(k)
	o(k)=dummy1

end subroutine

subroutine Substitute(a,o,n,b,x)

	integer::n
	integer::o(n)
	double precision::a(n,n),b(n),x(n),suma

	do i=2,n
		suma=b(o(i))
		do j=1,i-1
			suma=suma-a(o(i),j)*b(o(j))
		enddo
		b(o(i))=suma
	enddo
	x(n)=b(o(n))/a(o(n),n)
	do i=n-1,1,-1
		suma=0d0
		do j=i+1,n
			suma=suma+a(o(i),j)*x(j)
		enddo
		x(i)=(b(o(i))-suma)/a(o(i),i)
	enddo

end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine invertir_matriz(a,n,tol,ai)

integer::n,er=0
integer::o(n)
double precision::a(n,n),s(n),tol,b(n),x(n),ai(n,n)

	call Decompose(a,n,tol,o,s,er)
	if(er==0)then
		do i=1,n
			do j=1,n
				if(i==j)then
					b(j)=1.
				else
					b(j)=0.
				endif
			enddo
			call Substitute(a,o,n,b,x)
			do j=1,n
				ai(j,i)=x(j)
			enddo
		enddo
	!write(*,*)a
	else
	!write(*,*)'sistema mal condicionado'
	endif

end subroutine
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

