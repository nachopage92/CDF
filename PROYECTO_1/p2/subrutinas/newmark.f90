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
	call lu_pivoteado(M,nn,nn,indx,d)
	call lu_inversa(M,nn,nn,indx,M_inv)
		
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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine lu_pivoteado(a, n, np, indx, d)
	
	integer, parameter :: NMAX = 500
	integer i, j, k, imax, n, np, indx(n)
	double precision aamax, dum, suma, vv(NMAX), d, a(np,np)
	double precision, parameter :: min = 1.0d-20

!	Dada una matriz a(1:n,1:n), con dimensiones físicas np x np,
!	esta rutina reemplaza esta por la descomposición LU con 
!	permutaciones por fila en a. La matriz a y el entero n son 
!	entradas. Además a es salida conteniendo en sí a la descomposicion
!   LU omitiendo la la diagonal de 1 de la matriz L; indx(1:n) es un
!	vector de salida que registra las permutaciones por filas efectuadas
!	por el pivoteo parcial; d es una salida +1 o -1 dependiendo si el
!	número de filas intercambiado	es par o impar respectivamente.
!	El arreglo vv registra el escalamiento implicito en cada fila.

	d = 1.0d0

!	Loop sobre las filas para obtener la información del escalamiento 
!	implicito.

	do i = 1, n
	  aamax = 0.0d0
	  do j = 1, n
	    if ( abs(a(i, j)) > aamax ) aamax = abs( a(i, j) )
	  enddo
	  if ( aamax == 0.0d0 ) stop 'Matriz singular, caso 1' 
	  vv(i) = 1.0d0 / aamax
	enddo

!	Ahora se aplicará un loop para el método de Crout (o reducción de Gauss).

	do j = 1, n
	  do i = 1, j - 1
	    suma = a(i, j)
		do k = 1, i - 1
		  suma = suma - a(i, k) * a(k, j)
		enddo
		a(i, j) = suma
	  enddo
	  aamax = 0.0d0

!	Inicio de búsqueda de elemento pivote.

	  do i = j, n
	    suma = a(i, j)
		do k = 1, j - 1
		  suma = suma - a(i, k) * a(k, j)
		enddo
		a(i, j) = suma
		dum = vv(i) * abs(suma)
		if ( dum >= aamax ) then
		  imax = i
		  aamax = dum
		endif
	  enddo
	  if ( j /= imax ) then
	    do k = 1, n
		  dum = a(imax, k)
		  a(imax, k) = a(j, k)
		  a(j, k) = dum
		enddo
		d = - 1.0d0 * d
		vv(imax) =vv(j)
	  endif
	  indx(j) = imax
	  if ( a(j, j) == 0.0d0 ) a(j, j) = min
	  if ( j /= n ) then
	    dum = 1.0d0 / a(j, j)
		do i = j + 1, n
		  a(i, j) = a(i, j) * dum
		enddo
	  endif
	enddo
	
	end subroutine


!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------



	subroutine luevaluacion(a, n, np, indx, b)

	integer i, j, ii, ll, n, np, indx(n)
	double precision a(np, np), b(n), suma

!	Solver para un set de n ecuaciones lineales en la forma Ax = b, donde
!	A es la matriz correspondiente a la descomposición LU determinada en 
!	la subrutina lu_pivoteada. Aquí b(1:n) es una entrada correspondiente
!	al vector columna b del lado derecho del sistema y sale como el vector
!	columna x del sistema (esta entrada es la única que se modifica, por 
!	la que el resto puede ocuparse nuevamente con otro vector b de entrada).
!	El arreglo indx es ingresado como un vector que contiene las
!	permutaciones realizadas en la subrutina lu_pivoteada.

	ii = 0.0d0
	do i = 1, n
	  ll = indx(i)
	  suma = b(ll)
	  b(ll) = b(i)
	  if ( ii /= 0.0d0 ) then
	    do j = ii, i - 1
		  suma = suma - a(i, j) * b(j)
		enddo
	  elseif ( suma /= 0.0d0 ) then
	    ii = i
	  endif
	  b(i) = suma
	enddo

	do i = n, 1, -1
	  suma = b(i)
	  do j = i + 1, n
	    suma = suma - a(i, j) * b(j)
	  enddo
	  b(i) = suma / a(i, i)
	enddo

	end subroutine



!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------



	subroutine lu_inversa(a, n, np, indx, y)

	implicit none
	integer n, np, i, j, indx(np)
	double precision a(np, np), y(np, np), d

	do i = 1, n
	  do j = 1, n
	    y(i, j) = 0.0d0
	  enddo
	  y(i, i) = 1.0d0
	enddo
	
	do j = 1, n
	  call luevaluacion(a, n, np, indx, y(1, j))
	enddo

	end subroutine



!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------
!	----------------------------------------------------------------------



	subroutine numero_condicion(a, y, n, np, num)

!	Esta subrutina esta destinada al cálculo del número de condición de 
!	una matriz a(1:n,1:n) con dimensión física np x np. Aquí a e y son 
!	la matriz A y su inversa A^(-1) respectiva. num es la única salida 
!	que contiene el número de condición calculado mediante la norma 
!	infinito de matrices.

	integer n, np
	double precision a(np,np), y(np,np), num, suma, d, d1

	suma = 0.0d0
	d = 0.0d0
	d1 = 0.0d0

!	Cálculo tamaño de A utilizando norma infinito.

	do i = 1, n
	  suma =0.0d0
	  do j = 1, n
	    suma = suma + abs(a(i,j))
	  enddo
	  if ( abs(suma) >= d ) then
	    d = suma
	  endif
	enddo

!	Cálculo tamaño de A^(-1) utilizando norma infinito.

	do i = 1, n
	  suma =0.0d0
	  do j = 1, n
	    suma = suma + abs(y(i,j))
	  enddo
	  if ( abs(suma) >= d1 ) then
	    d1 = suma
	  endif
	enddo

!	Cálculo número de Condición

	num = d * d1

	end subroutine
