subroutine leap_frog(gama,lambda,t,dt,x,y_0,y_1,y_2)
!---------------------------
	use modulo_pregunta2
	implicit none
	!entrada
	real(kind=8),intent(in)	:: gama,lambda,t,dt
	real(kind=8),dimension(:),intent(in) :: x
	!entrada-salida
	real(kind=8),dimension(:),intent(inout) :: y_0,y_1,y_2
	!variables locales
	integer::n,i
	real(kind=8)::alfa
	real(kind=8),dimension(:),allocatable::f
!---------------------------
	n=size(x)
	alfa = lambda*gama**2._8
	allocate(f(n))
	f(:) = (/ ( (1._8+pi**2._8*gama**2._8)*exp(-t)*sin(pi*x(i)) , i=1,n) /)  * dt**2._8
	y_2(1) = 0._8
	y_2(2) = (2._8-2._8*alfa)*y_1(2) + alfa*y_1(3) - y_0(2) + f(2)
	do i=3,n-2
		y_2(i) = alfa*y_1(i-1) + (2._8-2._8*alfa)*y_1(i) + alfa*y_1(i+1) - y_0(i) + f(i)
	end do
	y_2(n-1) = (2._8-2._8*alfa)*y_1(n-1) + alfa*y_1(n-2) - y_0(n-1) + f(n-1)
	y_2(n) = 0._8
end subroutine
