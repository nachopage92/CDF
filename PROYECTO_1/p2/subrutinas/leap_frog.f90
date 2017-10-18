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
	real(kind=8),dimension(:),allocatable::f
!---------------------------
	n=size(x)
	allocate(f(n))
	f(:) = (/ ( (1._8-pi**2._8*gama**2._8)*exp(-t)*sin(pi*x(i)) , i=1,n) /)
	y_2(1) = 2._8-2._8*(lambda*gama)**2._8 * y_1(1) + (lambda*gama)**2._8 * y_1(2) - y_0(1) + f(1) + y_1(n)
	do i=2,n-1
		y_2(i) = (lambda*gama)**2._8*y_1(i-1) + 2._8-2._8*(lambda*gama)**2._8*y_1(i) + (lambda*gama)**2._8*y_1(i+1) - y_0(i) +f(i)
	end do
	y_2(n) = 2._8-2._8*(lambda*gama)**2._8*y_1(n) + (lambda*gama)**2._8*y_1(n-1) - y_0(n) + f(n) + y_1(1)
end subroutine
