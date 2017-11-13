subroutine euler_implicito(alfa,beta,gama,t,dt,x,y_0,y_1)
!---------------------------------------------
	use modulo_pregunta2
	implicit none
	!variables entrada
	real(kind=np),intent(in) :: alfa,beta,gama,t,dt,x
	!variables entrada-salida
	real(kind=np),dimension(2),intent(inout)	:: y_0,y_1
	!variables locales
	integer							:: i,j
	real(kind=np)					:: suma
	real(kind=np),dimension(2)		:: yn0,yn1,f_0
	real(kind=np),dimension(2,2)	:: Matrix_A
!---------------------------------------------
	Matrix_A(1,:) = (/ 1._8+beta*dt , dt /)
	Matrix_A(2,:) = (/ -alfa*dt , 1._8 /)
	Matrix_A = Matrix_A/(1._8+beta*dt+alfa*dt**2._8)
	f_0(:) = (/ 0._8 , x*dt*gama*dp*(a+b*cos(w*(t-dt))) /)
	do i=1,2
		suma=0d0
		do j=1,2
			suma = suma + Matrix_A(i,j)*(y_0(j)+f_0(j)) 
		enddo
		y_1(i)=suma
	end do
end subroutine
