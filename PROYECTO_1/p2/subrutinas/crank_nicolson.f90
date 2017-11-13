subroutine crank_nicolson(alfa,beta,gama,t,x,y_0,y_1)
!---------------------------------------------
	use modulo_pregunta2
	implicit none
	!variables entrada
	real(kind=np),intent(in) :: alfa,beta,gama,x
	real(kind=np),dimension(2),intent(in)::t
	!variables entrada-salida
	real(kind=np),dimension(2),intent(inout)	:: y_0,y_1
	!variables locales
	integer							:: i,j
	real(kind=np)					:: suma,dt
	real(kind=np),dimension(2)		:: yn0,yn1,f_0
	real(kind=np),dimension(2,2)	:: Matrix_A,Matrix_B
!---------------------------------------------
	dt = t(2) - t(1)
	
	!matriz A
	Matrix_A(1,:) = (/ 1._8 + beta*dt*0.5_8 - alfa*dt**2._8*0.25_8 ,  dt /)
	Matrix_A(2,:) = (/ -alfa*dt , 1._8 - beta*dt*0.5_8 - alfa*0.25_8*dt**2._8 /)
	Matrix_A = Matrix_A/(1._8 + beta*dt*0.5_8 + alfa*dt**2._8*0.25_8)
	
	!matriz B
	Matrix_B(1,:) = (/ 1._8 + beta*dt*0.5_8 , dt*0.5_8 /)
	Matrix_B(2,:) = (/ -alfa*dt*0.5_8 , 1._8 /)
	Matrix_B = Matrix_B/(1._8 + beta*dt*0.5_8 + alfa*dt**2._8*0.25_8)
	
	!vector f_0
	f_0(:) = (/ 0._8 , x*dt*gama*dp*(a+b*cos(w*t(1)))*0.5_8 + x*dt*gama*dp*(a+b*cos(w*t(2)))*0.5_8  /)
	
	!calcular y_(n+1) y dy_(n+1) 
	do i=1,2
		suma=0d0
		do j=1,2
			suma = suma + Matrix_A(i,j)*y_0(j) + Matrix_B(i,j)*f_0(j) 
		enddo
		y_1(i)=suma
	end do
	
end subroutine
