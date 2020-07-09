subroutine rk4(Pr,Ra,beta,t0,dt,y0,y4)
	!-----------------------------
	! esta subrutina implementa el método de 
	! integracion RK4 vectorial (dimension 3)
	! de la forma dy/dt = g(t,y)
	!-----------------------------
	implicit none
	!entrada
	real(kind=8),intent(in)	:: t0,dt,Pr,Ra,beta
	real(kind=8),dimension(3),intent(in)	:: y0
	!salida
	real(kind=8),dimension(3),intent(out)	:: y4
	!local
	real(kind=8)				:: t1,t2
	real(kind=8),dimension(3)	:: y1,y2,y3
	!-----------------------------
	!tiempo
	t1 = t0 + 0.5_8*dt	!medio paso
	t2 = t0 + dt		!un paso
	!pasos rk4	
	y1 = y0 + 0.5_8*dt * g(Pr,Ra,beta,t0,y0)
	y2 = y0 + 0.5_8*dt * g(Pr,Ra,beta,t1,y1)
	y3 = y0 + dt * g(Pr,Ra,beta,t1,y2)	
	y4 = y0 +  (dt/6._8) * &
		& ( g(Pr,Ra,beta,t0,y0) + &
		& 2._8*g(Pr,Ra,beta,t1,y1) + &
		& 2._8*g(Pr,Ra,beta,t1,y2) + &
		& g(Pr,Ra,beta,t2,y3) )
	!-----------------------------
	contains
		function g(Pr,Ra,beta,t,z) result(v)
		! la funcion g representa el lado derecho
		! de la ecuacion (rhs) de lorenz
			real(kind=8),intent(in)::Pr,Ra,beta,t,z(3)
			real(kind=8)::v(3)
			v(1) = Pr*(z(2)-z(1))
			v(2) = Ra*z(1)-z(2)-z(1)*z(3)
			v(3) = z(1)*z(2)-beta*z(3)
		end function g
	!-----------------------------
end subroutine
