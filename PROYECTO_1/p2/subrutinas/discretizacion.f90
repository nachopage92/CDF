subroutine discretizacion(k,nx,nt,dx,dt)
	implicit none
	!entrada
	integer,intent(in)::k
	!salida
	integer,intent(out)::nx,nt
	real(kind=8),intent(out)::dx,dt
!-----------------------------------
	dt = 1._8 / (2._8**k*10)
	dx = dt
	nt = nint(1._8/dt)
	nx = nint(1._8/dx)
end subroutine
