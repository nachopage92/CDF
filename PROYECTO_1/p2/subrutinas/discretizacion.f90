subroutine discretizacion(k,nx,dx,nt_lf,dt_lf,nt_nm,dt_nm)
	implicit none
	!entrada
	integer,intent(in)::k
	!salida
	integer,intent(out)::nx,nt_lf,nt_nm
	real(kind=8),intent(out)::dx,dt_lf,dt_nm
!-----------------------------------
	! espacial
	dx = 1._8 / (2._8**k*10)
	nx = nint(1._8/dx)
	! tiempo lf
	dt_lf = 0.25_8*dx
	nt_lf = nint(1._8/dt_lf)
	! tiempo nm
	dt_nm = dx
	nt_nm = nint(1._8/dt_nm)
end subroutine
