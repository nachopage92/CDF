module modulo_pregunta2
	implicit none
	integer,parameter::np=kind(1d0)
	real(kind=np),parameter	:: &
	&	L		= 5d-2 ,&
	&	R_0		= 5d-3 ,&
	&	rho_w	= 1d3 ,&
	&	H		= 3d-4 ,&
	&	E		= 9d5 ,&
	&	pi		= 3.1415926535897932d0 ,&
	&	b		= 133.32d0 ,&
	&	dp		= 0.25d0*b ,&
	&	a		= 10*b ,&
	&	w		= 2d0*pi/0.8d0 
end module
