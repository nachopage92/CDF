subroutine eigenval(alfa,beta,lamb1,lamb2)
!---------------------------------------------
	use modulo_pregunta2
	implicit none
	!variables entrada
	real(kind=np),intent(in)	:: alfa,beta
	!variables salida
	complex(kind=np),intent(out)	:: lamb1,lamb2
!---------------------------------------------
	if ( beta .ge. 2._8*sqrt(alfa) ) then 
		lamb1 = complex( 0.5_8*(-beta+sqrt(beta**2d0-4d0*alfa)) , 0._8 )
		lamb2 = complex( 0.5_8*(-beta-sqrt(beta**2d0-4d0*alfa)) , 0._8 )
	else
		lamb1 = complex( -0.5_8*beta , 0.5_8*sqrt(-beta**2d0+4d0*alfa) ) 
		lamb2 = complex( -0.5_8*beta , -0.5_8*sqrt(-beta**2d0+4d0*alfa) ) 
	end if
!---------------------------------------------	
end subroutine
