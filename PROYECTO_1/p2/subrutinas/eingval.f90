subroutine eigenval(alfa,beta,lamb1,lamb2)
!---------------------------------------------
	use modulo_pregunta2
	implicit none
	!variables entrada
	real(kind=pcs),intent(in)	:: alfa,beta
	!variables salida
	real(kind=pcs),intent(out)	:: lamb1,lamb2
!---------------------------------------------
	lamb1=5d-1*(-beta+sqrt(beta**2d0-4d0*alfa))
	lamb2=5d-1*(-beta-sqrt(beta**2d0-4d0*alfa))
!---------------------------------------------	
end subroutine
