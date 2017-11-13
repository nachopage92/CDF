!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! la subrutina matrix_mult realiza la multiplicacion
! de las matrices A(na,ma) y B(nb,mb). Si ma=nb ento-
! nces se obtiene como resultado una matriz C(na,mb)
! y l=0 . Si ma=/nb entonces l=1
! 
! entrada
! A: matriz lado izquierdo
! B: matriz lado derecho
! 
! salida
! C: matriz resultante
! l: indicador de error
!	l=0 operacion satisfactoria (m. consistentes)
!	l=1 fallo (matrices inconsistentes)
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine matrix_mult(A,B,C,l)
	!entrada
	real(kind=4),dimension(:,:),intent(in)::A,B
	!salida
	integer,intent(out)::l
	real(kind=4),dimension(:,:),allocatable,intent(out)::C
	!variables locales
	integer::i,j
	integer,dimension(2)::dim_a,dim_b,dim_c
	real(kind=4)::tmp
	!-------------------------------------------
	l=0
	dim_a=shape(A)
	dim_b=shape(B)
	if ( dim_a(2) .ne. dim_b(1) ) then
		l=1
		return
	end if
	allocate(C(dim_a(1),dim_b(2)))
	dim_c=(/dim_a(1),dim_b(2)/)
	do i=1,dim_c(1)
		do j=1,dim_c(2)
			tmp=0._8
			do k=1,dim_a(2) !do k=1,dim_b(1)
				tmp=tmp+A(i,k)*b(k,j)
			end do
			C(i,j)=tmp
		end do
	end do
	!-------------------------------------------
end subroutine
