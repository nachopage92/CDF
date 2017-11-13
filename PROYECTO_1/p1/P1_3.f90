!EL PROGRAMA EJECUTA DOS OPERACIONES DE EJEMPLO
!COMO SE MUESTRA A CONTINUACION

!PRIMER CASO AxB=C (DIMENSIONES CONSISTENTES)

!	A =	1	2	3	4
!		5	6	7	8
!		9	10	11	12
	
!	B =	1	2
!		3	4
!		5	6
!		7	8
	
!SEGUNDO CASO AxD (DIMENSIONES INCOSITENTES)

!	A =	1	2	3	4
!		5	6	7	8
!		9	10	11	12
	
!	D =	1	2
!		3	4
!		5	6

program EJERCICIOS_EN_FORTRAN_3
	!-------------------------------------------
	!INTERFACE
	implicit none
	interface
		subroutine matrix_mult(A,B,C,i)
		real(kind=4),dimension(:,:),intent(in)::A,B
		real(kind=4),dimension(:,:),intent(out)::C
		integer,intent(out)::i
		end subroutine
	end interface
	!-------------------------------------------
	!PREAMBULO
	
	!PREGUNTA 3
	integer::k,i,tmp,tmp_v(2)
	real(kind=4) :: A(3,4),B(4,2),D(3,2)
	real(kind=4),allocatable :: C(:,:)
	
	!-------------------------------------------
	
	!presentacion en pantalla
	write(*,*) '---------------------'
	write(*,*) 'EJERCICIOS_EN_FORTRAN_3'
	write(*,*) 'codigo		: P1_3.f90'
	write(*,*) 'ejecutable	: P1_3'
	write(*,*)
	write(*,*) 'Este programa implementan la'
	write(*,*) "subrutina 'matrix_mult' en"
	write(*,*) 'dos casos posibles:'
	write(*,*)
	write(*,*) '1. multiplicacion de matrices'
	write(*,*) 'con dimensiones consistentes'
	write(*,*) 
	write(*,*) '2. multiplicacion de matrices'
	write(*,*) 'con dimensions inconsistentes'
	write(*,*) '---------------------'
	
	!PRIMER CASO AxB=C (DIMENSIONES CONSISTENTES)
	write(*,*) '---------------------'
	write(*,*) 'PRIMER CASO AxB=C '//&
	&'(DIMENSIONES CONSISTENTES)'
	write(*,*)
	
	!se crean las matrices A, B y D
	A=reshape((/(1._8*k,k=1,12)/),(/3,4/))
	B=reshape((/(1._8*k,k=1,8)/),(/4,2/))
	D=reshape((/(1._8*k,k=1,6)/),(/3,2/))
	
	
	!se muestran las matrices A y B matrices
	write(*,*) 'matriz A'
	do k=1,3
		write(*,*) A(k,:)
	end do
	write(*,*)
	write(*,*) 'matriz B'
	do k=1,4
		write(*,*) B(k,:)
	end do
	write(*,*)
	
	!se ejecuta la subrutina matrix_mult
	call matrix_mult(A,B,C,i)

	!se muestra el resultado C en pantalla
	if (i.eq.0) then
		tmp_v=shape(C)
		tmp=tmp_v(1)
		do k=1,tmp
			write(*,*) C(k,:)
		end do
	else
		write(*,*) 'ERROR!'
		write(*,*) 'DIMENSIONES INCONSISTENTES'
	end if
	
	deallocate(C)
	!-------------------------------------------
	
	!SEGUNDO CASO AxD (DIMENSIONES INCOSITENTES)
	write(*,*) '---------------------'
	write(*,*) 'SEGUNDO CASO AxD '//&
	&'(DIMENSIONES INCONSISTENTES)'
	write(*,*)
	
	!se muestran las matrices A y D matrices
	write(*,*) 'matriz A'
	do k=1,3
		write(*,*) A(k,:)
	end do
	write(*,*)
	write(*,*) 'matriz D'
	do k=1,3
		write(*,*) D(k,:)
	end do
	write(*,*)
	
	!se ejecuta la subrutina matrix_mult
	call matrix_mult(A,D,C,i)
	
	!se muestra el resultado C en pantalla
	if (i.eq.0) then
		tmp_v=shape(C)
		tmp=tmp_v(1)
		do k=1,tmp
			write(*,*) C(k,:)
		end do
	else
		write(*,*) 'ERROR!'
		write(*,*) 'DIMENSIONES INCONSISTENTES'
	end if
	
end program
