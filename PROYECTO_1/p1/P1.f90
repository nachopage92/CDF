program Ejercicio_1
	!-------------------------------------------
	!INTERFACE
	implicit none
	interface
		subroutine matrix_mult(A,B,C)
		real(kind=8),dimension(:,:),intent(in)::A,B
		real(kind=8),dimension(:,:),intent(out)::C
		end subroutine
	end interface
	!-------------------------------------------
	!PREAMBULO
	! PREGUNTA 1
	integer::i,inf,contador
	real(kind=4)::A
	real(kind=8)::B
	! PREGUNTA 2
	integer::m,u0,u1,u2,j
	!PREGUNTA 3
	integer::k
	real(kind=8)::AA(3,4),BB(4,2),CC(3,2)

	!-------------------------------------------
!	!DESARROLLO
! 
	!PREGUNTA 1
	A=1._4
	B=1._8
	inf=100000000
	contador=0
	write(1,*)inf
	write(1,*)1,A,B
	do i=2,inf
		A=A+(1._4/i)
		B=B+(1._8/i)
		contador=contador+1
		if (contador.eq.1000) then
			write(1,*)i,A,B
			contador=0
		end if
	end do 

!	!PREGUNTA 2
!	write(*,*)"ingresar entero 'n'"
!	read(*,*)m
	m=100
	u0=0 ; u1=1
	write(2,*)m
	write(2,*)0,u0
	write(2,*)1,u1
	do j=2,m
		u2=u0+u1
		write(2,*)j,u2
		u0=u1
		u1=u2
	end do
	
	!PREGUNTA 3
	AA=reshape((/(1._8*k,k=1,12)/),(/3,4/))
	BB=reshape((/(1._8*k,k=1,8)/),(/4,2/))
	call matrix_mult(AA,BB,CC)
!	! quitar comentario para mostrar
!	write(*,*)'A'
!	do k=1,3
!		write(*,*)AA(k,:)
!	end do
!	write(*,*)
!	write(*,*)'B'
!	do k=1,4
!		write(*,*)BB(k,:)
!	end do
!	write(*,*)
!	do k=1,3
!		write(*,*)CC(k,:)
!	end do
	

	!-------------------------------------------
	call system("gnuplot plot.txt")
	!-------------------------------------------
end program

subroutine matrix_mult(A,B,C)
	!entrada
	real(kind=8),dimension(:,:),intent(in)::A,B
	!entrada
	real(kind=8),dimension(:,:),intent(out)::C
	!variables locales
	integer::i,j
	integer,dimension(2)::dim_a,dim_b,dim_c
	real(kind=8)::tmp
	!---------------------
	dim_a=shape(A)
	dim_b=shape(B)
	if ( dim_a(2) .ne.dim_b(1) ) then
		write(*,*) 'dimensiones inconsistentes'
		return
	end if
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
	!---------------------
end subroutine
