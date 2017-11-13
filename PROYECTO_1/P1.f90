program EJERCICIOS_EN_FORTRAN_1
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
	real(kind=4)::A !simple precision
	real(kind=8)::B !doble precision
	character(len=20)::datos_p1='./datos/datos_p1.dat'
	
	! PREGUNTA 2
	integer::n,u0,u1,u2,j
	character(len=20)::datos_p2='./datos/datos_p2.dat'
	
	!PREGUNTA 3
	integer::k
	real(kind=8)::AA(3,4),BB(4,2),CC(3,2)

	!-------------------------------------------
 
	!PREGUNTA 1
	A=1._4
	B=1._8
	inf=10000000
	contador=0
	
	!los datos generados los guarda en el
	!directorio ./datos/datos_p1.dat
	open(unit=10,file=datos_p1,action='write')
	write(10,*)inf
	write(10,*)1,A,B
	do i=2,inf
	
		!contador
		contador=contador+1
	
		!desarrollo de la serie
		A=A+(1._4/i)
		B=B+(1._8/i)
		
		!cada 1000 operaciones exportar datos
		if (contador.eq.1000) then
			write(10,*)i,A,B
			contador=0
		end if
		
	end do 
	close(unit=10)
	
	!-------------------------------------------

!	!PREGUNTA 2
	n=100
	u0=0 ; u1=1
	open(unit=10,file=datos_p2,action='write')
	write(10,*)n
	write(10,*)0,u0
	write(10,*)1,u1
	do j=2,n
		u2=u0+u1
		write(10,*)j,u2
		u0=u1
		u1=u2
	end do
	close(unit=10)
		
	!-------------------------------------------
	
	!PREGUNTA 3
	!1. se multiplican dos matrices consistentes
	AA=reshape((/(1._8*k,k=1,12)/),(/3,4/))
	BB=reshape((/(1._8*k,k=1,8)/),(/4,2/))
	call matrix_mult(AA,BB,CC)
	
	!2.se multiplican dos matrices no consistentes
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
	
end program

