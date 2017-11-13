program lorenz

!Parametros

implicit none

integer(kind=4),parameter::m=3
integer(kind=4),parameter::n=200000
character(len=255) command_filename
integer(kind=4) command_unit
character(len=255) data_filename
integer(kind=4) data_unit
real(kind=8) dt
integer(kind=4) i
integer(kind=4) j
external lorenz_rhs
real(kind=8)t(0:n)
real(kind=8)t_final
real(kind=8)x(m,0:n)

call timestamp ( )
write (*,'(a)') ''
write (*,'(a)') 'LORENZ_ODE'
write (*,'(a)') 'FORTRAN90'
write (*,'(a)') 'Calcular las soluciones del Sistema de Lorenz.'
  write (*,'(a)') 'Escribir data a un file para usar en Gnuplot.'

!Data

t_final=20.0D+00
dt=t_final/real(n,kind=8)

!Condiciones iniciales

do j=0,n
   t(j)=real(j,kind=8)*t_final/real(n,kind=8)
end do

x(1:m,0)=(/ 8.0D+00, 1.0D+00, 1.0D+00 /)

!Calcular la aproximacion de la solicion en iguales espacion de tiempo

do j=0, n-1
call rk4vec (t(j),m,x(1:m,j),dt,lorenz_rhs,x(1:m,j+1))
end do

!Crear el data file

call get_unit (data_unit)
data_filename = 'lorenz_edo_data.txt'
open (unit = data_unit, file = data_filename,status= 'replace')
do j=0, n, 50
write (data_unit, '(2x,g14.6,2x,g14.6,2x,g14.6,2x,g14.6)' ) t(j),x(1:m,j)
end do
 close (unit=data_unit)
write (*,'(a)') 'Crear un data file "' // trim (data_filename)// '".'

!Crear el comando file 

call get_unit (command_unit)
command_filename='lorenz_data.txt'
open (unit=command_unit, file = command_filename, status='replace')
write (command_unit,'(a)') '#' // trim (command_filename)
write (command_unit,'(a)') '#'
write (command_unit,'(a)') '# Usage:'
write (command_unit,'(a)') '# gnuplot < ' // trim ( command_filename )
write (command_unit,'(a)') '#'
write (command_unit,'(a)') 'set term png'
write (command_unit,'(a)') &
    'set output "xyz_time.png"'
write (command_unit,'(a)') 'set xlabel "<--- T --->"'
write (command_unit,'(a)') 'set ylabel "<--- X(T), Y(T), Z(T) --->"'
write (command_unit,'(a)') &
'set title "X(T), Y(T), Z(T) versus Time"'
write (command_unit,'(a)') 'set grid'
write (command_unit,'(a)') 'set style data lines'
write (command_unit,'(a)') 'plot "' // trim (data_filename) // &
'" using 1:2 lw 3 linecolor rgb "blue",' // &
' "" using 1:3 lw 3 linecolor rgb "red",' // &
' "" using 1:4 lw 3 linecolor rgb "green"'
write (command_unit,'(a)') &
'set output "xyz_3d.png"'
write (command_unit, '(a)') 'set xlabel "<--- X(T) --->"'
write (command_unit,'(a)') 'set ylabel "<--- Y(T) --->"'
write (command_unit,'(a)') 'set zlabel "<--- Z(T) --->"'
write (command_unit,'(a)') &
'set title "(X(T),Y(T),Z(T)) trajectory"'
write (command_unit,'(a)') 'set grid'
write (command_unit,'(a)') 'set style data lines'
write (command_unit,'(a)') 'splot "' // trim (data_filename) // &
'" using 2:3:4 lw 1 linecolor rgb "blue"'
 close (unit=command_unit)
write (*,'(a)') &
' Crear un command file "' // trim (command_filename) // '".'

!Terminar

write (*,'(a)') ''
write (*,'(a)') 'LORENZ_ODE:'
write (*,'(a)') 'Ejecucion normal.'
write (*,'(a)') ''
call timestamp ( )
stop
end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine get_unit ( iunit )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
integer(kind=4) i
integer(kind=4) ios
integer(kind=4) iunit
logical lopen
  iunit = 0
  do i = 1, 99
    if (i/=5.and.i/=6.and.i/=9)then
      inquire (unit=i,opened=lopen,iostat=ios)
      if(ios==0)then
        if(.not.lopen)then
          iunit=i
          return
        end if
      end if
    end if
  end do
return
end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine lorenz_rhs(t,m,x,dxdt)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Input, real el valor de la variable independiente 
!Input, integer en la dimension espacio
!Input, Valores de las variables en el tiempo T
!Output, Valores de las derivadas de las variables dependientes en el tiempo
implicit none
integer(kind=4) m
real(kind=8),parameter::beta=8.0D+00/3.0D+00
real(kind=8) dxdt(m)
real(kind=8),parameter::rho=28.0D+00
real(kind=8),parameter::sigma=10.0D+00
real(kind=8) t
real(kind=8) x(m)

dxdt(1)=sigma*(x(2)-x(1))
dxdt(2)=x(1)*(rho-x(3))-x(2)
dxdt(3)=x(1)*x(2)-beta*x(3)

return
end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine rk4vec (t0,m,u0,dt,f,u)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!rk4vec toma un paso de RK4 para un vector
!Output Runge-Kutta de 4° orden
!solucion estimada para el tiempo actual

implicit none

integer(kind=4) m
real (kind=8) dt
external f
real (kind=8) f0(m)
real (kind=8) f1(m)
real (kind=8) f2(m)
real (kind=8) f3(m)
real (kind=8) t0
real (kind=8) t1
real (kind=8) t2
real (kind=8) t3
real (kind=8) u(m)
real (kind=8) u0(m)
real (kind=8) u1(m)
real (kind=8) u2(m)
real (kind=8) u3(m)
!
!Obtener cuatro muestras de los valores de las derivadas.
!
call f(t0,m,u0,f0)
t1 = t0 + dt / 2.0D+00
u1(1:m) = u0(1:m) + dt * f0(1:m) / 2.0D+00
call f(t1,m,u1,f1)
t2 = t0 + dt / 2.0D+00
u2(1:m) = u0(1:m) + dt * f1(1:m) / 2.0D+00
call f(t2,m,u2,f2)
t3 = t0 + dt
u3(1:m) = u0(1:m) + dt * f2(1:m)
call f (t1, m, u1, f3)

!Combinarlos para estimar la solucion al tiempo T0+DT

u(1:m)=u0(1:m)+dt*(f0(1:m)+2.0D+00*f1(1:m)+2.0D+00*f2(1:m)+f3(1:m))/6.0D+00

return
end 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine timestamp ( )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none

character(len=8) ampm
integer(kind=4) d
integer(kind=4) h
integer(kind=4) m
integer(kind=4) mm
character(len=9),parameter,dimension(12)::month=(/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
integer(kind=4) n
integer(kind=4) s
integer(kind=4) values(8)
integer(kind=4) y

call date_and_time (values=values)

y=values(1)
m=values(2)
d=values(3)
h=values(5)
n=values(6)
s=values(7)
mm=values(8)
if (h<12) then
  ampm = 'AM'
  else if (h==12) then
  if (n==0.and.s==0) then
           ampm='Noon'
  else
      ampm='PM'
end if
  else
    h=h-12
  if (h<12) then
      ampm='PM'
  else if (h==12) then
  if (n==0.and.s==0) then
        ampm='Midnight'
   else
        ampm='AM'
   end if
 end if
end if

write(*,'(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)' ) &
d,trim (month(m)),y,h,':',n,':',s,'.',mm,trim (ampm)
return
end