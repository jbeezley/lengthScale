module derivative_module
use common_module
implicit none

contains

subroutine diffx(nx,ny,A,D)
implicit none
integer,intent(in)::nx,ny
real(kind=rk),dimension(nx,ny),intent(in)::A
real(kind=rk),dimension(nx,ny),intent(out)::D

integer::i,j
real(kind=rk)::h,h1,h2
h=domainResolution(nx,'x')
h1=1./h
h2=1./(2.*h)
do j=1,ny
  do i=2,nx-1
    D(i,j)=h2*(A(i-1,j)-A(i+1,j))
  enddo
  D(1 ,j)=h1*(A(2 ,j)-A(1   ,j))
  D(nx,j)=h1*(A(nx,j)-A(nx-1,j))
enddo
end subroutine diffx

subroutine diffy(nx,ny,A,D)
implicit none
integer,intent(in)::nx,ny
real(kind=rk),dimension(nx,ny),intent(in)::A
real(kind=rk),dimension(nx,ny),intent(out)::D

integer::i,j
real(kind=rk)::h,h1,h2
h=domainResolution(ny,'y')
h1=1./h
h2=1./(2.*h)
do j=2,ny-1
  do i=1,nx
    D(i,j)=h2*(A(i,j-1)-A(i,j+1))
  enddo
enddo
do i=1,nx
  D(i,1 )=h1*(A(i, 2)-A(i,   1))
  D(i,ny)=h1*(A(i,ny)-A(i,ny-1))
enddo
end subroutine diffy

subroutine diff2x(nx,ny,A,D)
implicit none
integer,intent(in)::nx,ny
real(kind=rk),dimension(nx,ny),intent(in)::A
real(kind=rk),dimension(nx,ny),intent(out)::D

integer::i,j
real(kind=rk)::h,h2
h=domainResolution(nx,'x')
h2=1./(h*h)
do j=1,ny
  do i=2,nx-1
    D(i,j)=h2*(A(i-1,j)-2.*A(i,j)+A(i+1,j))
  enddo
  D(1 ,j)=D(2,j)
  D(nx,j)=D(nx-1,j)
enddo
end subroutine diff2x

subroutine diff2y(nx,ny,A,D)
implicit none
integer,intent(in)::nx,ny
real(kind=rk),dimension(nx,ny),intent(in)::A
real(kind=rk),dimension(nx,ny),intent(out)::D

integer::i,j
real(kind=rk)::h,h2
h=domainResolution(ny,'y')
h2=1./(h*h)
do j=2,ny-1
  do i=1,nx
    D(i,j)=h2*(A(i,j-1)-2.*A(i,j)+A(i,j+1))
  enddo
enddo
do i=1,nx
  D(i,1 )=D(i,2)
  D(i,ny)=D(i,ny-1)
enddo
end subroutine diff2y

subroutine diffxy(nx,ny,A,D)
implicit none
integer,intent(in)::nx,ny
real(kind=rk),dimension(nx,ny),intent(in)::A
real(kind=rk),dimension(nx,ny),intent(out)::D
real(kind=rk),dimension(nx,ny)::T
call diffx(nx,ny,A,T)
call diffy(nx,ny,T,D)
end subroutine diffxy

end module derivative_module
