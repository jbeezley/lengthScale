module morphing_module
use common_module
use image_interp , only : fast_interp
use derivative_module
implicit none

type morphing
  integer::nx,ny
  real(kind=rk),dimension(:,:),pointer::Dx=>null(),Dy=>null()
  real(kind=rk),dimension(:,:),pointer::Dxdx=>null(),Dydx=>null(), &
                                        Dxdy=>null(),Dydy=>null()
  real(kind=rk),dimension(:,:),pointer::Dxdxx=>null(),Dydxx=>null(), &
                                        Dxdxy=>null(),Dydxy=>null(), &
                                        Dxdyx=>null(),Dydyx=>null(), &
                                        Dxdyy=>null(),Dydyy=>null()
end type morphing

contains

subroutine createMorphing(m,nx,ny,Dx,Dy)
implicit none
type(morphing),intent(inout)::m
integer,intent(in)::nx,ny
real(kind=rk),dimension(nx,ny),target,intent(in)::Dx,Dy
call setupMorphing(m,nx,ny)
m%Dx=>Dx
m%Dy=>Dy
call diffx(nx,ny,Dx,m%Dxdx)
call diffx(nx,ny,Dy,m%Dydx)
call diffy(nx,ny,Dx,m%Dxdy)
call diffy(nx,ny,Dy,m%Dydy)
call diff2x(nx,ny,Dx,m%Dxdxx)
call diff2x(nx,ny,Dy,m%Dydxx)
call diffxy(nx,ny,Dx,m%Dxdxy)
call diffxy(nx,ny,Dy,m%Dydxy)
call diff2y(nx,ny,Dx,m%Dxdyy)
call diff2y(nx,ny,Dy,m%Dydyy)
m%Dxdyx=>m%Dxdxy
m%Dydyx=>m%Dydxy
end subroutine createMorphing

subroutine setupMorphing(m,nx,ny)
implicit none
type(morphing),intent(inout)::m
integer,intent(in)::nx,ny
m%nx=nx
m%ny=ny
allocate(m%Dxdx(nx,ny),m%Dydx(nx,ny),m%Dxdy(nx,ny),m%Dydy(nx,ny))
allocate(m%Dxdxx(nx,ny),m%Dydxx(nx,ny),m%Dxdxy(nx,ny),m%Dydxy(nx,ny))
allocate(m%Dxdyy(nx,ny),m%Dydyy(nx,ny))
end subroutine setupMorphing

subroutine deformArray(m,nx,ny,A,B)
implicit none
type(morphing),intent(in)::m
integer,intent(in)::nx,ny
real(kind=rk),dimension(nx,ny),intent(in)::A
real(kind=rk),dimension(nx,ny),intent(out)::B
integer::n
real(kind=rk)::llx,lly,urx,ury
if(nx.ne.m%nx.or.ny.ne.m%ny)then
  call crash('deformArray: invalid input array size')
endif
llx=0.
lly=0.
urx=1.
ury=1.
call fast_interp(nx,ny,A,m%Dx,m%Dy,B,llx,lly,urx,ury)
end subroutine deformArray

subroutine destroyMorphing(m)
implicit none
type(morphing),intent(inout)::m
m%nx=0
m%ny=0
m%Dx=>null()
m%Dy=>null()
deallocate(m%Dxdx,m%Dydx,m%Dxdy,m%Dydy)
deallocate(m%Dxdxx,m%Dydxx,m%Dxdxy,m%Dydxy)
deallocate(m%Dxdyy,m%Dydyy)
m%Dxdyx=>null()
m%Dydyx=>null()
end subroutine destroyMorphing

end module morphing_module
