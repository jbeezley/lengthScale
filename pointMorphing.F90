module pointMorphing_module
use morphing_module , only : deformArray, &
                             setupMorphing, &
                             p_createMorphing=>createMorphing, &
                             p_destroyMorphing=>destroyMorphing, &
                             morphing
use point_module , only : pointFunction
use common_module
implicit none

contains

subroutine createMorphing(m,nx,ny,s,d)
implicit none
type(morphing),intent(inout)::m
integer,intent(in)::nx,ny
real(kind=rk),dimension(2),intent(in)::s,d
allocate(m%Dx(nx,ny),m%Dy(nx,ny))
call setupMorphing(m,nx,ny)
call pointFunction(nx,ny,s,d,m%Dx,m%Dy,'')
call pointFunction(nx,ny,s,d,m%Dxdx,m%Dydx,'x')
call pointFunction(nx,ny,s,d,m%Dxdy,m%Dydy,'y')
call pointFunction(nx,ny,s,d,m%Dxdxx,m%Dydxx,'xx')
call pointFunction(nx,ny,s,d,m%Dxdxy,m%Dydxy,'xy')
call pointFunction(nx,ny,s,d,m%Dxdyy,m%Dydyy,'yy')
m%Dxdyx=>m%Dxdxy
m%Dydyx=>m%Dydxy
end subroutine createMorphing

subroutine destroyMorphing(m)
implicit none
type(morphing),intent(inout)::m
deallocate(m%Dx,m%Dy)
call p_destroyMorphing(m)
end subroutine destroyMorphing

end module pointMorphing_module
