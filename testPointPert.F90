program testPointPert
use point_module
use fileio_module
use common_module
implicit none

integer,parameter::nx=1000,ny=500
real(kind=rk),dimension(nx,ny)::Ax,Ay
real(kind=rk),dimension(2)::orig,dest

orig=(/.5,.5/)
dest=(/.7,.7/)
call initialize('testPoint.nc',nx,ny)
call pointFunction(nx,ny,orig,dest,Ax,Ay,'')
call write_array(nx,ny,Ax,'Dx')
call write_array(nx,ny,Ay,'Dy')
call pointFunction(nx,ny,orig,dest,Ax,Ay,'x')
call write_array(nx,ny,Ax,'Dxdx')
call write_array(nx,ny,Ay,'Dydx')
call pointFunction(nx,ny,orig,dest,Ax,Ay,'y')
call write_array(nx,ny,Ax,'Dxdy')
call write_array(nx,ny,Ay,'Dydy')
call pointFunction(nx,ny,orig,dest,Ax,Ay,'xx')
call write_array(nx,ny,Ax,'Dxdxx')
call write_array(nx,ny,Ay,'Dydxx')
call pointFunction(nx,ny,orig,dest,Ax,Ay,'xy')
call write_array(nx,ny,Ax,'Dxdxy')
call write_array(nx,ny,Ay,'Dydxy')
call pointFunction(nx,ny,orig,dest,Ax,Ay,'yy')
call write_array(nx,ny,Ax,'Dxdyy')
call write_array(nx,ny,Ay,'Dydyy')
call finalize()

end program testPointPert
