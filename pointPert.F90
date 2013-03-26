
module point_module
use common_module
implicit none

integer,private,parameter::ndim=2

contains

  subroutine pointFunction(nx,ny,origin,destination,Ax,Ay,der)
  implicit none
  integer,intent(in)::nx,ny
  real(kind=rk),dimension(ndim),intent(in)::origin,destination
  real(kind=rk),dimension(nx,ny),intent(out)::Ax,Ay
  character(len=*),optional,intent(in)::der
  
  real(kind=rk),dimension(nx)::x
  real(kind=rk),dimension(ny)::y
  integer::i,j
  real(kind=rk)::hx,hy
  integer::dtx,dty

  if(origin(1).le.0..or.origin(1).ge.1. .or. &
     destination(1).le.0..or.destination(1).ge.1. .or. &
     origin(2).le.0..or.origin(2).ge.1. .or. &
     destination(2).le.0..or.destination(2).ge.1.)then 
    call crash('pointFunction: origin and destination must be in (0,1)')
  endif
  
  if(origin(1).ne..5.or.origin(2).ne..5)then
    call crash('pointFunction: origin must be .5')
  endif
  
  if(.not.present(der))then
    dtx=0
    dty=0
  elseif(trim(der).eq.'')then
    dtx=0
    dty=0
  elseif(trim(der).eq.'x')then
    dtx=1
    dty=0
  elseif(trim(der).eq.'y')then
    dtx=0
    dty=1
  elseif(trim(der).eq.'xx')then
    dtx=2
    dty=0
  elseif(trim(der).eq.'xy')then
    dtx=1
    dty=1
  elseif(trim(der).eq.'yy')then
    dtx=0
    dty=2
  else
    call crash('pointFunction: invalid der')
  endif

  hx=domainResolution(nx,'x')
  hy=domainResolution(ny,'y')
  if(dtx.eq.0)then
  do i=1,nx
    x(i)=pert_1d(origin(1),(i-1)*hx)
  enddo
  elseif(dtx.eq.1)then
  do i=1,nx
    x(i)=pert_1d_d1(origin(1),(i-1)*hx)
  enddo
  elseif(dtx.eq.2)then
  do i=1,nx
    x(i)=pert_1d_d2(origin(1),(i-1)*hx)
  enddo
  endif

  if(dty.eq.0)then
  do i=1,ny
    y(i)=pert_1d(origin(2),(i-1)*hy)
  enddo
  elseif(dty.eq.1)then
  do i=1,ny
    y(i)=pert_1d_d1(origin(2),(i-1)*hy)
  enddo
  elseif(dty.eq.2)then
  do i=1,ny
    y(i)=pert_1d_d2(origin(2),(i-1)*hy)
  enddo
  endif
  
  do j=1,ny
    do i=1,nx
      Ax(i,j)=(destination(1)-origin(1))*x(i)*y(j)
      Ay(i,j)=(destination(2)-origin(2))*x(i)*y(j)
    enddo
  enddo

  end subroutine pointFunction
  
  real(kind=rk) function pert_1d(orig,x0)
  real(kind=rk),intent(in)::orig,x0
  real(kind=rk)::t,a,b,c,d,x
  x=2.*abs(x0-.5)
  pert_1d=(2*x-3)*x*x+1
  end function pert_1d
  
  real(kind=rk) function pert_1d_d1(orig,x0)
  real(kind=rk),intent(in)::orig,x0
  real(kind=rk)::t,a,b,c,d,x
  integer::i
  x=2.*abs(x0-.5)
  if(x0.lt..5)then
    t=-2.
  else
    t=2.
  endif
  pert_1d_d1=6*t*(x-1)*x
  end function pert_1d_d1

  real(kind=rk) function pert_1d_d2(orig,x0)
  real(kind=rk),intent(in)::orig,x0
  real(kind=rk)::t,a,b,c,d,x
  integer::i
  x=2.*abs(x0-.5)
  pert_1d_d2=4*(12*x-6)
  end function pert_1d_d2

end module point_module
