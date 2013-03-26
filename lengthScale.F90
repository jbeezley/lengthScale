module lengthScale_module
use common_module
use pointMorphing_module
use correlation_module
use sparse_module
implicit none

type ls_coefs
  real(kind=rk),dimension(:,:),pointer::beta1=>null(), &
                                        beta2=>null(), &
                                        detJd=>null()
end type ls_coefs

contains

subroutine create_ls_coefs(coef,m)
implicit none
type(ls_coefs),intent(inout)::coef
type(morphing),intent(in)::m

integer::nx,ny
nx=m%nx
ny=m%ny
allocate(coef%beta1(nx,ny),coef%beta2(nx,ny),coef%detJd(nx,ny))
end subroutine create_ls_coefs

subroutine destroy_ls_coefs(coef)
implicit none
type(ls_coefs),intent(inout)::coef
deallocate(coef%beta1,coef%beta2,coef%detJd)
end subroutine destroy_ls_coefs

subroutine compute_ls_coefs(coef,m,cr)
implicit none
type(ls_coefs),intent(inout)::coef
type(morphing),intent(in)::m
type(correlation),intent(in)::cr

real(kind=rk)::a,b,c,g11,g12,g22,d1g11,d2g11,d1g12,d2g12,d1g22,d2g22
real(kind=rk)::dxdx,dxdy,dydx,dydy
real(kind=rk)::dxdxx,dxdxy,dxdyx,dxdyy
real(kind=rk)::dydxx,dydxy,dydyx,dydyy
integer::i,j,nx,ny

a=cr%a
b=cr%b
c=cr%c

nx=m%nx
ny=m%ny

do j=1,ny
  do i=1,nx
    dxdx=m%Dxdx(i,j) + 1._rk
    dxdy=m%Dxdy(i,j)
    dydx=m%Dydx(i,j)
    dydy=m%Dydy(i,j) + 1._rk

    dxdxx=m%Dxdxx(i,j)
    dxdxy=m%Dxdxy(i,j)
    dxdyx=m%Dxdyx(i,j)
    dxdyy=m%Dxdyy(i,j)

    dydxx=m%Dydxx(i,j)
    dydxy=m%Dydxy(i,j)
    dydyx=m%Dydyx(i,j)
    dydyy=m%Dydyy(i,j)

!    g11=a*dxdx*dxdx+b*dydx*dydx+2.*c*dxdx*dydx
!    g22=a*dxdy*dxdy+b*dydy*dydy+2.*c*dxdy*dydy
!    g12=a*dxdx*dxdy+b*dydx*dydy+c*(dydy*dxdx+dxdy*dydx)

    d1g11=2.*a*dxdx*dxdxx+2.*b*dydx*dydxx+2.*c*(dxdxx*dydx+dxdx*dydxx)
    d2g11=2.*a*dxdx*dxdxy+2.*b*dydx*dydxy+2.*c*(dxdxy*dydx+dxdx*dydxy)
    
    d1g22=2.*a*dxdy*dxdxy+2.*b*dydy*dydxy+2.*c*(dxdxy*dydy+dxdy*dydxy)
    d2g22=2.*a*dxdy*dxdyy+2.*b*dydy*dydyy+2.*c*(dxdyy*dydy+dxdy*dydyy)

    d1g12=a*(dxdxx*dxdy+dxdx*dxdxy)+b*(dydxx+dydx*dydxy)+ &
          c*(dydxy*dxdx+dydy*dxdxx+dxdxy*dydx+dxdy*dydxx)
    d2g12=a*(dxdxy*dxdy+dxdx*dxdyy)+b*(dydxy+dydx*dydyy)+ &
          c*(dydyy*dxdx+dydy*dxdxy+dxdyy*dydx+dxdy*dydxy)

    coef%detJd(i,j)=dxdx*dydy-dxdy*dydx
    coef%beta1(i,j)= d1g11+2.*d2g12-d1g22
    coef%beta2(i,j)=-d2g11+2.*d1g12+d2g22
  enddo
enddo

end subroutine compute_ls_coefs

pure integer function estimateNNZ(nx,ny)
implicit none
integer,intent(in)::nx,ny
estimateNNZ=nx*ny*20
end function estimateNNZ

subroutine constructA(nx,ny,coef,cr,A,rhs)
implicit none
integer,intent(in)::nx,ny
type(ls_coefs),intent(in)::coef
type(correlation),intent(in)::cr
type(sparseType),intent(inout)::A
real(kind=rk),dimension(nx,ny,2),intent(out)::rhs

integer::i,j,k,l,idx1,idy1
real(kind=rk)::b1,b2,gamma,ac,bc,cc,detJd,hx,hy,ihx,ihy,ihx2,ihy2,r1,r2

ac=cr%a
bc=cr%b
cc=cr%c

hx=domainResolution(nx,'x')
hy=domainResolution(ny,'y')
ihx=one/(two*hx)
ihy=one/(two*hy)
ihx2=one/(hx*hx)
ihy2=one/(hy*hy)

idx1=0
idy1=ny*nx

k=0
do j=1,ny
  do i=1,nx
    k=k+1

    ! spatial dependent coefficients
    b1=coef%beta1(i,j)
    b2=coef%beta2(i,j)
    detJd=coef%detJd(i,j)
    gamma=-one/( two*(ac*bc-cc*cc)*detJd )
    
    ! right hand side from substitution:
    !   Dx <- Dx + x
    !   Dy <- Dy + y
    rhs(i,j,1)=-gamma * ( b1*bc - b2*cc )
    rhs(i,j,2)=-gamma * (-b1*cc + b2*ac )
    
    ! diagonal terms
    r1=-two*ihx2-two*ihy2
    call fastset(A,k+idx1,k+idx1,r1)
    call fastset(A,k+idy1,k+idy1,r1)

    if(j.gt.1)then
      l=ravel_multi_index(i,j-1,nx,ny)

      call fastset(A,l+idx1,k+idx1,ihy2-gamma* b1*cc*ihy)
      call fastset(A,l+idy1,k+idx1,-gamma* b1*bc*ihy)
      call fastset(A,l+idx1,k+idy1, gamma* b1*ac*ihy)
      call fastset(A,l+idy1,k+idy1,ihy2+gamma* b1*cc*ihy)
    endif

    if(j.lt.ny)then
      l=ravel_multi_index(i,j+1,nx,ny)

      call fastset(A,l+idx1,k+idx1,ihy2+gamma* b1*cc*ihy)
      call fastset(A,l+idy1,k+idx1, gamma* b1*bc*ihy)
      call fastset(A,l+idx1,k+idy1,-gamma* b1*ac*ihy)
      call fastset(A,l+idy1,k+idy1,ihy2-gamma* b1*cc*ihy)
    endif

    if(i.gt.1)then
      l=ravel_multi_index(i-1,j,nx,ny)

      call fastset(A,l+idx1,k+idx1,ihx2+gamma* b2*cc*ihx)
      call fastset(A,l+idy1,k+idx1, gamma* b2*bc*ihx)
      call fastset(A,l+idx1,k+idy1,-gamma* b2*ac*ihx)
      call fastset(A,l+idy1,k+idy1,ihx2-gamma* b2*cc*ihx)
    endif

    if(i.lt.nx)then
      l=ravel_multi_index(i+1,j,nx,ny)

      call fastset(A,l+idx1,k+idx1,ihx2-gamma* b2*cc*ihx)
      call fastset(A,l+idy1,k+idx1,-gamma* b2*bc*ihx)
      call fastset(A,l+idx1,k+idy1, gamma* b2*ac*ihx)
      call fastset(A,l+idy1,k+idy1,ihx2+gamma* b2*cc*ihx)
    endif
  enddo
enddo

end subroutine constructA

end module lengthScale_module
