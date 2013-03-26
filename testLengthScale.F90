program testLengthScale
use common_module
use fileio_module
use sparse_module
use lengthScale_module
use pointMorphing_module
use correlation_module
use solver_module
implicit none

type(morphing)::m

! the domain size
integer,parameter::nx=500,ny=500

! the point morphing perturbation (s cannot be changed)
! components of d must be in (0.,1.)
real(kind=rk),parameter,dimension(2)::s=(/.5_rk,.5_rk/), &
                                      d=(/.7_rk,.5_rk/)

! correlation parameters                                    
real(kind=rk),parameter::a=-.01_rk,b=-.01_rk,c=0.0_rk

type(ls_coefs)::coef
type(correlation)::corr
type(sparseType)::As
integer::maxnnz
real(kind=rk),dimension(nx,ny,2)::rhs
real(kind=rk),dimension(nx*ny*2)::x

call initialize('testLengthScale.nc',nx,ny)
maxnnz=estimateNNZ(nx,ny)
call initializeMatrix(maxNNZ,As)
call createCorrelation(corr,a,b,c)
call createMorphing(m,nx,ny,s,d)

call write_array(nx,ny,m%Dx,'Dx')
call write_array(nx,ny,m%Dy,'Dy')
call write_array(nx,ny,m%Dxdx,'Dxdx')
call write_array(nx,ny,m%Dydx,'Dydx')
call write_array(nx,ny,m%Dxdy,'Dxdy')
call write_array(nx,ny,m%Dydy,'Dydy')
call write_array(nx,ny,m%Dxdxx,'Dxdxx')
call write_array(nx,ny,m%Dydxx,'Dydxx')
call write_array(nx,ny,m%Dxdyx,'Dxdyx')
call write_array(nx,ny,m%Dydyx,'Dydyx')
call write_array(nx,ny,m%Dxdxy,'Dxdxy')
call write_array(nx,ny,m%Dydxy,'Dydxy')
call write_array(nx,ny,m%Dxdyy,'Dxdyy')
call write_array(nx,ny,m%Dydyy,'Dydyy')

call create_ls_coefs(coef,m)
call compute_ls_coefs(coef,m,corr)

call write_array(nx,ny,coef%beta1,'beta1')
call write_array(nx,ny,coef%beta2,'beta2')
call write_array(nx,ny,coef%detJd,'detJd')

call constructA(nx,ny,coef,corr,As,rhs)

call write_sparse_array(As,'A')
call write_array(nx*ny*2,reshape(rhs,(/nx*ny*2/)),'b')

call solve(nx*ny*2,As,reshape(rhs,(/nx*ny*2/)),x)
call write_array(nx*ny*2,x,'x')

call destroy_ls_coefs(coef)
call destroyMorphing(m)
call destroyCorrelation(corr)
call destroyMatrix(As)
call finalize()
end program testLengthScale
