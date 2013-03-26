
module image_interp
!use EZspline_obj
use common_module
implicit none
private

!*** module containing stubs to pspline interpolation routines
!    possibly adding problem specific optimizations in the future
!    using low level pspline functions

!*** interpolation object, initialized once before computing values
!    requires deallocation through destroy subroutine
!type image2d_interp_data
!  private
!  type(EZspline2_r8)::spline_o
!  integer::nx,ny
!  real(kind=rk),dimension(:,:),pointer::img
!end type image2d_interp_data
!
!public::image2d_interp_data,init_interp_image2d,destroy_interp, &
!        interp_image2d,grad_image2d,interp_image2d_grid,        &
!        fast_interp
!

public::fast_interp

contains

!subroutine init_interp_image2d(imgd,nx,ny,imgin,llx,lly,urx,ury)
!use EZspline 
!implicit none
!
!!*** initialize/allocate interpolation data structure out of imgin data
!
!type(image2d_interp_data),intent(inout)::imgd       ! data structure
!integer,intent(in)::nx,ny                           ! image size
!real(kind=rk),dimension(nx,ny),target,intent(in)::imgin      ! image data
!real(kind=rk),intent(in)::llx,lly,urx,ury
!
!integer::ier
!integer,dimension(2)::bcond
!
!! todo:  ADD sub image specifications!!!
!bcond=(/0,0/)
!call EZspline_init2_r8(imgd%spline_o,nx,ny,bcond,bcond,ier)
!call interp_error(ier)
!
!call unif_grid_sub(nx,imgd%spline_o%x1,llx,urx)
!call unif_grid_sub(ny,imgd%spline_o%x2,lly,ury)
!
!call EZspline_setup(imgd%spline_o,imgin,ier)
!call interp_error(ier)
!
!imgd%nx=nx
!imgd%ny=ny
!imgd%img=>imgin
!
!end subroutine init_interp_image2d
!
!subroutine destroy_interp(imgd)
!use EZspline 
!implicit none
!
!!*** destroy/deallocate interpolation data structure
!
!type(image2d_interp_data),intent(inout)::imgd
!integer::ier
!call EZspline_free(imgd%spline_o,ier)
!call interp_error(ier)
!
!imgd%nx=-1
!imgd%ny=-1
!nullify(imgd%img)
!
!end subroutine destroy_interp
!
!subroutine interp_image2d(imgd,nxy,x,y,vals)
!use EZspline 
!implicit none
!
!!***  interpolate the image at nxy values
!
!type(image2d_interp_data),intent(in)::imgd  ! interpolation data structure
!                                            ! (should be initialized)
!integer,intent(in)::nxy                     ! number of data points
!real(kind=rk),dimension(nxy),intent(in)::x,y         ! location of requested data
!real(kind=rk),dimension(nxy),intent(out)::vals       ! interpolated data
!
!!*** all location information (x,y) is assumed to be relative to a [0,1]x[0,1]
!!    grid, all (x,y) should be in this interval or an error will be thrown.
!!    Perhaps I will change this behavior for the benefit of warp_image code, as
!!    it stands now data will have to be copied to contiguous arrays.
!
!integer::ier
!if(nxy.lt.1)return
!call EZspline_interp(imgd%spline_o,nxy,x,y,vals,ier)
!call interp_error(ier)
!end subroutine interp_image2d
!
!subroutine interp_image2d_grid(imgd,nx,ny,vals)
!use EZspline
!implicit none
!
!!***  interpolate an image to a structured grid
!
!type(image2d_interp_data),intent(in)::imgd  ! interpolation data structure
!integer,intent(in)::nx,ny                   ! size of output grid in x and y
!real(kind=rk),dimension(nx,ny),intent(out)::vals     ! output values
!
!real(kind=rk),dimension(nx)::x
!real(kind=rk),dimension(ny)::y
!integer::ier
!if(nx.lt.1.or.ny.lt.1)return
!call unif_grid(nx,x)
!call unif_grid(ny,y)
!call EZspline_interp(imgd%spline_o,nx,ny,x,y,vals,ier)
!call interp_error(ier)
!end subroutine interp_image2d_grid
!
!subroutine unif_grid_sub(n,x,sx,ex)
!implicit none
!integer,intent(in)::n
!real(kind=rk),dimension(n),intent(out)::x
!real(kind=rk),intent(in)::sx,ex
!real(kind=rk)::dx
!integer::i
!dx=(ex-sx)/real(n-1)
!x(1)=sx
!do i=2,n-1
!  x(i)=real(i-1)*dx+sx
!enddo
!x(n)=ex
!end subroutine unif_grid_sub
!
!subroutine unif_grid(n,x)
!implicit none
!
!!***  return a uniformly spaced grid i.e. x=linspace(0,1,n) in matlab
!
!integer,intent(in)::n
!real(kind=rk),dimension(n),intent(out)::x
!real(kind=rk)::dx
!integer::i
!dx=1./real(n-1)
!x(1)=0.
!do i=2,n-1
!  x(i)=real(i-1)*dx
!enddo
!x(n)=1.
!end subroutine unif_grid
!
!subroutine grad_image2d(imgd,nxy,x,y,dvals)
!use EZspline 
!implicit none
!
!integer,parameter::two=2
!
!!*** returns the gradient of the image represented by imgd, at (x,y).  
!!    dvals(:,1) contains x-direction partial
!!    dvals(:,2) contains y-direction partial
!
!type(image2d_interp_data),intent(in)::imgd
!integer,intent(in)::nxy
!real(kind=rk),dimension(nxy),intent(in)::x,y
!real(kind=rk),dimension(nxy,two)::dvals
!
!integer::ier
!if(nxy.lt.1)return
!call EZspline_gradient(imgd%spline_o,nxy,x,y,dvals,ier)
!call interp_error(ier)
!end subroutine grad_image2d
!
!subroutine interp_error(ier)
!use EZspline
!implicit none
!integer,intent(in)::ier
!call EZspline_error(ier)
!if(ier.ne.0)then
!  call crash('')
!endif
!end subroutine interp_error

subroutine fast_interp(nx,ny,img,x,y,vals,llx,lly,urx,ury)
implicit none
integer,intent(in)::nx,ny
real(kind=rk),dimension(nx,ny),intent(in)::img
real(kind=rk),dimension(nx,ny),intent(in)::x,y
real(kind=rk),dimension(nx,ny),intent(out)::vals
real(kind=rk),intent(in)::llx,lly,urx,ury
integer::i,j,sx,sy,ex,ey
real(kind=rk)::mx,my
do j=1,ny
  my=(y(i,j)-lly)*real(ny-1)/(ury-lly)+1.
  sy=max(min(floor(my),ny-1),1)
  ey=sy+1
  do i=1,nx
    mx=(x(i,j)-llx)*real(nx-1)/(urx-llx)+1.
    sx=max(min(floor(mx),nx-1),1)
    ex=sx+1
    vals(i,j)=(ex-mx)*(ey-my)*img(sx,sy)+&
              (mx-sx)*(ey-my)*img(ex,sy)+&
              (ex-mx)*(my-sy)*img(sx,ey)+&
              (mx-sx)*(my-sy)*img(ex,ey)
  enddo
enddo
end subroutine fast_interp

end module image_interp
