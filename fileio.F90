
module fileio_module
use netcdf
use common_module
use sparse_module
implicit none
private

public::initialize,finalize,write_array,write_sparse_array

integer,parameter::maxstring=64
character(len=*),parameter::xdim='x',ydim='y',ndim='N'
integer,parameter::xtype=nf90_double
character(len=maxstring),save::s_filename
integer,save::ncid,xid,yid,s_nx,s_ny,nid

interface write_array
  module procedure write_array_1,write_array_2
end interface

contains

  subroutine initialize(filename,nx,ny)
  implicit none
  character(len=*),intent(in)::filename
  integer,intent(in)::nx,ny
  call check(nf90_create(filename,nf90_clobber,ncid))
  call check(nf90_def_dim(ncid,xdim,nx,xid))
  call check(nf90_def_dim(ncid,ydim,ny,yid))
  call check(nf90_def_dim(ncid,ndim,nx*ny*2,nid))
  s_filename=filename
  s_nx=nx
  s_ny=ny
  call check(nf90_enddef(ncid))
  end subroutine initialize
  
  subroutine write_array_2(nx,ny,A,aname)
  implicit none
  integer,intent(in)::nx,ny
  real(kind=rk),dimension(nx,ny),intent(in)::A
  character(len=*),intent(in)::aname
  integer::varid
  if(nx.ne.s_nx.or.ny.ne.s_ny)then
    call crash('invalid dimensions of '//trim(aname))
  endif
  call check(nf90_sync(ncid))
  call check(nf90_redef(ncid))
  call check(nf90_def_var(ncid,trim(aname),xtype,(/xid,yid/),varid))
  call check(nf90_enddef(ncid))
  call check(nf90_put_var(ncid,varid,A))
  end subroutine write_array_2
  
  subroutine write_array_1(n,A,aname)
  implicit none
  integer,intent(in)::n
  real(kind=rk),dimension(n),intent(in)::A
  character(len=*),intent(in)::aname
  integer::varid
  if(n.ne.s_nx*s_ny*2)then 
    call crash('invalid dimension of '//trim(aname))
  endif
  call check(nf90_sync(ncid))
  call check(nf90_redef(ncid))
  call check(nf90_def_var(ncid,trim(aname),xtype,(/nid/),varid))
  call check(nf90_enddef(ncid))
  call check(nf90_put_var(ncid,varid,A))
  end subroutine write_array_1

  subroutine write_sparse_array(A,aname)
  implicit none
  type(sparseType),intent(in)::A
  character(len=*),intent(in)::aname
  integer::dimid,rowid,colid,dataid
  call check(nf90_sync(ncid))
  call check(nf90_redef(ncid))
  call check(nf90_def_dim(ncid,trim(aname)//'_NNZ',A%nnz,dimid))
  call check(nf90_def_var(ncid,trim(aname)//'_Row',nf90_int,(/dimid/),rowid))
  call check(nf90_def_var(ncid,trim(aname)//'_Col',nf90_int,(/dimid/),colid))
  call check(nf90_def_var(ncid,trim(aname)//'_Dat',xtype,(/dimid/),dataid))
  call check(nf90_put_att(ncid,rowid,'sparse',1))
  call check(nf90_put_att(ncid,colid,'sparse',1))
  call check(nf90_put_att(ncid,dataid,'sparse',1))
  call check(nf90_enddef(ncid))
  call check(nf90_put_var(ncid,rowid,A%row(1:A%nnz)))
  call check(nf90_put_var(ncid,colid,A%col(1:A%nnz)))
  call check(nf90_put_var(ncid,dataid,A%val(1:A%nnz)))
  end subroutine write_sparse_array

  subroutine finalize()
  implicit none
  call check(nf90_sync(ncid))
  call check(nf90_close(ncid))
  s_filename=''
  ncid=-1
  xid=-1
  yid=-1
  s_nx=-1
  s_ny=-1
  end subroutine finalize

  subroutine check(ierr)
  implicit none
  integer,intent(in)::ierr
  if(ierr /= nf90_noerr) then
    call crash(trim(nf90_strerror(ierr)))
  endif
  end subroutine check
end module fileio_module
