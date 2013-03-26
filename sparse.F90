module sparse_module
use common_module
!*** a simple (but inefficent) sparse matrix type
implicit none

type sparseType
    !*** derived type defining a sparse matrix
    !  non-zero entries are defined at
    !  A(row(k),col(k)) == val(k) for k=1,...,nnz

    integer::nnz,   &  ! number of non-zeros contained in the matrix
             maxNNZ    ! maximum number of non-zeros possible (defined at
                       ! initialization
    real(kind=rk),pointer,dimension(:) :: val ! array containing the values of
                                              ! non-zeros
    integer,pointer,dimension(:) :: col, & ! array containing column indices of non-zeros
                                    row    ! array containing row indices of non-zeros
end type sparseType

contains

subroutine initializeMatrix(maxNNZ,A)
!*** initialize a sparse matrix type to a given fixed
!    maximum number of non-zeros
!    (this allocates memory that must be deallocated when finished with the
!    object)
implicit none
integer,intent(in)::maxNNZ
type(sparseType),intent(inout)::A
A%maxNNZ=maxNNZ
A%nnz=0
allocate(A%val(maxNNZ),A%col(maxNNZ),A%row(maxNNZ))
end subroutine initializeMatrix

subroutine destroyMatrix(A)
!*** destroy a sparse matrix object and deallocate memory
implicit none
type(sparseType),intent(inout)::A
A%nnz=-1
A%maxNNZ=-1
deallocate(A%val,A%col,A%row)
end subroutine destroyMatrix

subroutine set(A,i,j,val)
!***  set A_ij to val
!
!     slow method that searches for and replaces old values
implicit none
type(sparseType),intent(inout)::A
integer,intent(in)::i,j
real(kind=rk),intent(in)::val

integer::k

k=findIdx(A,i,j)
if(k.lt.1)then
    A%nnz=A%nnz+1
    k=A%nnz
    if(k.gt.A%maxNNZ)then
        call crash('maxNNZ not large enough.')
    endif
    A%col(k)=i
    A%row(k)=j
endif
A%val(k)=val
end subroutine set

subroutine fastset(A,i,j,val)
!*** fast version of set that will just append
!    a new value rather than searching for old values...
!    this should only be used if A_ij is 0 before calling
implicit none
type(sparseType),intent(inout)::A
integer,intent(in)::i,j
real(kind=rk),intent(in)::val

integer::k

A%nnz=A%nnz+1
k=A%nnz
if(k.gt.A%maxNNZ)then
    call crash('maxNNZ not large enough.')
endif
A%col(k)=i
A%row(k)=j
A%val(k)=val
end subroutine fastset

real(kind=rk) function get(A,i,j)
!*** return the value of A_ij
implicit none
type(sparseType),intent(in)::A
integer,intent(in)::i,j
integer::k
k=findIdx(A,i,j)
if(k.lt.1) then
    get=0.
else
    get=A%val(k)
endif
end function get

integer function findIdx(A,i,j)
!*** helper function that searches A%row and A%col
!    for the i==A%col(k) and j==A%row(k)
!    returns k>0 if found or k=-1 if not found
implicit none
type(sparseType),intent(in)::A
integer,intent(in)::i,j
integer::k
do k=1,A%nnz
    if(A%col(k).eq.i.and.A%row(k).eq.j)then
        findIdx=k
        return
    endif
enddo
findIdx=-1
end function findIdx

subroutine printMatrix(A)
!*** write the matrix to 3 seperate (but equally sized) files
!    rows.txt containing the row indices of non-zeros
!    columns.txt containing the column indices of non-zeros
!    values.txt containing the values of non-zeros
implicit none
type(sparseType),intent(in)::A
integer,parameter::iu=26
integer::i
open(unit=iu,file='values.txt')
do i=1,A%nnz
    write(iu,*) A%val(i)
enddo
close(iu)
open(unit=iu,file='rows.txt')
do i=1,A%nnz
    write(iu,*) A%row(i)
enddo
close(iu)
open(unit=iu,file='columns.txt')
do i=1,A%nnz
    write(iu,*) A%col(i)
enddo
close(iu)
end subroutine printMatrix

end module sparse_module
