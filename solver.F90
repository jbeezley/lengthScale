module solver_module
use sparse_module
use common_module
implicit none

contains

subroutine solve(n,A,b,x)
implicit none
integer,intent(in)::n
type(sparseType),intent(in)::A
real(kind=rk),dimension(n),intent(in)::b
real(kind=rk),dimension(n),intent(out)::x

real(kind=rk),dimension(A%nnz)::val
integer,dimension(A%nnz)::irow,icol
real(kind=rk),dimension(:),pointer::values
integer,dimension(:),pointer::columns,rows

! solve the system of equations: A*x=b
! 
! where A is a non-symmetric sparse matrix of size n x n
!       b is a dense vector of length n

! Set up arrays of the sparse array A in triplet format
! Non-zero values of A
values=>A%val(1:A%nnz)

! Non-zero columns of A (indexed between 1 and n, fortran style)
columns=>A%col(1:A%nnz)
 
! Non-zero rows of A (indexed between 1 and n, fortran style)
rows=>A%row(1:A%nnz)

! call mumps here and output into 'x'
! for now just set x to zero:
x(:)=0.0_rk
end subroutine solve

end module solver_module
