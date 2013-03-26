module correlation_module
use common_module
implicit none

type correlation
  real(kind=rk)::a,b,c
end type correlation

contains

subroutine createCorrelation(CR,a,b,c)
implicit none
type(correlation),intent(inout)::CR
real(kind=rk),intent(in)::a,b,c
CR%a=a
CR%b=b
CR%c=c
end subroutine createCorrelation

subroutine destroyCorrelation(CR)
implicit none
type(correlation),intent(inout)::CR
CR%a=0
CR%b=0
CR%c=0
end subroutine destroyCorrelation

end module correlation_module
