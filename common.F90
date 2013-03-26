
module common_module
implicit none

integer,parameter::rk=selected_real_kind(12)
real(kind=rk)::one=1.0_rk,zero=0.0_rk,two=2.0_rk

contains

subroutine crash(msg)
implicit none
character(len=*),intent(in)::msg
write(*,*) trim(msg)
call abort()
end subroutine crash

real(kind=rk) function domainResolution(n,axis)
implicit none
integer,intent(in)::n
character(len=*),optional,intent(in)::axis
domainResolution=1.0_rk/dble(n-1)
end function domainResolution

pure integer function ravel_multi_index(i,j,nx,ny)
integer,intent(in)::i,j,nx,ny
ravel_multi_index=i+(j-1)*nx
end function ravel_multi_index

end module common_module
