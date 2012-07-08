program test
        use, intrinsic :: iso_c_binding
        implicit none
        include "fftw3.f03"
        type(C_PTR) :: forward, inverse 
        integer, parameter :: N=4
        integer, parameter :: Nr=(N/2)+1
!        real(kind=8), dimension(N,N) :: v
        real(C_DOUBLE), dimension(N,N) :: v
        complex(kind=8), parameter :: ii=(0,1)
        complex(C_DOUBLE_COMPLEX), dimension(Nr,N) :: out
        real(kind=8), parameter :: pi = 3.1415926535897938_8
        real(kind=8), parameter :: L=30.0_8
        real(kind=8), dimension(:,:), allocatable :: kx,ky
        integer :: i,j 
        forward = fftw_plan_dft_r2c_2d(N,N,v,out,FFTW_ESTIMATE)
        inverse = fftw_plan_dft_c2r_2d(N,N,out,v,FFTW_ESTIMATE)
        allocate(kx(N,Nr))
        allocate(ky(N,Nr))
        forall(j=1:Nr,i=1:N) kx(i,j)=real(j-1,8)
        forall(j=1:Nr,i=1:Nr-1) ky(i,j)=real(i-1,8) !possible oneliner
        forall(j=1:Nr,i=Nr:N) ky(i,j)=real(i-N-1,8)        
        forall(j=1:N,i=1:N) v(i,j)=1.0/real(i+j-1,8)
        do i=1,N-1
                v(i,i+1)=real(i+1,8)
        end do
        v=transpose(v)
        call fftw_execute_dft_r2c(forward, v, out)      
        !out=out*kx*ii
       print *, out
end program test


subroutine print_sq_matrix(A,n)
        implicit none
        integer, intent(in) :: n
        real(kind=8), intent(in) :: A(n,n)
        integer :: i
        do i=1,n
                print*, A(i,:)
        end do
end 

subroutine print_fft_matrix(A,n)
        implicit none
        integer, intent(in) :: n
        real(kind=8), intent(in) :: A(n,(n/2+1))
        integer :: i,j
        do i=1,n
                print*, A(i,:)
        end do
end 
