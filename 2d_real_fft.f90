program test
        use, intrinsic :: iso_c_binding
        implicit none
        include "fftw3.f03"
        
        type(C_PTR) :: forward, inverse 

        integer, parameter :: N=32
        integer, parameter :: Nr=(N/2)+1
        integer :: i,j 

        real(kind=8), dimension(N,N+2) :: v
        real(kind=8), dimension(:,:), allocatable :: kx,ky

        complex(kind=8), parameter :: ii=(0,1)
        complex(C_DOUBLE_COMPLEX), dimension(Nr,N) :: out

        forward = fftw_plan_dft_r2c_2d(N,N,v,out,FFTW_ESTIMATE)
        inverse = fftw_plan_dft_c2r_2d(N,N,out,v,FFTW_ESTIMATE)
        
        allocate(kx(Nr,N))
        allocate(ky(Nr,N))
        forall(i=1:Nr,j=1:N/2) kx(i,j)=real(j-1,8)
        forall(i=1:Nr,j=N/2+2:N) kx(i,j)=real(j-N-1,8)
        forall(i=1:Nr-1,j=1:N) ky(i,j)=real(i-1,8)
        ky(Nr,:)=0.0_8
        forall(i=1:N,j=1:N) v(i,j)=1.0_8/real(i+j-1,8) 
        do i=1,N-1
                v(i,i+1)=real(i+1,8)
        end do
        do i=1,N
                print '(50f14.6)', v(i,1:N)
        end do
   !     do i=1,Nr
   !             print '(50f10.6)', kx(i,1:N)
   !     end do
   !     print '()'
   !     do i=1,Nr
   !             print '(50f10.6)', ky(i,1:N)
   !     end do
        call fftw_execute_dft_r2c(forward, v, out)      
        out= out*ii*kx
        call fftw_execute_dft_c2r(inverse,out,v)
        print '()'
        forall(i=1:N,j=1:N) v(i,j)=v(i,j)/real(N*N,8)
        do i=1,N
                print '(50f15.6)', v(i,1:N)
        end do
end program test


subroutine print_sq_matrix(A,n)
        implicit none
        integer, intent(in) :: n
        real(kind=8), intent(in) :: A(n,n)
        integer :: i
        do i=1,n
                print '(5f10.6)', A(i,1:n)
        end do
end 

subroutine print_fft_matrix(A,n)
        implicit none
        integer, intent(in) :: n
        real(kind=8), intent(in) :: A((n/2+1),n)
        integer :: i,j
        do i=1,n
                print*, A(i,:)
        end do
end 
