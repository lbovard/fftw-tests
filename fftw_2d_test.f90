program fftw_test
        use, intrinsic :: iso_c_binding
        implicit none
        include "fftw3.f03"
        type(C_PTR) :: forward, inverse
        !integer, parameter :: M=10, L=10
        integer, parameter :: N=10
        integer, parameter :: Nr=(N/2)+1
        real(C_DOUBLE), pointer :: arr(:,:)
        type(C_PTR) :: p
        real(kind=8), dimension(:,:), allocatable :: kx,ky

        complex(kind=8), parameter :: ii=(0,1)
        complex(C_DOUBLE_COMPLEX), dimension(Nr,N) :: out 
        p = fftw_alloc_real(int(N*N, C_SIZE_T))

        call c_f_pointer(p, arr, [N,N])
        forward = fftw_plan_dft_r2c_2d(N,N,arr,out,FFTW_ESTIMATE)
        inverse = fftw_plan_dft_c2r_2d(N,N,out,arr,FFTW_ESTIMATE)

end program fftw_test
