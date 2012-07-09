program test
        use, intrinsic :: iso_c_binding
        implicit none
        include 'fftw3.f03'
        integer, parameter :: N=128
        real(C_DOUBLE), pointer :: y(:,:)
        type(C_PTR) :: p
        p=fftw_alloc_real(int(N,C_SIZE_T))
end program test

