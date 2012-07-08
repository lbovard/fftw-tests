program test
implicit none
include 'fftw3.f'

integer :: n, i
integer*8 :: plan
double complex, allocatable :: in(:), out(:)

write(*,*) 'Input data:'
n = 10
allocate(in(n))
allocate(out(n))
do i=1,n
  in(i) = cmplx(i,0.0)
  write(*,*) in(i)
enddo

! Forward Fourier transform
call dfftw_plan_dft_1d(plan,n,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)

write(*,*) 'Fourier transform of the input data:'
do i=1,n
  write(*,*) out(i)
enddo

! Inverse Fourier transform
call dfftw_plan_dft_1d(plan,n,out,in,FFTW_BACKWARD,FFTW_ESTIMATE)
call dfftw_execute(plan)
call dfftw_destroy_plan(plan)

write(*,*) 'Recovered data from inverse Fourier transform:'
do i=1,n
  write(*,*) real(in(i)/n)
enddo

end program test
