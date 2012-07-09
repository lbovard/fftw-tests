program test_1d_spectral            ! fftw3.3 1d r2hc and hc2r test code
      use, intrinsic :: iso_c_binding
      implicit none 
      include 'fftw3.f03'      

      real(C_DOUBLE), parameter :: pi = 3.1415926535897932384626433D+0
      integer, parameter :: N = 128 

      type(C_PTR) :: plan_r2hc,plan_hc2r,pin,pout,pu,puprime
      real(C_DOUBLE), pointer :: in(:),out(:)
      real(C_DOUBLE) :: dx,xlen,t1
      real(C_DOUBLE), dimension(N) :: x,fmx,u,uprime
      integer :: i
      
      pin = fftw_alloc_real(int(N,C_SIZE_T))
      call c_f_pointer(pin,in,[N])
      pout = fftw_alloc_real(int(N,C_SIZE_T))
      call c_f_pointer(pout,out,[N])

      write(6,*) 'planning FFTW_R2HC with FFTW_PATIENT'
      plan_r2hc = fftw_plan_r2r_1d(N,in,out,FFTW_R2HC,FFTW_PATIENT) 
      write(6,*) 'planning FFTW_HC2R with FFTW_PATIENT'
      plan_hc2r = fftw_plan_r2r_1d(N,out,in,FFTW_HC2R,FFTW_PATIENT)
 
      xlen = 2.0D+0*pi
      dx = xlen/DBLE(N)
      do i = 1,N
         x(i) = DBLE(i-1)*dx    ! DBLE(i)*dx should also work
         u(i)      = DEXP(DSIN(x(i)))
         uprime(i) = DCOS(x(i))*u(i)
         in(i)     = u(i)
      end do
      fmx(0+1) = 0.0D+0
      do i = 1,N/2-1
         fmx(i  +1) = -2.0D+0*pi*DBLE(i)/xlen
         fmx(N-i+1) =  2.0D+0*pi*DBLE(i)/xlen
      end do
      fmx(N/2+1) = 0.0D+0
      
      write(6,*) 'r2hc transform'
      call fftw_execute_r2r(plan_r2hc,in,out)
      
      write(6,*) 'd/dx via fourier multiples'
      in(1) = fmx(1)*out(1) 
      do i = 1,N-1
        in(i+1) = fmx(i+1)*out(N-i+1)
      end do
 
      write(6,*) 'hc2r transform'
      call fftw_execute_r2r(plan_hc2r,in,out)
      t1 = 1.0D+0/DBLE(N)
      out(1:N) = t1*out(1:N)
 
      write(6,*) ' '
      write(6,*) 'N = ',N,'  error = ',MAXVAL(DABS(out(1:N)-uprime(1:N)))
      write(6,*) ' '
      write(6,*) 'cleaning up'
      call fftw_free(pin)
      call fftw_free(pout)
      call fftw_destroy_plan(plan_r2hc)
      call fftw_destroy_plan(plan_hc2r)
      
      stop
      end

