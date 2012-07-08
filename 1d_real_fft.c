/*
 Computes the first spectral derivative of a real function 
 see Spectral Methods in Matlab pg 25
*/
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fftw3.h>

#define pi 3.141592653589793
#define max(x,y) ((x) > (y) ? x : y)

void print_fft(fftw_complex *matrix, int N) {
        int i;
        for(i=0;i<N/2+1;i++) {
		if(matrix[i]<0)
		printf("%lf%lfi\t", matrix[i][0],matrix[i][1]);
		else
		printf("%lf+%lfi\t", matrix[i][0],matrix[i][1]);
                printf("\n");
        }
}
int main() {
        int i,N=8;
        /* k co-efficients only needs to be N/2+1 sized array because of symmetry*/ 
        double v[N],x[N],k[N/2+1],vd[N];
        double dx=2*pi/N,tmp=0;
        fftw_complex *out;
        fftw_plan forward,inverse;
        /* allocate fft array which, because of symmetry is only N/2+1 elements long */
        out= (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N/2+1)); 
        /* r2c_1d always is forward, while c2r_1d is always inverse */
        forward=fftw_plan_dft_r2c_1d(N,v,out,FFTW_ESTIMATE);
        inverse=fftw_plan_dft_c2r_1d(N,out,vd,FFTW_ESTIMATE);
 
        /* initialise arrays */
        for(i=0;i<N;i++) {          
                x[i]=dx*(i+1);
                v[i]=max(0,1-0.5*fabs(x[i]-pi));
		printf("%lf ", v[i]);
        }   
	printf("\n");
        /* initialise co-efficient array */ 
        for(i=0;i<N;i++) {
                        k[i]=i;
        }
        k[N/2]=0;
        /* compute fft */ 
        fftw_execute(forward); 
	print_fft(out,N);
        /* multiply by ik */ 
        for(i=0;i<(N/2+1);i++) {
                tmp=out[i][0];
                out[i][0]=-k[i]*out[i][1];
                out[i][1]=k[i]*tmp;
        }
        /* compute ifft*/
        fftw_execute(inverse);
        fftw_destroy_plan(forward);
        fftw_destroy_plan(inverse);
        fftw_free(out);
        return 0;
      
}

