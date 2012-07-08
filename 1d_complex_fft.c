#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fftw3.h>

int main() {
        /* make an array of 20 elements */
        int N=20;
        int i;
        /* completely random! */ 
        unsigned int seed = 50085;
        srand( seed );
        /* define special fftw complex datatypes */
        fftw_complex *in, *out;
        /* define fftw_plan which contains all the data needed
           to compute the fft */
        fftw_plan forward;
        fftw_plan backward;
        /* allocate memory */
        in = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
        out = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*N);
        
        /* define the plan for a 1d complex dft
           function prototype is 
           fftw_plan_dft_1d(int n, fftw_complex *in, fftw_complex *out, int sign, unsigned flags);
                                
           int n is the number of elements
           fftw_complex *in,*out are the pointers to the in array and the out array, if the same this defines an
                in-place transform, i.e. one that overwrites the input data             
           int sign indicates whether one is doing a forward or inverse transform
           unsigned flags refer to various ways in which the fft is computed. FFTW_ESTIMATE does not optmise and is best
                best for doing a transform only a few times. FFTW_MEASURE does initialisation tests to optimise ffts and is good
                if running many many ffts

        */
        forward=fftw_plan_dft_1d(N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
        backward=fftw_plan_dft_1d(N,out,in,FFTW_BACKWARD,FFTW_ESTIMATE);
        
        for(i=0;i<N;i++) {
                in[i][0] = (double)rand()/RAND_MAX;
                in[i][1] = (double)rand()/RAND_MAX;
        }
        printf("\n Initial Array\n");
        for(i=0;i<N;i++) {
                printf("%10f + %10f i\n", in[i][0],in[i][1]);
        }

        /* compute FFT which stores output in out since it is not an inplace transformation*/
        fftw_execute(forward);

        printf("\n FFT Array\n");
        for(i=0;i<N;i++) {
                /* make sure to print the out arrays */
                printf("%10f + %10f i\n", out[i][0],out[i][1]);
        }
        
        /* compute IFFT which stores outpit in in since it is not an inplace transformation */
        fftw_execute(backward);
        
        printf("\n Original Array\n");
        for(i=0;i<N;i++) {
                /* make sure to print the in arrays */
                printf("%10f + %10f i\n", in[i][0]/(double)N,in[i][1]/(double)N);
        }

        /* clean up */
        fftw_destroy_plan(forward);
        fftw_destroy_plan(backward);
        fftw_free(in);
        fftw_free(out);
        
        return 0;
}

