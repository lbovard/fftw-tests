#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>

// get A_{ij} w/ ind_sqices i,j=1,...,N
#define ind_sq(i,j) (i-1)*N+(j-1)
#define ind_fft(i,j) (i-1)*(N/2+1)+(j-1)
#define nl printf("\n");

void print_matrix(double *matrix, int mtype, int N) {
       int i,j;
       switch(mtype) {
        case 0:
               for(i=1;i<=N;i++) {
                       for(j=1;j<=N;j++) {
                               printf("%lf\t ",matrix[ind_sq(i,j)]);
                       }
                       printf("\n"); 
               }        
               break;
        case 1:
               for(i=1;i<=N;i++) {
                       for(j=1;j<=(N/2+1);j++) {
                               printf("%lf ",matrix[ind_fft(i,j)]);
                       }
                       printf("\n"); 
               }        
               break;
        default: 
                printf("wrong matrix type\n");
                break;
        }
}
/* normalise because FFTW does not*/
void normalise(double *matrix, int N) {
        int i,j;
        for(i=1;i<=N;i++) {
                for(j=1;j<=N;j++) {
                        matrix[ind_sq(i,j)]/=(double)(N*N);
                }
        }
}
void print_fft(fftw_complex *matrix, int N) {
        int i,j;
        for(i=1;i<=N;i++) {
                for(j=1;j<=(N/2+1);j++) {
                        if(matrix[ind_fft(i,j)][1]<0)
                        printf("%lf%lfi\t", matrix[ind_fft(i,j)][0],matrix[ind_fft(i,j)][1]);
                        else
                        printf("%lf+%lfi\t", matrix[ind_fft(i,j)][0],matrix[ind_fft(i,j)][1]);
                }
                printf("\n");
        }
}
void print_fft_1d(fftw_complex *matrix, int N) {
        int i;
        for(i=0;i<N*(N/2+1);i++) {
		if(matrix[i][1]<0)
		printf("%lf%lfi\t", matrix[i][0],matrix[i][1]);
		else
		printf("%lf+%lfi\t", matrix[i][0],matrix[i][1]);
                printf("\n");
        }
}

/* compute the spectral derivative */
void derivative(fftw_complex *matrix,double *wn, int N) {
        int i,j;
        double tmp;
        for(i=1;i<=N;i++) {
                for(j=1;j<=(N/2+1);j++) {       
                        tmp = matrix[ind_fft(i,j)][0];
                        matrix[ind_fft(i,j)][0]=-wn[ind_fft(i,j)]*matrix[ind_fft(i,j)][1];
                        matrix[ind_fft(i,j)][1]=wn[ind_fft(i,j)]*tmp; 
                        }
        }

}
int main() {
        int N=4,i,j,Nr=N/2+1;
        double *v;
        double *k_x,*k_y, tmp;
        fftw_complex *out;
        fftw_plan forward, inverse;
        /* memory allocation */ 
       // v=fftw_malloc(N*N*sizeof(double));   
	v=fftw_alloc_real(N*N);
        k_x=fftw_malloc(N*(N/2+1)*sizeof(double));
        k_y=fftw_malloc(N*(N/2+1)*sizeof(double));
        out=(fftw_complex*)fftw_malloc(N*(N/2+1)*sizeof(fftw_complex));
        /* create plans */ 
        forward=fftw_plan_dft_r2c_2d(N,N,v,out,FFTW_ESTIMATE);
        inverse=fftw_plan_dft_c2r_2d(N,N,out,v,FFTW_ESTIMATE); 
        /* set up the k_x,k_y vectors */
        for(i=0;i<N;i++) {
                for(j=0;j<(N/2)+1;j++) {
                        k_x[ind_fft(i+1,j+1)]=(double)(j);        
                        if(i<N/2) {
                                k_y[ind_fft(i+1,j+1)]=(double)(i);
                        }
                        else if(i==N/2) {
                                k_y[ind_fft(i+1,j+1)]=(double)(0.0);
                        }
                        else  {
                                k_y[ind_fft(i+1,j+1)]=(double)(i-N);
                        }
                }
        }
        /* initialise matrix */
        
        for(i=1;i<=N;i++) {
                for(j=1;j<=N;j++) {
                     //   v[ind_sq(i,j)]=(double)(i*i-j);
                        v[ind_sq(i,j)]=1.0/(double)(i+j-1);
                }
        }
        for(i=1;i<N;i++) {
                v[ind_sq(i,i+1)]=(double)(i+1);
       } 
//        print_matrix(v,0,N); 

        fftw_execute(forward); 
//        print_fft(out,N);nl
	print_fft_1d(out,N); nl
        derivative(out,k_x,N); 
 //       print_fft(out,N);
        fftw_execute(inverse);

        nl
        normalise(v,N);
//	print_matrix(v,0,N);
        /* clean up */
        fftw_destroy_plan(forward);
        fftw_destroy_plan(inverse);
        fftw_free(v);
        fftw_free(out); 
        return 0;
}
