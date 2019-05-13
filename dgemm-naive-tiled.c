/* 
    Please include compiler name below (you may also include any other modules you would like to be loaded)

COMPILER= gnu

    Please include All compiler flags and libraries as you want them run. You can simply copy this over from the Makefile's first few lines
 
CC = cc
OPT = -O3
CFLAGS = -Wall -std=gnu99 $(OPT)
MKLROOT = /opt/intel/composer_xe_2013.1.117/mkl
LDLIBS = -lrt -Wl,--start-group $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_sequential.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -lpthread -lm

*/

//Average percentage of Peak = 26.2667
//Grade = 39.4001
#include <stdint.h>

const char* dgemm_desc = "Naive, three-loop dgemm.";

/* This routine performs a dgemm operation
 *  C := C + A * B
 * where A, B, and C are lda-by-lda matrices stored in column-major format.
 * On exit, A and B maintain their input values. */    
void square_dgemm (int n, double* A, double* B, double* C)
{
//transpose A
  double *tA;
  tA = (double *) malloc(n*n * sizeof(double));
  for (int i=0; i<n; ++i)
   for (int j=0; j<n; ++j) {
    tA[j*n+i] = A[i*n+j];
   }
//Pad Matrices
 int step = 6;
 int p = n%step;
 int s = n + (step-p);
 double *pA;
 pA = (double *) malloc(s*s * sizeof(double));
 double *pB;
 pB = (double *) malloc(s*s * sizeof(double));
 double *pC;
 pC = (double *) malloc(s*s * sizeof(double));
 for (int i=0; i<s; ++i) {
		for (int j=0; j<s; ++j) {
			if(j>=n || i >=n) {
				pA[i*s+j] = 0;
				pB[i*s+j] = 0;
				pC[i*s+j] = 0;
				continue;
			}
				pA[i*s+j] = tA[i*n+j];
				pB[i*s+j] = B[i*n+j];
				pC[i+j*s] = C[i+j*n];  
		}
 }	
 free(tA);

	/* For each row i of A */
	for (int i = 0; i < s; i+=6 )
    { 
	    /* For each column j of B */
	    for (int j = 0; j < s; j+=3  ) 
	    {
	      /* Compute C(i,j) */
	      double cij = pC[i+j*s];
		  double cij1 = pC[i+(j+1)*s];
          double cij2 = pC[i+(j+2)*s];

		  double cijI = pC[(i+1)+j*s];
		  double cij1I = pC[(i+1)+(j+1)*s];
          double cij2I = pC[(i+1)+(j+2)*s];

          double cijI2 = pC[(i+2)+j*s];
		  double cij1I2 = pC[(i+2)+(j+1)*s];
          double cij2I2 = pC[(i+2)+(j+2)*s];

          double cijI3 = pC[(i+3)+j*s];
		  double cij1I3 = pC[(i+3)+(j+1)*s];
          double cij2I3 = pC[(i+3)+(j+2)*s];

          double cijI4 = pC[(i+4)+j*s];
		  double cij1I4 = pC[(i+4)+(j+1)*s];
          double cij2I4 = pC[(i+4)+(j+2)*s];

          double cijI5 = pC[(i+5)+j*s];
		  double cij1I5 = pC[(i+5)+(j+1)*s];
          double cij2I5 = pC[(i+5)+(j+2)*s];

	      for( int k = 0; k < s; k++ )
	      {
			cij += pA[i*s+k] * pB[k+j*s];
			cij1 += pA[i*s+k] * pB[k+(j+1)*s];
            cij2 += pA[i*s+k] * pB[k+(j+2)*s];

			cijI += pA[(i+1)*s+k] * pB[k+j*s];
			cij1I += pA[(i+1)*s+k] * pB[k+(j+1)*s];
            cij2I += pA[(i+1)*s+k] * pB[k+(j+2)*s];

            cijI2 += pA[(i+2)*s+k] * pB[k+j*s];
			cij1I2 += pA[(i+2)*s+k] * pB[k+(j+1)*s];
            cij2I2 += pA[(i+2)*s+k] * pB[k+(j+2)*s];

            cijI3 += pA[(i+3)*s+k] * pB[k+j*s];
			cij1I3 += pA[(i+3)*s+k] * pB[k+(j+1)*s];
            cij2I3 += pA[(i+3)*s+k] * pB[k+(j+2)*s];

            cijI4 += pA[(i+4)*s+k] * pB[k+j*s];
			cij1I4 += pA[(i+4)*s+k] * pB[k+(j+1)*s];
            cij2I4 += pA[(i+4)*s+k] * pB[k+(j+2)*s];

            cijI5 += pA[(i+5)*s+k] * pB[k+j*s];
			cij1I5 += pA[(i+5)*s+k] * pB[k+(j+1)*s];
            cij2I5 += pA[(i+5)*s+k] * pB[k+(j+2)*s];
	      }
	      pC[i+j*s] = cij;
		  pC[i+(j+1)*s] = cij1;
          pC[i+(j+2)*s] = cij2;

		  pC[(i+1)+j*s] = cijI;
		  pC[(i+1)+(j+1)*s] = cij1I;
          pC[(i+1)+(j+2)*s] = cij2I;

          pC[(i+2)+j*s] = cijI2;
		  pC[(i+2)+(j+1)*s] = cij1I2;
          pC[(i+2)+(j+2)*s] = cij2I2;

          pC[(i+3)+j*s] = cijI3;
		  pC[(i+3)+(j+1)*s] = cij1I3;
          pC[(i+3)+(j+2)*s] = cij2I3;

          pC[(i+4)+j*s] = cijI4;
		  pC[(i+4)+(j+1)*s] = cij1I4;
          pC[(i+4)+(j+2)*s] = cij2I4;

          pC[(i+5)+j*s] = cijI5;
		  pC[(i+5)+(j+1)*s] = cij1I5;
          pC[(i+5)+(j+2)*s] = cij2I5;
	    }
    }

	  //Unpad C
	  for (int i=0; i<s; ++i){
			   for (int j=0; j<s; ++j) {
				if(j>=n || i >=n) {
					continue;
				} 
				C[i+j*n] = pC[i+j*s]; 
			}
	  }
	  free(pA);
	  free(pB);
	  free(pC);
}
