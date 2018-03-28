"""
Linear equation solver
* A is a square matrix and b is a vector.
* The solution x is returned in vector b
* 
* Reference:
* 	Golub and Van Loan, "Matrix Computations", 1999, Ch 3

void linsolve(REAL **A, REAL *b, int N){

  int i,j,k;
  REAL sumi;
  
  // Algorithm to find LU decomp - See page 99
  for(k=0;k<N-1;k++){
    for(i=k+1;i<N;i++){
      A[i][k] = A[i][k]/A[k][k];
      for(j=k+1;j<N;j++){
	A[i][j] = A[i][j] - A[k][j]*A[i][k];
      }
    }
  }

  // Solve L.y=b via forward substitution (Alg 3.1.1);
  b[0] = b[0];
  for(i=1;i<N;i++){
    sumi=0.0;
    for(j=0;j<i;j++){
      sumi = sumi + A[i][j]*b[j];
    }
    b[i]=b[i]-sumi;
  }
  
  //Solve U.x=y via backward substitution (Alg 3.1.2)
  
  b[N-1] = b[N-1]/A[N-1][N-1];
  for(i=N-2;i>=0;i--){
    sumi=0.0;
    for(j=i+1;j<N;j++){
      sumi = sumi + A[i][j]*b[j];
    }
    b[i] = (b[i] - sumi)/A[i][i];
  }  
} // End of linsolve
"""

import numpy as np 
from numba import jit

#void linsolve(REAL **A, REAL *b, int N){
@jit(nopython=False, error_model='numpy')
def linsolve(A, b, N):
    #// Algorithm to find LU decomp - See page 99
    for k in range(0,N-1):
        for i in range(k+1,N):
            if A[k,k]!=0.:
                A[i,k] = A[i,k]/A[k,k]

            for j in range(k+1,N):
                A[i,j] = A[i,j] - A[k,j]*A[i,k]
  
    # // Solve L.y=b via forward substitution (Alg 3.1.1);
    b[0] = b[0]
    for i in range(1,N):
        sumi = 0.
        for j in range(0,i):
            sumi = sumi + A[i,j]*b[j]
  
        b[i]=b[i]-sumi
 
    #//Solve U.x=y via backward substitution (Alg 4.1.2)
    if A[N-1,N-1] != 0.:
        b[N-1] = b[N-1]/A[N-1,N-1]

    for i in range(N-2,-1,-1):
        sumi = 0.
        for j in range(i+1,N):
            sumi = sumi + A[i,j]*b[j]

        if A[i,i] != 0.:
            b[i] = (b[i] - sumi)/A[i,i]

    return b


 
#def linsolve(A, b, N):
#    #// Algorithm to find LU decomp - See page 99
#    for k in range(0,N-1):
#        for i in range(k+1,N):
#            A[i,k] = A[i,k]/A[k,k]
#            for j in range(k+1,N):
#                A[i,j] = A[i,j] - A[k,j]*A[i,k]
#  
#    # // Solve L.y=b via forward substitution (Alg 3.1.1);
#    b[0] = b[0];
#    for i in range(N):
#        sumi = 0.
#        for j in range(0,i):
#            sumi = sumi + A[i,j]*b[j]
#  
#        b[i]=b[i]-sumi
# 
#    #//Solve U.x=y via backward substitution (Alg 3.1.2)
#
#    b[N-1] = b[N-1]/A[N-1,N-1]
#    for i in range(N-2,0):
#        sumi = 0.
#        for j in range(i+1,N):
#            sumi = sumi + A[i,j]*b[j]
#        b[i] = (b[i] - sumi)/A[i,i]
#
#    return b

  
