#include"linealgbr.h"
#include"BTDMA.h"

// BLOCK TDMA solver
vector<vector<double>> BTDMA(int N , int n , const vector<vector<vector<double>>> &A ,const vector<vector<vector<double>>> &B , const vector<vector<vector<double>>> &C ,const vector<vector<double>> &D){
   vector<vector<vector<double>>> C_h(N, vector<vector<double>>(n, vector<double>(n, 0.0f))); 
    vector<vector<double>> D_h(N, vector<double>(n, 0.0f));
    vector<vector<double>> del_U(N, vector<double>(n, 0.0f));
    vector<vector<double>> M(n, vector<double>(n,0.0f));
    vector<vector<double>> B_STAR(n, vector<double>(n,0.0f));
    vector<vector<double>> B_inv(n, vector<double>(n,0.0f));;
    vector<double> D_STAR(n, 0.0f);

    // FORWARD ELIMINATION
    C_h[0]= matrix_multiply (matrix_inversion(B[0]),C[0]);
    D_h[0]= matrix_vector_multiply (matrix_inversion(B[0]),D[0]);
  
    
    for (int i=1 ;i<N-1 ;i++){
        M = matrix_multiply(A[i] ,C_h[i-1]);
        B_STAR = matrix_subtraction(B[i],M);
        B_inv = matrix_inversion(B_STAR);
        D_STAR = vector_subtraction(D[i],matrix_vector_multiply(A[i],D_h[i-1]));
        C_h[i] = matrix_multiply(B_inv,C[i]);
        D_h[i] = matrix_vector_multiply(B_inv,D_STAR);


    }


    M = matrix_multiply(A[N-1],C_h[N-2]);
    B_STAR = matrix_subtraction(B[N-1],M);
    D_STAR = vector_subtraction(D[N-1],matrix_vector_multiply(A[N-1],D_h[N-2]));
    

    // BACK SUBSTITUTION
    B_inv = matrix_inversion(B_STAR);
    del_U[N-1] = matrix_vector_multiply(B_inv,D_STAR);

    for(int i=N-2 ; i>=0 ;i--){
        del_U[i] = vector_subtraction(D_h[i],matrix_vector_multiply(C_h[i],del_U[i+1]));
    }

    return del_U;

}