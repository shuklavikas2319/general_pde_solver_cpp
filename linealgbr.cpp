#include "linealgbr.h"

// matrix addition
vector<vector<double>> matrix_addition( const vector<vector<double>> & A , const vector<vector<double>> & B )
{
   int rows= A.size();
   int column= A[0].size();
   vector<vector<double>>C(rows , vector<double>(column,0));
   for (int x=0; x<rows ; x++){
    for( int y=0 ; y<column ; y++){
        C[x][y]=A[x][y]+B[x][y];
    }
   }
   return C;
}
// matrix subtraction
vector<vector<double>> matrix_subtraction( const vector<vector<double>> & A , const vector<vector<double>> & B ){
   int rows= A.size();
   int column= A[0].size();
   vector<vector<double>>C(rows , vector<double>(column,0.0f));
   for (int x=0; x<rows ; x++){
    for( int y=0 ; y<column ; y++){
        C[x][y]=A[x][y]-B[x][y];
    }
   }
   return C;
}
// 2-D TRANSPOSE 
vector<vector<double>> matrix_TRANSPOSE( const vector<vector<double>> & A  ){
   int rows = A.size();
   int column = A[0].size();
   vector<vector<double>>C(column , vector<double>(rows,0.0f));
   for(int i=0; i<rows ;i++){
    for(int j=0 ;j<column ;j++){
        C[j][i]=A[i][j];
    }
   }
   return C;
}
// matrix to vector multiplication
vector<double> matrix_vector_multiply(
    const vector<vector<double>>& A,
    const vector<double>& x)
{
    int rows = A.size();
    int cols = A[0].size();

    vector<double> result(rows, 0.0f);

    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            result[i] += A[i][j] * x[j];
        }
    }

    return result;
}
// vector subtraction
vector<double> vector_subtraction(
    const vector<double>& A,
    const vector<double>& B)
{
    vector<double> C(A.size(), 0.0f);

    for (int i = 0; i < A.size(); i++)
        C[i] = A[i] - B[i];

    return C;
}
// matrix to matrix multiply
vector<vector<double>> matrix_multiply(
    const vector<vector<double>>& A,
    const vector<vector<double>>& B)
{
    int rows_A = A.size();
    int cols_A = A[0].size();
    int cols_B = B[0].size();

    vector<vector<double>> C(
        rows_A,
        vector<double>(cols_B, 0.0f)
    );

    for (int i = 0; i < rows_A; i++)
    {
        for (int j = 0; j < cols_B; j++)
        {
            for (int k = 0; k < cols_A; k++)
            {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }

    return C;
}
// matrix inversion
vector<vector<double>> matrix_inversion(
    const vector<vector<double>>& A)
{
    int N = A.size();

    vector<vector<double>> X(N, vector<double>(N, 0.0f));

    // Solve AX = I, one column at a time
    for (int p = 0; p < N; p++)
    {
        vector<vector<double>> C = A;

        // Current column of identity matrix
        vector<double> d(N, 0.0f);
        d[p] = 1.0f;

        vector<double> x(N, 0.0f);
        
        // Forward elimination

        for (int k = 0; k < N - 1; k++)
        {
            for (int i = k + 1; i < N; i++)
            {
                double factor = C[i][k] / C[k][k];

                for (int j = k; j < N; j++)
                {
                    C[i][j] -= factor * C[k][j];
                }

                d[i] -= factor * d[k];
            }
        }
        // Back substitution

        for (int i = N - 1; i >= 0; i--)
        {
            double sum = 0.0f;

            for (int j = i + 1; j < N; j++)
            {
                sum += C[i][j] * x[j];
            }

            x[i] = (d[i] - sum) / C[i][i];
        }

        // Store solution as column p
        for (int i = 0; i < N; i++)
        {
            X[i][p] = x[i];
        }
    }

    return X;
}
