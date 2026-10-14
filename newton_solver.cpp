#include <iostream>
#include "jacobian_maker.h"
#include "BTDMA.h"
#include "linealgbr.h"
using std::cout;
using std::endl;

vector<vector<double>> newton_raphson(
    const Jacobian& J,
    const vector<vector<double>>& R)
{
    int N = R.size();
    int n = R[0].size();

    if (J.A.size() != N ||
        J.B.size() != N ||
        J.C.size() != N)
    {
        cout << "Error: Jacobian and residual dimensions do not match."
             << endl;
        return {};
    }

    vector<vector<double>> D(
        N, vector<double>(n, 0.0)
    );

    // RHS = -R
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < n; j++)
        {
            D[i][j] = -R[i][j];
        }
    }

    vector<vector<double>> del_U =
        BTDMA(N, n, J.A, J.B, J.C, D);

    return matrix_TRANSPOSE(del_U);
}