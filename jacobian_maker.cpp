#include <iostream>
#include "jacobian_maker.h"
#include "residual_maker.h"
#include "BTDMA.h"
#include "linealgbr.h"
using std::cout;
using std::endl;

Jacobian jacobian_maker(int n, int  N , int num,  double del_x, double dt, const vector<vector<double>> &x,const vector<vector<double>> &x_old, const vector<BCSpec> & bc_left , const vector<BCSpec> & bc_right){
    if (num != n)
    {
        cout << "Error: Newton-BTDMA requires num == n." << endl;
        return {};
    }
    double epsilon =1.e-9;
    Jacobian J;

    J.A.resize(N, vector<vector<double>>(n, vector<double>(n, 0.0)));
    J.B.resize(N, vector<vector<double>>(n, vector<double>(n, 0.0)));
    J.C.resize(N, vector<vector<double>>(n, vector<double>(n, 0.0)));
    
    vector<double> R_local(num);
    vector<double> R_perturbed(num);
    vector<double> R_localleft(num);
    vector<double> R_localright(num);
    vector<double> R_perturbedleft(num);
    vector<double> R_perturbedright(num);
    // creating perturbed residual for numerically calculating jacobian
    
    
  
// creating  Jacobians

for (int i = 0; i < N; i++)
{
    // =========================================================
    // LEFT BOUNDARY
    // =========================================================
    if (i == 0)
    {
        LocalState_left s =
            make_local_state_left(x, x_old, i, n);

        R_localleft =
            residual_left(s, bc_left, del_x);

        vector<vector<vector<double>>> R_pl(
            2,
            vector<vector<double>>(n,
                vector<double>(num, 0.0))
        );

        for (int k = 0; k < 2 * n; k++)
        {
            LocalState_left s_perturbed = s;

            int neighbour;
            int variable;

            if (k < n)
            {
                // central node i
                neighbour = 0;
                variable = k;

                s_perturbed.i[variable] += epsilon;
            }
            else
            {
                // right node i+1
                neighbour = 1;
                variable = k - n;

                s_perturbed.ip[variable] += epsilon;
            }

            R_perturbedleft =
                residual_left(
                    s_perturbed,
                    bc_left,
                    del_x
                );

            for (int r = 0; r < num; r++)
            {
                R_pl[neighbour][variable][r] =
                    R_perturbedleft[r];
            }
        }

        // A[0] does not exist
        for (int j = 0; j < num; j++)
        {
            for (int k = 0; k < n; k++)
            {
                J.B[0][j][k] =
                    (R_pl[0][k][j] - R_localleft[j])
                    / epsilon;

                J.C[0][j][k] =
                    (R_pl[1][k][j] - R_localleft[j])
                    / epsilon;
            }
        }
    }

    // =========================================================
    // RIGHT BOUNDARY
    // =========================================================
    else if (i == N - 1)
    {
        LocalState_right s =
            make_local_state_right(x, x_old, i, n);

        R_localright =
            residual_right(s, bc_right, del_x);

        vector<vector<vector<double>>> R_pr(
            2,
            vector<vector<double>>(n,
                vector<double>(num, 0.0))
        );

        for (int k = 0; k < 2 * n; k++)
        {
            LocalState_right s_perturbed = s;

            int neighbour;
            int variable;

            if (k < n)
            {
                // left node i-1
                neighbour = 0;
                variable = k;

                s_perturbed.im[variable] += epsilon;
            }
            else
            {
                // central node i
                neighbour = 1;
                variable = k - n;

                s_perturbed.i[variable] += epsilon;
            }

            R_perturbedright =
                residual_right(
                    s_perturbed,
                    bc_right,
                    del_x
                );

            for (int r = 0; r < num; r++)
            {
                R_pr[neighbour][variable][r] =
                    R_perturbedright[r];
            }
        }

        // C[N-1] does not exist
        for (int j = 0; j < num; j++)
        {
            for (int k = 0; k < n; k++)
            {
                J.A[N - 1][j][k] =
                    (R_pr[0][k][j] - R_localright[j])
                    / epsilon;

                J.B[N - 1][j][k] =
                    (R_pr[1][k][j] - R_localright[j])
                    / epsilon;
            }
        }
    }

    // =========================================================
    // INTERIOR
    // =========================================================
    else
    {
        LocalState s =
            make_local_state(x, x_old, i, n);

        R_local =
            residual(s, del_x ,dt);

        vector<vector<vector<double>>> R_p(
            3,
            vector<vector<double>>(n,
                vector<double>(num, 0.0))
        );

        for (int k = 0; k < 3 * n; k++)
        {
            LocalState s_perturbed = s;

            int neighbour;
            int variable;

            if (k < n)
            {
                neighbour = 0;
                variable = k;

                s_perturbed.im[variable] += epsilon;
            }
            else if (k < 2 * n)
            {
                neighbour = 1;
                variable = k - n;

                s_perturbed.i[variable] += epsilon;
            }
            else
            {
                neighbour = 2;
                variable = k - 2 * n;

                s_perturbed.ip[variable] += epsilon;
            }

            R_perturbed =
                residual(s_perturbed, del_x ,dt);

            for (int r = 0; r < num; r++)
            {
                R_p[neighbour][variable][r] =
                    R_perturbed[r];
            }
        }

        for (int j = 0; j < num; j++)
        {
            for (int k = 0; k < n; k++)
            {
                J.A[i][j][k] =
                    (R_p[0][k][j] - R_local[j])
                    / epsilon;

                J.B[i][j][k] =
                    (R_p[1][k][j] - R_local[j])
                    / epsilon;

                J.C[i][j][k] =
                    (R_p[2][k][j] - R_local[j])
                    / epsilon;
            }
        }
    }
}

return J;
}

