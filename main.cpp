#include<iostream>
#include <cmath>
#include <algorithm>
#include"BTDMA.h"
#include"jacobian_maker.h"
#include"linealgbr.h"
#include"residual_maker.h"
#include"newton_solver.h"
using namespace std;
double residual_norm(const vector<vector<double>>& R)
{
    double max_R = 0.0f;

    for (int i = 0; i < R.size(); i++)
    {
        for (int j = 0; j < R[i].size(); j++)
        {
            max_R = max(max_R, fabs(R[i][j]));
        }
    }

    return max_R;
}
int main (){
    int n , N ,num;
    double time_target ;
    cout<<"give the value of number of variables"<<endl;
    cin>>n;
    cout<<"give the value of number of nodes"<<endl;
    cin>>N;
    cout<<"give the value of number of residual equations"<<endl;
    cin>>num;
    cout<<"give the value of target time time_target"<<endl;
    cin>>time_target;

    double L=1 ,del_x=L/(N-1) , dt = 0.01 ,current_time=0;
    vector<vector<double>>x(n, vector<double>(N,0));// store value of the variables at each node with initialization of 0
    vector<vector<double>>x_old(n, vector<double>(N,0));// store value of the variables at each node with initialization of 0
  
  

    vector<vector<double>> R(N, vector<double>(num, 0.0f));

    vector<vector<double>> del_U;
    while ( current_time< time_target){
        double tolerance = 1e-6;
        double R_max ;
        R = res_func(num,N,n,del_x ,dt,x,x_old);
        R_max = residual_norm(R);
        int iter = 0;
        int max_iter = 100;

        while (R_max > tolerance && iter < max_iter)
        {       vector<BCSpec> bc_left(n);
            vector<BCSpec> bc_right(n);
            bc_left[0]  = {1.0, 0.0, 0.0};
            bc_left[1]  = {1.0, 0.0, 0.0};
            bc_right[0] = {1.0, 0.0, 1.0};
            bc_right[1] = {1.0, 0.0, 1.0};
            Jacobian J = jacobian_maker(n, N , num, del_x, dt, x, x_old, bc_left , bc_right);
            del_U = newton_raphson(J , R);// suggested update from newton raphson BTDMA

            x = matrix_addition(x, del_U); //updating current variables

            R = res_func(num,N,n,del_x,dt,x,x_old);

            R_max = residual_norm(R);//CALCULATING MAXIMUM IN ALL RESIDUALS FROM ALL EQAUTION

            iter++;

            cout << "Iteration = " << iter
                << "   Residual = " << R_max << endl;
            
        }

        if (R_max <= tolerance){
            cout << "Newton-Raphson converged!" << endl;
        }
        else{
            cout << "Newton-Raphson did not converge." << endl;
        }
        current_time+=dt;
        x_old=x;
    }
    for(int i=0; i<n ;i++){
        for(int j=0 ;j<N ;j++){
            cout<<"value of updated x at row"<<i<<"&column"<<j<<"is="<<x[i][j]<<endl;
        }
    }

    return 0;
}