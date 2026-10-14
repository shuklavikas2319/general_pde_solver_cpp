#pragma once
#include <vector>
#include "residual_maker.h"
using std::vector;
struct Jacobian
{
    vector<vector<vector<double>>> A;
    vector<vector<vector<double>>> B;
    vector<vector<vector<double>>> C;
};
Jacobian jacobian_maker(int n, int  N , int num,  double del_x,double dt, const vector<vector<double>> &x,const vector<vector<double>> &x_old, const vector<BCSpec> & bc_left , const vector<BCSpec> & bc_right);
