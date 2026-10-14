#pragma once
#include <vector>
using std::vector;

struct LocalState
{
    vector<double> im, i, ip, i_old;
};
struct LocalState_left
{
    vector<double> i, ip, i_old;
};
struct LocalState_right
{
    vector<double> im, i, i_old;
};
struct BCSpec
{
    double alpha, beta, theta;
};

LocalState make_local_state(const vector<vector<double>>& x, const vector<vector<double>>& x_old, int node, int n);
LocalState_left make_local_state_left(const vector<vector<double>>& x, const vector<vector<double>>& x_old, int node, int n);
LocalState_right make_local_state_right(const vector<vector<double>>& x, const vector<vector<double>>& x_old, int node, int n);

vector<double> residual(const LocalState& s, double del_x ,double dt);
vector<double> residual_left(const LocalState_left& s, const vector<BCSpec>& bc_left, double del_x);
vector<double> residual_right(const LocalState_right& s, const vector<BCSpec>& bc_right, double del_x);

vector<vector<double>> res_func(int num, int N, int n, double del_x , double dt,const vector<vector<double>>& x,const vector<vector<double>>& x_old);