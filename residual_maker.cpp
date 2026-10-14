
#include "residual_maker.h"
// for each node i
LocalState make_local_state(
    const vector<vector<double>>& x,const vector<vector<double>>& x_old,
    int node,
    int n)
{
    LocalState state;

    state.im.resize(n);
    state.i.resize(n);
    state.ip.resize(n);
    state.i_old.resize(n);


    for (int i = 0; i < n; i++)
    {
        state.im[i] = x[i][node - 1];
        state.i[i]   = x[i][node];
        state.ip[i] = x[i][node + 1];
        state.i_old[i]   = x_old[i][node];
    }

    return state;
}

LocalState_left make_local_state_left(
    const vector<vector<double>>& x,const vector<vector<double>>& x_old,
    int node,
    int n)
{
    LocalState_left state_left;

    state_left.i.resize(n);
    state_left.ip.resize(n);
    state_left.i_old.resize(n);


    for (int i = 0; i < n; i++)
    {
        
        state_left.i[i]   = x[i][node];
        state_left.ip[i] = x[i][node+1];
        state_left.i_old[i]   = x_old[i][node];
    }

    return state_left;
}



LocalState_right make_local_state_right(
    const vector<vector<double>>& x,const vector<vector<double>>& x_old,
    int node,
    int n)
{
    LocalState_right state_right;

    state_right.im.resize(n);
    state_right.i.resize(n);
    state_right.i_old.resize(n);


    for (int i = 0; i < n; i++)
    {
        
        state_right.im[i]   = x[i][node-1];
        state_right.i[i] = x[i][node];
        state_right.i_old[i]   = x_old[i][node];
    }

    return state_right;
}
// bc residual making for left boundary
vector<double> residual_left(
    const LocalState_left& s,
    const vector<BCSpec>& bc_left,
    double del_x)
{
    vector<double> R_local(2, 0.0);

    for (int j = 0; j < 2; j++)
    {
        R_local[j] =
            bc_left[j].alpha * s.i[j]
            + bc_left[j].beta *
              (s.ip[j] - s.i[j]) / del_x
            - bc_left[j].theta;
    }

    return R_local;
}
// bc residual making for right boundary
vector<double> residual_right(
    const LocalState_right& s,
    const vector<BCSpec>& bc_right,
    double del_x)
{
    vector<double> R_local(2, 0.0);

    for (int j = 0; j < 2; j++)
    {
        R_local[j] =
            bc_right[j].alpha * s.i[j]
            + bc_right[j].beta *
              (s.i[j] - s.im[j]) / del_x
            - bc_right[j].theta;
    }

    return R_local;
}

// for local residual keep in mind it has to be hard coded for the physics of the problem like for 2 in rlocal is for 2 residuals
vector<double> residual(
    const LocalState& s,
    double del_x ,double dt)
{
    vector<double> R_local(2, 0.0);

    // Residual equation 1
    R_local[0] =
    (s.i[0] - s.i_old[0]) / dt
    - (s.im[0] - 2.0*s.i[0] + s.ip[0])/(del_x*del_x)
    - s.i[1];

    R_local[1] =
    (s.i[1] - s.i_old[1]) / dt
    - (s.im[1] - 2.0*s.i[1] + s.ip[1])/(del_x*del_x)
    - s.i[0];

    return R_local;
}

vector<vector<double>> res_func(
    int num,
    int N,
    int n,
    double del_x,
    double dt,
    const vector<vector<double>>& x,
    const vector<vector<double>>& x_old)
{
    vector<BCSpec> bc_left(n);
    vector<BCSpec> bc_right(n);
    bc_left[0]  = {1.0, 0.0, 0.0};
    bc_left[1]  = {1.0, 0.0, 0.0};
    bc_right[0] = {1.0, 0.0, 1.0};
    bc_right[1] = {1.0, 0.0, 1.0};
    //this has to be problem specific
    
    
    vector<vector<double>> R(
        N,
        vector<double>(num, 0.0)
    );

    
    // INTERIOR
    for (int i = 0; i < N ; i++)
    {
        if (i==0){
            
            LocalState_left s = make_local_state_left(x, x_old, i, n);
            vector<double> R_local = residual_left(s ,bc_left ,del_x);
            
            for (int j = 0; j < num; j++)
            {
                R[i][j] = R_local[j];
            }
            }
        
        else if ( i==N-1){
            
            LocalState_right s = make_local_state_right(x, x_old, i, n);
            vector<double> R_local = residual_right(s ,bc_right ,del_x);
            
             for (int j = 0; j < num; j++)
            {
                R[i][j] = R_local[j];
            }
        }
        else{
            LocalState s = make_local_state(x, x_old, i, n);

            vector<double> R_local = residual(s ,del_x ,dt);

            for (int j = 0; j < num; j++)
            {
                R[i][j] = R_local[j];
            }
        }
    }

    

    return R;
}