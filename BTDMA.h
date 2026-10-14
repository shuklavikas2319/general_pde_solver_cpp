#pragma once
#include <vector>
using std::vector;
vector<vector<double>> BTDMA(
    int N, int n,
    const vector<vector<vector<double>>>& A,
    const vector<vector<vector<double>>>& B,
    const vector<vector<vector<double>>>& C,
    const vector<vector<double>>& D);