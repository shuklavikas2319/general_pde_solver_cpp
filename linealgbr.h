#pragma once
#include<vector>
using std::vector;

vector<vector<double>> matrix_addition( const vector<vector<double>> & A , const vector<vector<double>> & B );
vector<vector<double>> matrix_subtraction( const vector<vector<double>> & A , const vector<vector<double>> & B );
vector<vector<double>> matrix_TRANSPOSE( const vector<vector<double>> & A  );
vector<double> matrix_vector_multiply(const vector<vector<double>>& A, const vector<double>& x);
vector<double> vector_subtraction(const vector<double>& A, const vector<double>& B);
vector<vector<double>> matrix_multiply(const vector<vector<double>>& A, const vector<vector<double>>& B);
vector<vector<double>> matrix_inversion (const vector<vector<double>>& A);