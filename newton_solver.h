#pragma once

#include <vector>
#include "jacobian_maker.h"

using std::vector;

vector<vector<double>> newton_raphson(
    const Jacobian& J,
    const vector<vector<double>>& R);