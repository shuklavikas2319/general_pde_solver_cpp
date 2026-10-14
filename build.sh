#!/bin/bash
g++ -std=c++17 -Wall -o solver main.cpp BTDMA.cpp jacobian_maker.cpp linealgbr.cpp newton_solver.cpp residual_maker.cpp