#pragma once 
#include "FluidSim.hpp"
#include "CoefficientMatrixA.hpp"
#include <iostream>
#include <string>
#include <format>

void testMatrixASetup(int Nx, int Ny, double dt, double rho, double dx);

void printDiag(std::vector<double> aDiag, int Nx, int Ny);

std::string doub_to_string(std::string val);