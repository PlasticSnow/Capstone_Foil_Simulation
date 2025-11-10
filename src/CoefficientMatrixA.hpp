#pragma once
#include <vector>
#include <iostream>



class MatrixA {

    public:;
        // 2D Representation of 2D Matrix
        std::vector<double> aDiag, aPlusI, aPlusJ;
        int Nx_, Ny_;

        MatrixA();
        

        MatrixA(int Nx, int Ny, double dt, double rho, double dx, std::vector<bool> fluidState);

        // Index helper function
        inline int idx(int i, int j) const { return i + Nx_*j; }

        // Comput y = A * x
        void applyA(const std::vector<double>& givenPressure, std::vector<double>& resultPressure, const std::vector<bool>& fluidState);


        void createAMatrices(int Nx, int Ny, double dt, double rho, double dx,const std::vector<bool>& fluidState);
};


