#include "CoefficientMatrixA.hpp"


MatrixA::MatrixA(int Nx, int Ny, double dt, double rho, double dx, std::vector<bool> fluidState){
    Nx_ = Nx;
    Ny_ = Ny;
    
    std::cout << "Making MatrixA" << std::endl;

    aDiag = std::vector<double>(Nx*Ny, 0.0);
    aPlusI = std::vector<double>(Nx*Ny, 0.0);
    aPlusJ = std::vector<double>(Nx*Ny, 0.0);

    std::cout << "Initialized A vectors" << std::endl;

    createAMatrices(Nx, Ny, dt, rho, dx, fluidState);

    std::cout << "Created MatrixA " << std::endl;


}


void MatrixA::applyA(std::vector<double>& givenPs, std::vector<double>& resultPs, std::vector<bool>& fluidState){
    for (int j = 0; j < Ny_; j++){
        for (int i = 0; i < Nx_; i++){
            
            // Check to see if solid or empty, if so skip
            if (!fluidState[idx(i,j)])
                continue;


            double sum = 0.0;

            // Diagonal term: multiply by itself
            sum += aDiag[idx(i,j)] * givenPs[idx(i,j)];

            // Right neighbor
            if (i < Nx_ - 1 && fluidState[idx(i+1,j)]){
                sum += aPlusI[idx(i,j)] * givenPs[idx(i+1,j)];
            }

            // Left neighbor 
            if (i > 0 && fluidState[idx(i-1,j)]) {
                sum += aPlusI[idx(i-1,j)] * givenPs[idx(i-1,j)];
            }  

            // Top neighbor 
            if (j < Ny_ - 1 && fluidState[idx(i,j+1)]) {
                sum += aPlusJ[idx(i,j)] * givenPs[idx(i,j+1)];   
            }

            // Bottom neighbor 
            if (j > 0 && fluidState[idx(i,j-1)]) {
                sum += aPlusJ[idx(i,j-1)] * givenPs[idx(i,j-1)]; 
            }

            // Store result
            resultPs[idx(i,j)] = sum;

        }
    }
}


void MatrixA::createAMatrices(int Nx, int Ny, double dt, double rho, double dx, std::vector<bool> fluidState){
    double scale = dt / (rho*dx*dx);

    std::cout << "About to create A Matrices" << std::endl;

    for (int j = 0; j < Ny; j++){
        for (int i = 0; i < Nx; i++){

            
            if (i < Nx - 1){

                // If current cell is FLUID and right neighbor is FLUID
                if (fluidState[idx(i,j)] && fluidState[idx(i+1,j)]){
                    aDiag[idx(i,j)] += scale;
                    aDiag[idx(i+1,j)] += scale;
                    aPlusI[idx(i,j)] = -scale;
                } else if (fluidState[idx(i,j)] && !fluidState[idx(i+1,j)]) {
                    aDiag[idx(i,j)] += scale;
                }

            }
            
            if (j < Ny - 1){

                // If current cell is FLUID and bottom neighbor is FLUID
                if (fluidState[idx(i,j)] && fluidState[idx(i,j+1)]){
                    aDiag[idx(i,j)] += scale;
                    aDiag[idx(i,j+1)] += scale;
                    aPlusJ[idx(i,j)] = -scale;
                } else if (fluidState[idx(i,j)] && !fluidState[idx(i,j+1)]){
                    aDiag[idx(i,j)] += scale;
                }

            }

        }
    }
}