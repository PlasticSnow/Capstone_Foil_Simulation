#include "CoefficientMatrixA.hpp"


MatrixA::MatrixA(){
    Nx_ = 1;
    Ny_ = 1;
    aDiag = std::vector<double>(Nx_*Ny_, 0.0);
    aPlusI = std::vector<double>(Nx_*Ny_, 0.0);
    aPlusJ = std::vector<double>(Nx_*Ny_, 0.0);
}

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


void MatrixA::applyA(const std::vector<double>& givenPs, std::vector<double>& resultPs, const std::vector<bool>& fluidState){
    for (int j = 0; j < Ny_; j++){
        for (int i = 0; i < Nx_; i++){
            
            int index = idx(i, j);

            // Check to see if solid or empty, if so skip
            if (!fluidState[index]){
                resultPs[index] = 0.0;
                continue;
            }


            double sum = 0.0;

            // Diagonal term: multiply by itself
            sum += aDiag[index] * givenPs[index];

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
            resultPs[index] = sum;

        }
    }
}



void MatrixA::createAMatrices(int Nx, int Ny, double dt, double rho, double dx, const std::vector<bool>& fluidState) {
    double scale = 1.0 / (dx * dx);  // Poisson scale

    std::fill(aDiag.begin(), aDiag.end(), 0.0);
    std::fill(aPlusI.begin(), aPlusI.end(), 0.0);
    std::fill(aPlusJ.begin(), aPlusJ.end(), 0.0);

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idxC = idx(i, j);
            if (!fluidState[idxC]) continue;

            // X-direction neighbor (right)
            if (i < Nx - 1) {
                int idxR = idx(i + 1, j);
                if (fluidState[idxR]) {
                    aDiag[idxC] += scale;
                    aDiag[idxR] += scale;
                    aPlusI[idxC] = -scale;
                } else {
                    // solid or air → treat as Dirichlet wall
                    aDiag[idxC] += scale;
                }
            }

            // Y-direction neighbor (top)
            if (j < Ny - 1) {
                int idxT = idx(i, j + 1);
                if (fluidState[idxT]) {
                    aDiag[idxC] += scale;
                    aDiag[idxT] += scale;
                    aPlusJ[idxC] = -scale;
                } else {
                    // solid or air → treat as Dirichlet wall
                    aDiag[idxC] += scale;
                }
            }
        }
    }
}
