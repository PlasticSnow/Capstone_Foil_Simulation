#include "Tests.hpp"

void testMatrixASetup(int Nx, int Ny, double dt, double rho, double dx){

    std::vector<bool> fluidState = std::vector<bool>(Nx*Ny, true);

    MatrixA matrixA = MatrixA(Nx, Ny, dt, rho, dx, fluidState);

    std::string matrixStringDiag = ""; 
    std::string matrixStringPlusI = "";
    std::string matrixStringPlusJ = "";   

    std::cout << "About to print" << std::endl;

    for (int j = 0; j < Ny; j++){
        matrixStringDiag.append("(");
        matrixStringPlusI.append("(");
        matrixStringPlusJ.append("(");

        double val;
        for (int i = 0; i < Nx; i++){
            val = matrixA.aDiag[matrixA.idx(i,j)];
            matrixStringDiag.append(" " + std::to_string(val) + " ");

            val = matrixA.aPlusI[matrixA.idx(i,j)];
            matrixStringPlusI.append(" " + std::to_string(val) + " ");

            val = matrixA.aPlusJ[matrixA.idx(i,j)];
            matrixStringPlusJ.append(" " + std::to_string(val) + " ");

        }
        matrixStringDiag.append(")\n");
        matrixStringPlusI.append(")\n");
        matrixStringPlusJ.append(")\n");
    }

    std::cout << matrixStringDiag + "\n\n" + matrixStringPlusI + "\n\n" + matrixStringPlusJ << std::endl;

}
