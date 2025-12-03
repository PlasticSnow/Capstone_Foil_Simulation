#include "Tests.hpp"


void testMatrixASetup(int Nx, int Ny, double dt, double rho, double dx){

    std::vector<bool> fluidState = std::vector<bool>(Nx*Ny, true);

    for (int j = 1; j < Ny; j++){
        for (int i = 0; i < Nx; i++){
            if (i == 0 && j < Ny - 1){
                
                fluidState[i + Nx*j] = false;
            }
        }
    }

    for (int j = 0; j < Ny; j++){
        std::cout << "[ ";
        for (int i = 0; i < Nx; i++){
            std::cout << fluidState[i + Nx * j] << " ";
        }
        std::cout << "]\n";
    }
    

    MatrixA matrixA = MatrixA(Nx, Ny, dt, rho, dx, fluidState);

    std::string matrixStringDiag = ""; 
    std::string matrixStringPlusI = "";
    std::string matrixStringPlusJ = "";   

    std::cout << "About to print" << std::endl;

    for (int j = 0; j < Ny; j++){
        matrixStringDiag.append("[");
        matrixStringPlusI.append("[");
        matrixStringPlusJ.append("[");

        double val;
        for (int i = 0; i < Nx; i++){
            val = matrixA.aDiag[matrixA.idx(i,j)];
            matrixStringDiag.append(" " + std::to_string(val) + " ");

            val = matrixA.aPlusI[matrixA.idx(i,j)];
            matrixStringPlusI.append(" " + std::to_string(val) + " ");

            val = matrixA.aPlusJ[matrixA.idx(i,j)];
            matrixStringPlusJ.append(" " + std::to_string(val) + " ");

        
            matrixStringDiag.append(")\n");
            matrixStringPlusI.append(")\n");
            matrixStringPlusJ.append(")\n");
        }
    }


    printDiag(matrixA.aDiag, Nx, Ny);

    // std::cout << matrixStringDiag + "\n\n" + matrixStringPlusI + "\n\n" + matrixStringPlusJ << std::endl;

}



void printDiag(std::vector<double> aDiag, int Nx, int Ny){
    std::string aDiagString = "";
    int columns = Nx * Ny;
    
    for (int c = 0; c < columns; c++){
        aDiagString.append("[ ");
        for (int i = 0; i < columns; i++){
            if (i == c){
                aDiagString.append(doub_to_string(std::to_string(aDiag[c])) + " ");
            } else {
                aDiagString.append("0 ");
            }
        }
        aDiagString.append("]\n");
    }
    
    std::cout << aDiagString << "\n\n";
}



std::string doub_to_string(std::string val){
    std::string formattedString = "";
    size_t s_size = val.size();
    int i = 0;

    for (int i = 0; i < s_size; i++){
        if (val.at(i) != '.'){
            formattedString += val.at(i);
        } else {
            break;
        }
    }

    return formattedString; 
}