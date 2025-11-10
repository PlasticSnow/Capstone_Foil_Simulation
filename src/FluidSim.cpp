#include "FluidSim.hpp"

inline double maxAbs(const std::vector<double>& v) {
    double maxVal = 0.0;
    for (double x : v)
        maxVal = std::max(maxVal, std::fabs(x));
    return maxVal;
}


FluidSim::FluidSim(int Nx_, int Ny_, double dx_){
        Nx = Nx_;
        Ny = Ny_; 
        dx = dx_; 
        rho = 1.0; 
        gravity = -9.8;
        cfl = 1.0;
        dt = 0.1;
        poisson_iters = 500;
        tol = 1e-4;

        u = std::vector<double>((Nx+1)*Ny, 0.0);
        v = std::vector<double>(Nx*(Ny+1), 0.0);
        
        u_tmp = std::vector<double>((Nx+1)*Ny, 0.0);
        v_tmp = std::vector<double>(Nx*(Ny+1), 0.0);
        
        p = std::vector<double>(Nx*Ny, 0.0);
        rhs = std::vector<double>(Nx*Ny, 0.0);

        pGuess = std::vector<double>(Nx*Ny, 0.0);
        pTemp = std::vector<double>(Nx*Ny, 0.0);
        auxZ = std::vector<double>(Nx*Ny, 0.0);
        searchS = std::vector<double>(Nx*Ny, 0.0);
        residual = std::vector<double>(Nx*Ny, 0.0);
        
        dye = std::vector<double>(Nx*Ny, 0.0); 
        dye_tmp = std::vector<double>(Nx*Ny, 0.0);

        cellState = std::vector<bool>(Nx*Ny, true);

        matrixA = MatrixA(Nx, Ny, dt, rho, dx, cellState);

        calculateTimeStep();
}

double FluidSim::clamp(double x, double a, double b) {
        return std::max(a, std::min(b, x));
}

double FluidSim::sampleScalar(const std::vector<double>& s, double x, double y) const {
        double gx = x/dx - 0.5;
        double gy = y/dx - 0.5;
        gx = clamp(gx, 0.0, Nx-1.001);
        gy = clamp(gy, 0.0, Ny-1.001);
        int i0 = (int)std::floor(gx);
        int j0 = (int)std::floor(gy);
        int i1 = std::min(i0+1, Nx-1);
        int j1 = std::min(j0+1, Ny-1);
        double tx = gx - i0;
        double ty = gy - j0;
        double s00 = s[idxP(i0,j0)];
        double s10 = s[idxP(i1,j0)];
        double s01 = s[idxP(i0,j1)];
        double s11 = s[idxP(i1,j1)];
        double a = s00*(1-tx) + s10*tx;
        double b = s01*(1-tx) + s11*tx;
        return a*(1-ty) + b*ty;
}

double FluidSim::sampleU(const std::vector<double>& uf, double x, double y) const {
        double gx = x/dx;
        double gy = y/dx - 0.5;
        gx = clamp(gx, 0.0, Nx - 0.001);
        gy = clamp(gy, 0.0, Ny-1 - 0.001);
        int i0 = (int)std::floor(gx);
        int j0 = (int)std::floor(gy);
        int i1 = std::min(i0+1, Nx);
        int j1 = std::min(j0+1, Ny-1);
        double tx = gx - i0;
        double ty = gy - j0;
        double s00 = uf[idxU(i0,j0)];
        double s10 = uf[idxU(i1,j0)];
        double s01 = uf[idxU(i0,j1)];
        double s11 = uf[idxU(i1,j1)];
        double a = s00*(1-tx) + s10*tx;
        double b = s01*(1-tx) + s11*tx;
        return a*(1-ty) + b*ty;
}


double FluidSim::sampleV(const std::vector<double>& vf, double x, double y) const {
        double gx = x/dx - 0.5;
        double gy = y/dx;
        gx = clamp(gx, 0.0, Nx-1 - 0.001);
        gy = clamp(gy, 0.0, Ny - 0.001);
        int i0 = (int)std::floor(gx);
        int j0 = (int)std::floor(gy);
        int i1 = std::min(i0+1, Nx-1);
        int j1 = std::min(j0+1, Ny);
        double tx = gx - i0;
        double ty = gy - j0;
        double s00 = vf[idxV(i0,j0)];
        double s10 = vf[idxV(i1,j0)];
        double s01 = vf[idxV(i0,j1)];
        double s11 = vf[idxV(i1,j1)];
        double a = s00*(1-tx) + s10*tx;
        double b = s01*(1-tx) + s11*tx;
        return a*(1-ty) + b*ty;
}


void FluidSim::sampleVelocity(double x, double y, double& ux, double& vy) const {
        ux = sampleU(u, x, y);
        vy = sampleV(v, x, y);
}


void FluidSim::applyBoundary() {
        for (int j = 0; j < Ny; ++j) {
                u[idxU(0, j)]  = 0.0;
                u[idxU(Nx,j)]  = 0.0;
        }
        for (int i = 0; i < Nx; ++i) {
                v[idxV(i, 0)]  = 0.0;
                v[idxV(i, Ny)] = 0.0;
        }
}


void FluidSim::computeRHS() {
    const double scale = 1.0 / dx;  // divergence scale
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            if (!isFluid(i,j)) {
                rhs[idxP(i,j)] = 0.0;
                continue;
            }

            double div = 0.0;
                if(isFluid(i,j)){
                        // Right face
                        if (i + 1 < Nx) {
                                if (isFluid(i+1,j)) div += u[idxU(i+1,j)];
                        }
                        // Left face
                        if (i >= 0) {
                                if (isFluid(i-1,j)) div -= u[idxU(i,j)];
                        }

                        // Top face
                        if (j + 1 < Ny) {
                                if (isFluid(i,j+1)) div += v[idxV(i,j+1)];
                        }
                        // Bottom face
                        if (j >= 0) {
                                if (isFluid(i,j-1)) div -= v[idxV(i,j)];
                        }

                        rhs[idxP(i,j)] = div * -scale;
                }
            
        }
    }
}


// void FluidSim::cgPressureSolver(){
//         double scale = dt / (rho*dx*dx);
//         int iteration = 0;
//         double alpha;

//         // Set initial guess p = 0
//         std::fill(pGuess.begin(), pGuess.end(), 0);

//         // Set residual vector
//         computeResiduals(pGuess, rhs);

//         // Set auxilery vector z and search vector s = z
//         auxZ = residual;
//         searchS = auxZ;

//         //
//         double sigma = dotProduct(auxZ, residual);

//         // Loop until solved or iterations exceeded
//         for (int iters = 0; iters < poisson_iters; iters++){

//                 // Set auxiliary vector z to applyA(s)
//                 matrixA.applyA(searchS, auxZ, cellState);

//                 // Calculate denom for alpha and make sure no divide by zero
//                 double denom = dotProduct(searchS, auxZ);
//                 if (fabs(denom) < 1e-20) break; 

//                 // Calculate alpha
//                 alpha = sigma/denom;

//                 int index;
//                 double maxR = 0;
//                 for (int j = 0; j < Ny; j++){
//                         for (int i = 0; i < Nx; i++){
//                                 index = idxP(i,j);

//                                 // If solid skip
//                                 if (!cellState[index]) continue;

//                                 // Update guess pressures
//                                 pGuess[index] += alpha * searchS[index];
                                

//                                 // Update residuals
//                                 residual[index] -= alpha * auxZ[index];

//                                 maxR = std::max(maxR, fabs(residual[index]));
//                         }
//                 }

//                 // If residual lower than tol update pressure with pGuess and end
//                 if (maxR <= tol){
//                         std::cout << "CG converged in " << iters << " iterations\n";
//                         p = pGuess;
//                         return;
//                 }
                
//                 auxZ = residual;

//                 double newSigma = dotProduct(auxZ, residual);

//                 double beta = newSigma / sigma;

//                 // Set search vector s = z + Bs
//                 for (int i = 0; i < searchS.size(); i++){
//                         searchS[i] = auxZ[i] + (searchS[i] * beta);
//                 }

//                 // Update sigma
//                 sigma = newSigma;

//                 // Reset maxR
//                 maxR = 0.0;
//         }

//         // If max iterations exceeded
//         std::cout << "Max Iterations reached " + std::to_string(poisson_iters) << std::endl;
//         p = pGuess;
        
// }


void FluidSim::cgPressureSolver() {
    const double scale = dt / (rho * dx * dx);
    const int N = Nx * Ny;
    double alpha, beta, sigma, newSigma;

    // --- Early-out if nearly divergence-free ---
    if (maxAbs(rhs) < 1e-8) return;

    // --- Initial setup (pGuess preallocated) ---
    std::fill(pGuess.begin(), pGuess.end(), 0.0);

    // Compute initial residual (r = b - A*p)
    computeResiduals(pGuess, rhs);

    // Initialize auxiliary vectors
    auxZ = residual;    // Preconditioner could go here if used
    searchS = auxZ;

    sigma = dotProduct(auxZ, residual);

    double relTol = 1e-4;
    double absTol = 1e-8;
    double r0 = std::sqrt(sigma);
    double maxR = 0.0;

    for (int iters = 0; iters < poisson_iters; ++iters) {
        // --- A*s step ---
        matrixA.applyA(searchS, auxZ, cellState);

        // --- Compute alpha ---
        double denom = dotProduct(searchS, auxZ);
        if (fabs(denom) < 1e-20) break;
        alpha = sigma / denom;

        // --- Update p and residual ---
        maxR = 0.0;

        // Unroll inner loop, no branching in hot path
        #pragma omp simd reduction(max:maxR)
        for (int i = 0; i < N; ++i) {
            if (!cellState[i]) continue;
            pGuess[i] += alpha * searchS[i];
            residual[i] -= alpha * auxZ[i];
            maxR = std::max(maxR, fabs(residual[i]));
        }

        // --- Convergence check ---
        double relRes = maxR / (r0 + 1e-20);
        if (maxR <= absTol || relRes <= relTol) {
            std::cout << "CG converged in " << iters << " iterations\n";
            p = pGuess;
            return;
        }

        // --- Compute new direction ---
        auxZ = residual;  // (if preconditioned, apply M^-1 here)
        newSigma = dotProduct(auxZ, residual);
        beta = newSigma / sigma;

        // Update search direction: s = z + beta*s
        #pragma omp simd
        for (int i = 0; i < N; ++i) {
            searchS[i] = auxZ[i] + beta * searchS[i];
        }

        sigma = newSigma;
    }

    std::cout << "Max iterations reached " << poisson_iters << std::endl;
    p = pGuess;
}




void FluidSim::projectVelocity() {
        const double grad_scale = dt / (rho * dx);


        // U velocities (vertical faces)
        for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i <= Nx; ++i) {
                // Skip boundaries if not fluid on both sides
                if (i == 0 || i == Nx) continue;
                if (isFluid(i-1,j) && isFluid(i,j)) {
                u[idxU(i,j)] -= grad_scale * (p[idxP(i,j)] - p[idxP(i-1,j)]);
                }
        }
        }

        // V velocities (horizontal faces)
        for (int j = 0; j <= Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
                if (j == 0 || j == Ny) continue;
                if (isFluid(i,j-1) && isFluid(i,j)) {
                v[idxV(i,j)] -= grad_scale * (p[idxP(i,j)] - p[idxP(i,j-1)]);
                }
        }
        }

        applyBoundary();

}


void FluidSim::advectVelocity() {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i <= Nx; ++i) {
                double x = i*dx;
                double y = (j + 0.5)*dx;
                double ux, vy;sampleVelocity(x, y, ux, vy);
                double x0 = clamp(x - dt*ux, 0.0, Nx*dx);
                double y0 = clamp(y - dt*vy, 0.0, Ny*dx);
                u_tmp[idxU(i,j)] = sampleU(u, x0, y0);
            }
        }
        for (int j = 0; j <= Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = (i + 0.5)*dx;
                double y = j*dx;
                double ux, vy; sampleVelocity(x, y, ux, vy);
                double x0 = clamp(x - dt*ux, 0.0, Nx*dx);
                double y0 = clamp(y - dt*vy, 0.0, Ny*dx);
                v_tmp[idxV(i,j)] = sampleV(v, x0, y0);
            }
        }
        u.swap(u_tmp);
        v.swap(v_tmp);
        applyBoundary();
    }


void FluidSim::advectDye() {
        std::fill(dye_tmp.begin(), dye_tmp.end(), 0.0);
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                if (isFluid(i,j)){
                        double x = (i + 0.5)*dx;
                        double y = (j + 0.5)*dx;
                        double ux, vy; sampleVelocity(x, y, ux, vy);
                        double x0 = clamp(x - dt*ux, 0.0, Nx*dx);
                        double y0 = clamp(y - dt*vy, 0.0, Ny*dx);
                        dye_tmp[idxP(i,j)] = sampleScalar(dye, x0, y0);
                } else {
                        dye[idxP(i,j)] = 0.0;
                }
            }
        }
        dye.swap(dye_tmp);
        
        for (double &s : dye) s = clamp(s, 0.0, 1.0);
}


void FluidSim::addForces() {
        for (int j = 0; j <= Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                v[idxV(i,j)] += dt * gravity;
            }
        }
        applyBoundary();
}


void FluidSim::addInflow() {
        int i_max = 2;
        for (int j = 1; j < Ny - 1; ++j){
                for (int i = 0; i < i_max; ++i) {
                        dye[idxP(i,j)] = 1.0;
                }
                if (Nx >= 1) {
                        u[idxU(1, j)] = 1.0;
                }
        }
}


void FluidSim::calculateTimeStep() {
        double maxu = 1e-8;
        for (double val : u) maxu = std::max(maxu, std::abs(val));
        for (double val : v) maxu = std::max(maxu, std::abs(val));
        double dt_cfl = cfl * dx / maxu;
        dt = clamp(dt_cfl, 0.001, 0.02);
}


void FluidSim::step() {
        addInflow();
        addForces();
        
        enforceWallFaces();
        advectVelocity();
        enforceWallFaces();

        computeRHS();

        cgPressureSolver();

        projectVelocity();

        printMaxDivergence();

        enforceWallFaces();
        advectDye();
}



bool FluidSim::isFluid(int i, int j){
        if(!inBounds(i,j)){
                return false;
        }
        return cellState[i + (Nx * j)];
}



void FluidSim::setWallFaces(int i, int j){

        if(!inBounds(i,j)) {return;}

        if(!isFluid(i,j)){
                u[idxU(i, j)] = 0.0;
                u[idxU(i + 1,j)] = 0.0;
                v[idxV(i, j)] = 0.0;
                v[idxV(i, j + 1)] = 0.0;
        }

        // Set Left Face
        if (!isFluid(i-1,j)) {u[idxU(i, j)] = 0.0;} 

        // Set Right Face
        if (!isFluid(i+1, j)) {u[idxU(i + 1,j)] = 0.0;}

        // Set Bottom Face
        if (!isFluid(i, j-1)) {v[idxV(i, j)] = 0.0;}

        // Set Top Face
        if (!isFluid(i, j+1)) {v[idxV(i, j + 1)] = 0.0;}
 
}


bool FluidSim::inBounds(int i, int j) {
    return (i >= 0 && i < Nx && j >= 0 && j < Ny);
}



void FluidSim::enforceWallFaces(){
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                setWallFaces(i,j);
            }
        }
}

// Stub pre conditioner function currentlly just returns itself
std::vector<double> FluidSim::preConditioner(std::vector<double> residuals){
        for (double r : residuals){
                r = r * 1;
        }
        return residuals;
}



double FluidSim::dotProduct(const std::vector<double>& a, const std::vector<double>& b) const {
    double sum = 0.0;
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int index = idxP(i, j);
            if (!cellState[index]) continue; // Skip solids
            sum += a[index] * b[index];
        }
    }
    return sum;
}

void FluidSim::computeResiduals(const std::vector<double>& pGuess, const std::vector<double>& rhs){
        matrixA.applyA(pGuess, pTemp, cellState);
        for (int i = 0; i < rhs.size(); i++){
                residual[i] = rhs[i] - pTemp[i];
        }
}



void FluidSim::printMaxDivergence() {
    const double scale = 1.0 / dx;
    double maxDiv = 0.0;
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            if (!isFluid(i,j)) continue;
            double du_dx = u[idxU(i+1,j)] - u[idxU(i,j)];
            double dv_dy = v[idxV(i,j+1)] - v[idxV(i,j)];
            double div = -scale * (du_dx + dv_dy);
            maxDiv = std::max(maxDiv, std::abs(div));
        }
    }
    std::cout << "Max divergence: " << maxDiv << std::endl;
}




