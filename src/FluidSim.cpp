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
        rho = 0.5; 
        gravity = 9.8;
        inflowVelocity = 260.0; // Average top end jet
        // inflowVelocity = 10.0;
        cfl = 1.0;
        dt = cfl * dx / (inflowVelocity / 2);
        std::cout << "dt: " << dt << std::endl;
        poisson_iters = 500;
        tol = 1e-8;
        tuningConst = 0.999;
        safetyConst = 0.30;
        

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
        z = std::vector<double>(Nx*Ny, 0.0);
        q = std::vector<double>(Nx*Ny, 0.0);
        precon = std::vector<double>(Nx*Ny, 0.0);
        
        dye = std::vector<double>(Nx*Ny, 0.0); 
        dye_tmp = std::vector<double>(Nx*Ny, 0.0);

        cellState = std::vector<bool>(Nx*Ny, true);

        matrixA = MatrixA(Nx, Ny, dt, rho, dx, cellState);
        createMICPreconditioner();

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


// void FluidSim::applyBoundary() {
//         for (int j = 0; j <= Ny; ++j) {
//             for (int i = 0; i < Nx; ++i) {

//                 if (i == Nx){
//                     u[idxU(Nx,j)] = u[idxU(Nx-1,j)];
//                     v[idxV(Nx,j)] = v[idxV(Nx-1,j)];
//                 } else if (i == 0){
//                     u[idxU(0, j)]  = 0.0;
//                 } else if (j == 0) {
//                     v[idxV(i, 0)]  = 0.0;
//                 } else if (j == Ny){
//                     v[idxV(i, Ny)] = 0.0;
//                 } else {
//                     setWallFaces(i,j);
//                 }
                
//             }
//         }
        
// }

void FluidSim::applyBoundary() {
    // --- U velocities (vertical faces) ---
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i <= Nx; ++i) {
            // Left wall (no-slip)
            if (i == 0) {
                u[idxU(i,j)] = 0.0;
            }
            // Right outflow (zero-gradient)
            else if (i == Nx) {
                u[idxU(i,j)] = u[idxU(i-1,j)];
            }
            // Internal faces next to solids
            else {
                if (!isFluid(i-1,j) || !isFluid(i,j)) u[idxU(i,j)] = 0.0;
            }
        }
    }

    // --- V velocities (horizontal faces) ---
    for (int j = 0; j <= Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            // Bottom wall (no-slip)
            if (j == 0) {
                v[idxV(i,j)] = 0.0;
            }
            // Top wall (no-slip)
            else if (j == Ny) {
                v[idxV(i,j)] = 0.0;
            }
            // Internal faces next to solids
            else {
                if (!isFluid(i,j-1) || !isFluid(i,j)) v[idxV(i,j)] = 0.0;
            }
        }
    }

}


void FluidSim::computeRHS() {
    const double scale = 1.0 / dx;  // divergence scale
    for (int j = 0; j <= Ny; ++j) {
        for (int i = 0; i <= Nx; ++i) {
            if (!isFluid(i,j)) {
                rhs[idxP(i,j)] = 0.0;
                continue;
            }

            double div = 0.0;
                if(isFluid(i,j)){
                        // Right face
                        if (i == Nx) {
                            // Free surface: treat outside as air → no neighbor fluid velocity
                            div += u[idxU(i-1,j)]; // u(i+1,j) = 0 because air has no resistance
                        }
                        else if (isFluid(i+1,j)) {
                            div += u[idxU(i+1,j)];
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



void FluidSim::cgPressureSolver() {
    const double scale = dt / (rho * dx * dx);
    const int N = Nx * Ny;
    double alpha, beta, sigma, newSigma;

    // Before CG solve:
    for (int i = 0; i < rhs.size(); ++i)
        rhs[i] *= (rho / dt); 

    // --- Early-out if nearly divergence-free ---
    if (maxAbs(rhs) < 1e-10) return;

    // --- Initial setup (pGuess preallocated) ---
    std::fill(pGuess.begin(), pGuess.end(), 0.0);

    // Compute initial residual (r = b - A*p)
    computeResiduals(pGuess, rhs);

    // Initialize auxiliary vectors
    applyPreConditioner(residual, auxZ);
    searchS = auxZ;

    sigma = dotProduct(auxZ, residual);

    double relTol = 1e-6;
    double absTol = 1e-10;
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

        // --- Compute L2 norm of residual ---
        double resNorm = std::sqrt(dotProduct(residual, residual));
        double relRes = resNorm / (r0 + 1e-20);

        // --- Convergence check ---
        if (resNorm <= absTol || relRes <= relTol) {
        std::cout << "CG converged in " << iters << "\n";
        p = pGuess;
        return;
        }

        // --- Compute new direction ---
        applyPreConditioner(residual, auxZ);
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


// void FluidSim::projectVelocity() {
//     const double grad_scale = dt / (rho * dx);

//     // --- Solve pressure Poisson ---
//     computeRHS();
//     matrixA.createAMatrices(Nx, Ny, dt, rho, dx, cellState);
//     createMICPreconditioner();
//     cgPressureSolver();

//     // --- U velocities (vertical faces) ---
//     for (int j = 0; j < Ny; ++j) {
//         for (int i = 0; i <= Nx; ++i) {
//             bool leftFluid  = (i > 0) ? isFluid(i-1,j) : false;
//             bool rightFluid = (i < Nx) ? isFluid(i,j) : false;

//             if (leftFluid && rightFluid) {
//                 u[idxU(i,j)] -= grad_scale * (p[idxP(i,j)] - p[idxP(i-1,j)]);
//             } else if (leftFluid && !rightFluid) {
//                 u[idxU(i,j)] -= grad_scale * (0.0 - p[idxP(i-1,j)]); // right is solid
//             } else if (!leftFluid && rightFluid) {
//                 u[idxU(i,j)] -= grad_scale * (p[idxP(i,j)] - 0.0); // left is solid
//             }

//             // Enforce solid velocity
//             if (!leftFluid || !rightFluid) u[idxU(i,j)] = 0.0;
//         }
//     }

//     // --- V velocities (horizontal faces) ---
//     for (int j = 0; j <= Ny; ++j) {
//         for (int i = 0; i < Nx; ++i) {
//             bool bottomFluid = (j > 0) ? isFluid(i,j-1) : false;
//             bool topFluid    = (j < Ny) ? isFluid(i,j) : false;

//             if (bottomFluid && topFluid) {
//                 v[idxV(i,j)] -= grad_scale * (p[idxP(i,j)] - p[idxP(i,j-1)]);
//             } else if (bottomFluid && !topFluid) {
//                 v[idxV(i,j)] -= grad_scale * (0.0 - p[idxP(i,j-1)]);
//             } else if (!bottomFluid && topFluid) {
//                 v[idxV(i,j)] -= grad_scale * (p[idxP(i,j)] - 0.0);
//             }

//             if (!bottomFluid || !topFluid) v[idxV(i,j)] = 0.0;
//         }
//     }

//     applyBoundary();
// }


void FluidSim::projectVelocity() {
    applyBoundary();
    const double grad_scale = dt / (rho * dx);

    computeRHS();
    matrixA.createAMatrices(Nx, Ny, dt, rho, dx, cellState);
    createMICPreconditioner();
    cgPressureSolver();

    // --- U velocities ---
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i <= Nx; ++i) {
            bool leftFluid  = (i > 0) ? isFluid(i-1,j) : false;
            bool rightFluid = (i < Nx) ? isFluid(i,j) : false;

            if (i == Nx) {
                // Outflow: keep u at zero-gradient
                u[idxU(i,j)] = u[idxU(i-1,j)];
            } else if (leftFluid && rightFluid) {
                u[idxU(i,j)] -= grad_scale * (p[idxP(i,j)] - p[idxP(i-1,j)]);
            } else if (!leftFluid || !rightFluid) {
                u[idxU(i,j)] = 0.0; // stationary solid
            }
        }
    }

    // --- V velocities ---
    for (int j = 0; j <= Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            bool bottomFluid = (j > 0) ? isFluid(i,j-1) : false;
            bool topFluid    = (j < Ny) ? isFluid(i,j) : false;

            if (bottomFluid && topFluid) {
                v[idxV(i,j)] -= grad_scale * (p[idxP(i,j)] - p[idxP(i,j-1)]);
            } else if (!bottomFluid || !topFluid) {
                v[idxV(i,j)] = 0.0; // solid wall
            }
        }
    }

    applyBoundary();
}



// void FluidSim::projectVelocity() {
//     const double grad_scale = dt / (rho * dx);

//     // --- Solve pressure Poisson ---
//     computeRHS();
//     matrixA.createAMatrices(Nx, Ny, dt, rho, dx, cellState);
//     createMICPreconditioner();
//     cgPressureSolver();

//     // --- Update U velocities (vertical faces) ---
//     for (int j = 0; j < Ny; ++j) {
//         for (int i = 0; i <= Nx; ++i) {
//             // Check neighboring cells
//             bool leftFluid  = (i > 0)     ? isFluid(i - 1, j) : false;
//             bool rightFluid = (i < Nx)    ? isFluid(i, j)     : false;

//             if (leftFluid && rightFluid) {
//                 // Face between two fluid cells → apply gradient
//                 u[idxU(i,j)] -= grad_scale * (p[idxP(i,j)] - p[idxP(i-1,j)]);
//             } else {
//                 // Solid boundary or outside domain → enforce solid velocity (stationary)
//                 u[idxU(i,j)] = 0.0;
//             }
//         }
//     }

//     // --- Update V velocities (horizontal faces) ---
//     for (int j = 0; j <= Ny; ++j) {
//         for (int i = 0; i < Nx; ++i) {
//             bool bottomFluid = (j > 0)     ? isFluid(i, j-1) : false;
//             bool topFluid    = (j < Ny)    ? isFluid(i, j)   : false;

//             if (bottomFluid && topFluid) {
//                 // Face between two fluid cells → apply gradient
//                 v[idxV(i,j)] -= grad_scale * (p[idxP(i,j)] - p[idxP(i,j-1)]);
//             } else {
//                 // Solid boundary or outside domain → enforce solid velocity
//                 v[idxV(i,j)] = 0.0;
//             }
//         }
//     }

//     // --- Enforce domain boundaries / walls ---
//     applyBoundary();
// }



// void FluidSim::advectVelocity() {
//     applyBoundary();
//     // U velocities (vertical faces)
//     for (int j = 0; j < Ny; ++j) {
//         for (int i = 0; i <= Nx; ++i) {
//             // Skip faces that are fully inside solids
//             if (i > 0 && i < Nx && (!isFluid(i-1,j) && !isFluid(i,j))) {
//                 u_tmp[idxU(i,j)] = 0.0;
//                 continue;
//             }

//             double x = i*dx;
//             double y = (j + 0.5)*dx;

//             double ux, vy;
//             sampleVelocity(x, y, ux, vy);

//             // Backtrace position
//             double x0 = x - dt*ux;
//             double y0 = y - dt*vy;

//             // Clamp to domain
//             x0 = clamp(x0, 0.0, Nx*dx);
//             y0 = clamp(y0, 0.0, Ny*dx);

//             // If backtrace lands in solid, reset to original position
//             int ix = std::min(int(x0/dx), Nx-1);
//             int iy = std::min(int(y0/dx), Ny-1);
//             if (!isFluid(ix, iy)) {
//                 x0 = x;
//                 y0 = y;
//             }

//             u_tmp[idxU(i,j)] = sampleU(u, x0, y0);
//         }
//     }

//     // V velocities (horizontal faces)
//     for (int j = 0; j <= Ny; ++j) {
//         for (int i = 0; i < Nx; ++i) {
//             if (j > 0 && j < Ny && (!isFluid(i,j-1) && !isFluid(i,j))) {
//                 v_tmp[idxV(i,j)] = 0.0;
//                 continue;
//             }

//             double x = (i + 0.5)*dx;
//             double y = j*dx;

//             double ux, vy;
//             sampleVelocity(x, y, ux, vy);

//             double x0 = x - dt*ux;
//             double y0 = y - dt*vy;

//             x0 = clamp(x0, 0.0, Nx*dx);
//             y0 = clamp(y0, 0.0, Ny*dx);

//             int ix = std::min(int(x0/dx), Nx-1);
//             int iy = std::min(int(y0/dx), Ny-1);
//             if (!isFluid(ix, iy)) {
//                 x0 = x;
//                 y0 = y;
//             }

//             v_tmp[idxV(i,j)] = sampleV(v, x0, y0);
//         }
//     }

//     u.swap(u_tmp);
//     v.swap(v_tmp);

//     // Apply solid/no-slip boundaries
//     applyBoundary();
// }



void FluidSim::advectVelocity() {
    applyBoundary(); // enforce walls before advection

    // --- U velocities (vertical faces) ---
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i <= Nx; ++i) {
            // Skip faces that are fully inside solids
            if (i > 0 && i < Nx && (!isFluid(i-1,j) && !isFluid(i,j))) {
                u_tmp[idxU(i,j)] = 0.0;
                continue;
            }

            double x = i * dx;
            double y = (j + 0.5) * dx;

            double ux, vy;
            sampleVelocity(x, y, ux, vy);

            // Backtrace
            double x0 = x - dt * ux;
            double y0 = y - dt * vy;

            // Clamp to domain
            if (i == Nx) {
                // Right boundary → OUTFLOW
                x0 = clamp(x0, 0.0, Nx * dx);
            } else {
                // Normal solid border
                x0 = clamp(x0, 0.0, Nx * dx);
            }

            // Y clamping (still solid)
            y0 = clamp(y0, 0.0, Ny * dx);

            // If backtrace lands in solid, reset to original position
            int ix = std::min(int(x0 / dx), Nx - 1);
            int iy = std::min(int(y0 / dx), Ny - 1);
            if (i < Nx && !isFluid(ix, iy)) {
                double t = 0.25; // 0 = full backtrace, 1 = reset completely
                x0 = (1.0 - t) * x0 + t * x;
                y0 = (1.0 - t) * y0 + t * y;
            }

            u_tmp[idxU(i,j)] = sampleU(u, x0, y0);
        }
    }

    // --- V velocities (horizontal faces) ---
    for (int j = 0; j <= Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            if (j > 0 && j < Ny && (!isFluid(i,j-1) && !isFluid(i,j))) {
                v_tmp[idxV(i,j)] = 0.0;
                continue;
            }

            double x = (i + 0.5) * dx;
            double y = j * dx;

            double ux, vy;
            sampleVelocity(x, y, ux, vy);

            // Backtrace
            double x0 = x - dt * ux;
            double y0 = y - dt * vy;

            // Clamp to domain
            x0 = clamp(x0, 0.0, Nx * dx);
            y0 = clamp(y0, 0.0, Ny * dx);

            int ix = std::min(int(x0 / dx), Nx - 1);
            int iy = std::min(int(y0 / dx), Ny - 1);
            if (!isFluid(ix, iy)) {
                double t = 0.25; // 0 = full backtrace, 1 = reset completely
                x0 = (1.0 - t) * x0 + t * x;
                y0 = (1.0 - t) * y0 + t * y;
            }   

            v_tmp[idxV(i,j)] = sampleV(v, x0, y0);
        }
    }

    // Swap temporary arrays
    u.swap(u_tmp);
    v.swap(v_tmp);

    // Apply solid/no-slip boundaries
    applyBoundary();
}




void FluidSim::advectDye() {
    std::fill(dye_tmp.begin(), dye_tmp.end(), 0.0);

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            if (!isFluid(i,j)) {
                dye_tmp[idxP(i,j)] = 0.0;
                continue;
            }

            double x = (i + 0.5) * dx;
            double y = (j + 0.5) * dx;

            double ux, vy;
            sampleVelocity(x, y, ux, vy);

            // Backtrace position
            double x0 = x - dt * ux;
            double y0 = y - dt * vy;

            // Clamp to domain
            x0 = clamp(x0, 0.0, Nx * dx);
            y0 = clamp(y0, 0.0, Ny * dx);

            // Determine which cell backtrace lands in
            int ix = std::min(int(x0 / dx), Nx - 1);
            int iy = std::min(int(y0 / dx), Ny - 1);

            // If backtrace hits a solid cell, sample current dye instead
            if (!isFluid(ix, iy)) {
                x0 = x;
                y0 = y;
            }

            dye_tmp[idxP(i,j)] = sampleScalar(dye, x0, y0);
        }
    }

    dye.swap(dye_tmp);

    // Clamp dye values to [0,1]
    for (double &s : dye) s = clamp(s, 0.0, 1.0);
}



void FluidSim::addForces() {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                v[idxV(i,j)] += dt * gravity;
            }
        }
        applyBoundary();
}


void FluidSim::addInflow() {
    int i0 = 0;
    int i1 = 1;
    int j_min = Ny/4;
    int j_max = Ny-Ny/4;

    for (int j = j_min; j < j_max; ++j) {
        for (int i = i0; i <= i1; ++i) {
            int pIdx = idxP(i, j);
            dye[pIdx] = 1.0;

            // Add horizontal inflow velocity to multiple u-cells
            if (i < Nx - 1)
                u[idxU(i, j)] = inflowVelocity;
                // v[idxV(i, j)] = 0.00005;
        }
    }
    applyBoundary();
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

    projectVelocity();
    printMaxDivergence();

    advectVelocity();
    advectDye();
}



bool FluidSim::isFluid(int i, int j){
        if(!inBounds(i,j)){
                return false;
        }
        return cellState[i + (Nx * j)];
}



void FluidSim::setWallFaces(int i, int j){

        if(!inBounds(i,j)) return;

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


void FluidSim::createMICPreconditioner() {
    const int N = Nx * Ny;
    double Ad, e, term, val;
    double A_x, A_yL, pL;
    double A_y, A_xD, pD;
    std::fill(precon.begin(), precon.end(), 0.0);

    // Build precon in row-major order (i=0..Nx-1, j=0..Ny-1)
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int idx = idxP(i, j);
            if (!cellState[idx]) {
                precon[idx] = 0.0;
                continue;
            }

            Ad = matrixA.aDiag[idx];    // diag A(i,i)
            e = Ad;

            // Left neighbor (i-1, j)
            if (i > 0) {
                int idxL = idxP(i - 1, j);
                if (cellState[idxL]) {
                    A_x = matrixA.aPlusI[idxL];   // A(i, i-1) stored at left cell
                    pL  = precon[idxL];           // precon(i-1,j)
                    // subtract (A_x * pL)^2
                    term = A_x * pL;
                    e -= term * term;

                    // cross-term using left cell's AplusJ
                    A_yL = matrixA.aPlusJ[idxL];   // Aplusj(i-1,j)
                    if (A_yL != 0.0) {
                        e -= tuningConst * (A_x * A_yL) * (pL * pL);
                    }
                }
            }

            // Down neighbor (i, j-1)
            if (j > 0) {
                int idxD = idxP(i, j - 1);
                if (cellState[idxD]) {
                    A_y = matrixA.aPlusJ[idxD];   // A(i, i-nx) stored at down cell
                    pD  = precon[idxD];
                    term = A_y * pD;
                    e -= term * term;

                    // cross-term using down cell's AplusI
                    A_xD = matrixA.aPlusI[idxD];
                    if (A_xD != 0.0) {
                        e -= tuningConst * (A_y * A_xD) * (pD * pD);
                    }
                }
            }

            // Safety clamp: ensure e stays positive and not too small
            if (!(e > 0.0) || e < safetyConst * Ad) {
                e = Ad;
            }

            // store 1/sqrt(e)
            val = 1.0 / std::sqrt(e);
            if (!std::isfinite(val)) {
                // fallback: avoid NaN/Inf
                val = 1.0 / std::sqrt(std::max(Ad, 1e-12));
            }
            precon[idx] = val;
        }
    }
}



void FluidSim::applyPreConditioner(const std::vector<double>& r, std::vector<double>& z) {
    int idx;
    double t;
    double L_li, L_di;
    double LT_ir, LT_iu;
    // q temporary (forward solve result)
    std::fill(q.begin(), q.end(), 0.0);
    std::fill(z.begin(), z.end(), 0.0);

    // Forward solve: L * q = r
    // row-major: j = 0..Ny-1, i = 0..Nx-1
    int idxL, idxD;
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            idx = idxP(i, j);

            if (!cellState[idx]) continue;

            t = r[idx];

            // contribution from left neighbor (i-1, j)
            if (i > 0) {
                idxL = idxP(i - 1, j);
                if (cellState[idxL]) {
                    // L(i, i-1) = A(i-1 -> i) * precon(i-1)
                    L_li = matrixA.aPlusI[idxL] * precon[idxL];
                    t -= L_li * q[idxL];
                }
            }

            // contribution from down neighbor (i, j-1)
            if (j > 0) {
                idxD = idxP(i, j - 1);
                if (cellState[idxD]) {
                    L_di = matrixA.aPlusJ[idxD] * precon[idxD];
                    t -= L_di * q[idxD];
                }
            }

            // diagonal scaling
            q[idx] = t * precon[idx];
        }
    }

    // Backward solve: L^T * z = q
    // reverse row-major: j = Ny-1..0, i = Nx-1..0
    int idxR, idxU;
    for (int jj = Ny - 1; jj >= 0; --jj) {
        for (int ii = Nx - 1; ii >= 0; --ii) {
            idx = idxP(ii, jj);
            if (!cellState[idx]) continue;

            t = q[idx];

            // contribution from right neighbor (ii+1, jj)
            if (ii + 1 < Nx) {
                idxR = idxP(ii + 1, jj);
                if (cellState[idxR]) {
                    // L(ii, idxR) in transpose uses L(ii, idx) stored at idx
                    // L(ii,idxR) = A(ii -> ii+1) * precon[ii]
                    LT_ir = matrixA.aPlusI[idx] * precon[idx];
                    t -= LT_ir * z[idxR];
                }
            }

            // contribution from up neighbor (ii, jj+1)
            if (jj + 1 < Ny) {
                idxU = idxP(ii, jj + 1);
                if (cellState[idxU]) {
                    LT_iu = matrixA.aPlusJ[idx] * precon[idx];
                    t -= LT_iu * z[idxU];
                }
            }

            z[idx] = t * precon[idx];
        }
    }
}




double FluidSim::dotProduct(const std::vector<double>& a, const std::vector<double>& b){
    double sum = 0.0;
    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            int index = idxP(i,j);
            if (isFluid(i,j)){
                sum += a[index] * b[index];
            }            
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

    // Determine a reference velocity scale (could be max of u/v or inflow)
    double Uref = 0.0;
    for (double val : u) Uref = std::max(Uref, std::abs(val));
    for (double val : v) Uref = std::max(Uref, std::abs(val));
    Uref = std::max(Uref, 1e-8); // avoid division by zero

    for (int j = 0; j < Ny; ++j) {
        for (int i = 0; i < Nx; ++i) {
            if (!isFluid(i,j)) continue;

            double du_dx = u[idxU(i+1,j)] - u[idxU(i,j)];
            double dv_dy = v[idxV(i,j+1)] - v[idxV(i,j)];

            double div = (du_dx + dv_dy) * scale; // raw divergence

            // normalized divergence: non-dimensional
            double normDiv = std::abs(div) * dx / Uref;

            maxDiv = std::max(maxDiv, normDiv);
        }
    }

    std::cout << "Max normalized divergence: " << maxDiv << std::endl;
}

