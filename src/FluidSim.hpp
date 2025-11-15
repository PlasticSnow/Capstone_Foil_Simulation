#pragma once
#include <vector>
#include <numeric>
#include <algorithm>
#include <cmath>
#include "CoefficientMatrixA.hpp"

enum Operation {
    Add,
    Sub,
    Mult,
    Divide
};


class FluidSim{

    public:
    int Nx, Ny;         // grid resolution
    double dx;          // cell size
    double rho;         // density
    double gravity;     // gravity (acts on v)
    double cfl;         // CFL safety factor
    double dt;          // current dt (time step)
    int poisson_iters;  // iterations for pressure solve
    double tol;         // tolerance for divergence
    double tuningConst, safetyConst;
    double inflowVelocity;


    // Is Cell Solid
    std::vector<bool> cellState;


    // Staggered grid velocity components (MAC grid)
    // u = horizontal velocity (size: (Nx+1) * Ny)
    // v = vertical velocity (size: Nx * (Ny+1))
    std::vector<double> u, v, u_tmp, v_tmp;
    std::vector<double> p, pGuess, pTemp, auxZ, searchS, residual, precon, q, z;
    std::vector<double> rhs, dye, dye_tmp;
    std::vector<double> a_diag, a_plusI, a_plusJ;

    MatrixA matrixA;


    inline int idxP(int i, int j) const { return i + Nx*j; }
    inline int idxU(int i, int j) const { return i + (Nx+1)*j; }
    inline int idxV(int i, int j) const { return i + Nx*j; }

    FluidSim(int Nx_, int Ny_, double dx_);


    void setGravity(double grav){this->gravity = grav;}


    static double clamp(double x, double a, double b);


    double sampleScalar(const std::vector<double>& s, double x, double y) const;


    double sampleU(const std::vector<double>& uf, double x, double y) const;


    double sampleV(const std::vector<double>& vf, double x, double y) const;


    // Helper: sample full velocity vector at (x,y)
    void sampleVelocity(double x, double y, double& ux, double& vy) const;


    // Enforce no-slip boundary conditions (walls) ---------------------
    void applyBoundary();


    // Compute divergence of velocity field (for incompressibility) ------------    
    void computeRHS();


    // Solve pressure Poisson equation ∇²p = (ρ/Δt) ∇·u ---------------
    void cgPressureSolver();

    // Subtract pressure gradient from velocity field (projection) -----
    void projectVelocity();


    // Semi-Lagrangian advection of velocity fields ---------------------------
    void advectVelocity();


    // Semi-Lagrangian advection of dye scalar --------------------------------
    void advectDye();


    // Add external forces (gravity on v) --------------------------------------
    void addForces();


    // Constant inflow of dye + velocity at left boundary ----------------------
    void addInflow();


    // Adapt time step based on max velocity (CFL condition) -------------------
    void calculateTimeStep();

    /** 
     * @brief Contains all the functions needed in order to update the simulation
     * @details This function contains all the functions needed in order
     *          to update the simulation, it would be called for each time step.
     *          Functions in order as follows:
     *              addInflow();
     *              adaptDt();
     *              advectVelocity();
     *              addForces();
     *              computeDivergence();
     *              solvePressure();
     *              projectVelocity();
     *              advectDye();
     */
    void step();


    void setCellState(int i, int j, bool state){this->cellState[i + (Nx * j)] = state;}

    bool isFluid(int i, int j);

    void setWallFaces(int i, int j);

    bool inBounds(int i, int j);

    void enforceWallFaces();

    void createMICPreconditioner();

    void applyPreConditioner(const std::vector<double>& r, std::vector<double>& z) ;

    double dotProduct(const std::vector<double>& a, const std::vector<double>& b) const;

    // void applyScalar(std::vector<double>& v, double scalar, Operation scalarOp);

    void computeResiduals(const std::vector<double>& pGuess, const std::vector<double>& rhs);

    void printMaxDivergence();

};