#include <vector>

#pragma include


class FluidSim{

    public:
    int Nx, Ny;         // grid resolution
    double dx;          // cell size
    double rho;         // density
    double gravity;     // gravity (acts on v)
    double cfl;         // CFL safety factor
    double dt;          // current dt
    int poisson_iters;  // iterations for pressure solve

    // Is Cell Solid
    std::vector<bool> isSolid;


    // Staggered grid velocity components (MAC grid)
    // u = horizontal velocity (size: (Nx+1) * Ny)
    // v = vertical velocity (size: Nx * (Ny+1))
    std::vector<double> u, v, u_tmp, v_tmp;
    std::vector<double> p, div, dye, dye_tmp;

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
    void computeDivergence();


    // Solve pressure Poisson equation ∇²p = (ρ/Δt) ∇·u ---------------
    void solvePressure();


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
    void adaptDt();

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


    void setSolid(int i, int j, bool state){this->isSolid[i + (Nx * j)] = state;}

    bool getSolid(int i, int j) const {return isSolid[i + (Nx * j)];}

    void setWallFaces(int i, int j);

    bool inBounds(int i, int j) const {return }

};