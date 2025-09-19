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

    std::vector<double> u, v, u_tmp, v_tmp;
    std::vector<double> p, div, dye, dye_tmp;

    inline int idxP(int i, int j) const { return i + Nx*j; }
    inline int idxU(int i, int j) const { return i + (Nx+1)*j; }
    inline int idxV(int i, int j) const { return i + Nx*j; }

    FluidSim(int Nx_, int Ny_, double dx_);

    static double clamp(double x, double a, double b);

    double sampleScalar(const std::vector<double>& s, double x, double y) const;

    double sampleU(const std::vector<double>& uf, double x, double y) const;

    double sampleV(const std::vector<double>& vf, double x, double y) const;

    void sampleVelocity(double x, double y, double& ux, double& vy) const;

    void applyBoundary();

    void computeDivergence();

    void solvePressure();

    void projectVelocity();

    void advectVelocity();

    void advectDye();

    void addForces();

    void addInflow();

    void adaptDt();

    void step();

};