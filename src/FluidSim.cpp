#include "FluidSim.hpp"


FluidSim::FluidSim(int Nx_, int Ny_, double dx_){
        Nx = Nx_;
        Ny = Ny_; 
        dx = dx_; 
        rho = 1.0; 
        gravity = -9.8;
        cfl = 0.5;
        dt = 0.01;
        poisson_iters = 200;
        u = std::vector<double>((Nx+1)*Ny, 0.0);
        v = std::vector<double>(Nx*(Ny+1), 0.0);
        u_tmp = std::vector<double>((Nx+1)*Ny, 0.0);
        v_tmp = std::vector<double>(Nx*(Ny+1), 0.0);
        p = std::vector<double>(Nx*Ny, 0.0);
        div = std::vector<double>(Nx*Ny, 0.0);
        dye = std::vector<double>(Nx*Ny, 0.0); 
        dye_tmp = std::vector<double>(Nx*Ny, 0.0);
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


void FluidSim::computeDivergence() {
        const double invdx = 1.0/dx;
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double du_dx = (u[idxU(i+1,j)] - u[idxU(i,j)]) * invdx;
                double dv_dy = (v[idxV(i,j+1)] - v[idxV(i,j)]) * invdx;
                div[idxP(i,j)] = du_dx + dv_dy;
            }
        }
}


void FluidSim::solvePressure() {
        const double scale = rho / dt;
        std::fill(p.begin(), p.end(), 0.0);
        for (int it = 0; it < poisson_iters; ++it) {
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    double pL = (i > 0   ) ? p[idxP(i-1,j)] : p[idxP(i,j)];
                    double pR = (i < Nx-1) ? p[idxP(i+1,j)] : p[idxP(i,j)];
                    double pB = (j > 0   ) ? p[idxP(i,j-1)] : p[idxP(i,j)];
                    double pT = (j < Ny-1) ? p[idxP(i,j+1)] : p[idxP(i,j)];
                    double rhs = scale * div[idxP(i,j)];
                    p[idxP(i,j)] = (pL + pR + pB + pT - rhs * dx*dx) * 0.25;
                }
            }
        }
}


void FluidSim::projectVelocity() {
        const double grad_scale = dt / (rho * dx);
        for (int j = 0; j < Ny; ++j) {
            for (int i = 1; i < Nx; ++i) {
                double g = p[idxP(i,j)] - p[idxP(i-1,j)];
                u[idxU(i,j)] -= grad_scale * g;
            }
        }
        for (int j = 1; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double g = p[idxP(i,j)] - p[idxP(i,j-1)];
                v[idxV(i,j)] -= grad_scale * g;
            }
        }
        applyBoundary();
    }


void FluidSim::advectVelocity() {
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i <= Nx; ++i) {
                double x = i*dx;
                double y = (j + 0.5)*dx;
                double ux, vy; sampleVelocity(x, y, ux, vy);
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
        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double x = (i + 0.5)*dx;
                double y = (j + 0.5)*dx;
                double ux, vy; sampleVelocity(x, y, ux, vy);
                double x0 = clamp(x - dt*ux, 0.0, Nx*dx);
                double y0 = clamp(y - dt*vy, 0.0, Ny*dx);
                dye_tmp[idxP(i,j)] = sampleScalar(dye, x0, y0);
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
        int i_max = std::min(3, Nx-1);
        for (int j = Ny/3; j < 2*Ny/3; ++j) {
                for (int i = 0; i < i_max; ++i) {
                dye[idxP(i,j)] = 1.0;
                }
                if (Nx >= 2) {
                u[idxU(1, j)] = 2.0;
                }
        }
}


void FluidSim::adaptDt() {
        double maxu = 1e-8;
        for (double val : u) maxu = std::max(maxu, std::abs(val));
        for (double val : v) maxu = std::max(maxu, std::abs(val));
        double dt_cfl = cfl * dx / maxu;
        dt = clamp(dt_cfl, 0.001, 0.02);
}


void FluidSim::step() {
        addInflow();
        adaptDt();
        advectVelocity();
        addForces();
        computeDivergence();
        solvePressure();
        projectVelocity();
        advectDye();
}