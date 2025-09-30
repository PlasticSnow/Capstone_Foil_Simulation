// 2D Incompressible Fluid Simulation (Projection Method / MAC Grid) with SFML Visualization
// ---------------------------------------------------------------
// Features
//  - Marker-and-Cell (MAC) staggered grid
//  - Semi-Lagrangian advection for velocity and dye
//  - Pressure Poisson solve (Gauss–Seidel)
//  - Projection step to enforce incompressibility
//  - Simple inflow + dye injection
//  - Interactive SFML window visualization

#include <SFML/Graphics.hpp>
#include <cmath>
#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cassert>
#include "FluidSim.hpp"

int main() {
    int Nx = 50, Ny = 50;
    double dx = 1.0 / Nx;
    FluidSim sim(Nx, Ny, dx);
    sim.setGravity(0);
    sim.setSolid(Nx/2, Ny/2, true);

    int scale = 5; // pixel size per cell
    sf::RenderWindow window(sf::VideoMode(Nx*scale, Ny*scale), "2D Fluid Simulation (Projection)");
    window.setFramerateLimit(60);

    sf::Image img;
    img.create(Nx, Ny, sf::Color::Black);
    sf::Texture tex;
    sf::Sprite spr;
    tex.create(Nx, Ny);
    spr.setScale((float)scale, (float)scale);

    while (window.isOpen()) {
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        sim.step();

        for (int j = 0; j < Ny; ++j) {
            for (int i = 0; i < Nx; ++i) {
                double s = sim.dye[sim.idxP(i,j)];
                unsigned char c = (unsigned char)std::round(sim.clamp(s,0.0,1.0) * 255.0);
                img.setPixel(i, Ny-1-j, sf::Color(c,c,c));
            }
        }
        tex.update(img);
        spr.setTexture(tex);

        window.clear();
        window.draw(spr);
        window.display();
    }
    return 0;
}
