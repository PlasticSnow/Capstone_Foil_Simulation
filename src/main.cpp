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

struct Point {
    double x;
    double y;

    Point(double _x, double _y){
        x = _x;
        y = _y;
    }

};


std::vector<Point> createRectangle(double startX, double startY, int with, int height);

std::vector<Point> createDiagLine(double startX, double startY, int length);

void setMultState(FluidSim &_sim, std::vector<Point>_points);



int main() {
    int Nx = 75, Ny = 75;
    double dx = 1.0 / Nx;

    std::vector<Point> rectangleOne = createRectangle((Nx/2) - 5, (Ny/2) - 5,  5 , 5);
    for (Point point : rectangleOne){
        std::cout << "(" + std::to_string(point.x) + "," + std::to_string(point.y) + ")" << std::endl;
    }

    std::vector<Point> rectangleTwo = createRectangle((Nx/2) - (Nx/4), (Ny/4),  5 , 5);

    // std::vector<Point> diagLine = createDiagLine(Nx/2, Ny/2, Nx/2);
    // for (Point point : diagLine){
    //     std::cout << "(" + std::to_string(point.x) + "," + std::to_string(point.y) + ")" << std::endl;
    // }

    FluidSim sim(Nx, Ny, dx);
    sim.setGravity(0);
    setMultState(sim, rectangleOne);
    setMultState(sim, rectangleTwo);
    // setMultState(sim, diagLine);

    int scale = 3; // pixel size per cell
    sf::RenderWindow window(sf::VideoMode(Nx*scale, Ny*scale), "2D Fluid Simulation (Projection)");
    window.setFramerateLimit(60);

    sf::Image img;
    img.create(Nx, Ny, sf::Color::Black);
    sf::Texture tex;
    sf::Sprite spr;
    spr.setColor(sf::Color::Yellow);
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
                if (s == 0){
                    img.setPixel(i, Ny-1-j, sf::Color::Red);
                } else {
                    unsigned char c = (unsigned char)std::round(sim.clamp(s,0.0,1.0) * 255.0);
                    img.setPixel(i, Ny-1-j, sf::Color(c,c,c));
                }
                
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


std::vector<Point> createRectangle(double startX, double startY, int width, int height){
    std::vector<Point> points(width*height, Point(0,0));
    for (int j = 0; j < height; j++){
        for (int i = 0; i < width; i++){
            points[i + width*j].x = i + startX;
            points[i + width*j].y = j + startY;
        }
    }

    return points;
}


std::vector<Point> createDiagLine(double startX, double startY, int length){
    std::vector<Point> points(length, Point(0,0));
    for (int i = 0; i < length; i++){
        points[i].x = i + startX;
        points[i].y = i + startY;
    }

    return points;
}


void setMultState(FluidSim &_sim, std::vector<Point>_points){
    for (Point point : _points){
        _sim.setCellState(point.x, point.y, false);
    }
}
