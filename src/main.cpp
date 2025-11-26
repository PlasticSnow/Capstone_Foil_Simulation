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
#include "Tests.hpp"


struct Point {
    double x;
    double y;

    Point(double _x, double _y){
        x = _x;
        y = _y;
    }

};


std::vector<Point> createRectangle(double startX, double startY, int with, int height);

std::vector<Point> createDiagLinePosY(double startX, double startY, int length);


std::vector<Point> createDiagLineNegY(double startX, double startY, int length);

void setMultState(FluidSim &_sim, std::vector<Point>_points);

std::vector<double> getSpectrumColor(double val, bool fluid);


int main(int argc, char* argv[]) {

    if (argc == 2){
        testMatrixASetup(4, 4, 1.0, 1.0, 1.0);
    } else {

        int Nx = 100, Ny = 100;
        double dx = .1;

        std::vector<Point> rectangleOne = createRectangle((Nx/2) - 2, (Ny/2) - 2,  2 , 2);
        for (Point point : rectangleOne){
            std::cout << "(" + std::to_string(point.x) + "," + std::to_string(point.y) + ")" << std::endl;
        }

        std::vector<Point> rectangleTwo = createRectangle((Nx/2) - 3, (Ny/2) -4,  2 , 2);

        std::vector<Point> rectangleThree = createRectangle(74, 26, 1, (Nx/2)-1);

        std::vector<Point> diagLineOne = createDiagLinePosY(Nx/4, Ny/2, Nx/10);
        for (Point point : diagLineOne){
            std::cout << "(" + std::to_string(point.x) + "," + std::to_string(point.y) + ")" << std::endl;
        }

        std::vector<Point> diagLineTwo = createDiagLineNegY(Nx/4, Ny/2, Nx/5);
        for (Point point : diagLineTwo){
            std::cout << "(" + std::to_string(point.x) + "," + std::to_string(point.y) + ")" << std::endl;
        }

        FluidSim sim(Nx, Ny, dx);
        sim.setGravity(0);
        setMultState(sim, rectangleOne);
        setMultState(sim, rectangleTwo);
        // setMultState(sim, rectangleThree);
        // setMultState(sim, diagLineOne);
        // setMultState(sim, diagLineTwo);

        int scale = 4; // pixel size per cell
        sf::RenderWindow window(sf::VideoMode(Nx*scale, Ny*scale), "2D Fluid Simulation (Projection)");
        window.setFramerateLimit(60);

        sf::Image img;
        img.create(Nx, Ny, sf::Color::Black);
        sf::Texture tex;
        sf::Sprite spr;
        spr.setColor(sf::Color::White);
        tex.create(Nx, Ny);
        spr.setScale((float)scale, (float)scale);
        bool continueSim = true;

        double cellIndex = 0;
        while (window.isOpen()) {

            std::vector<double> color;

            sf::Event event;
            while (window.pollEvent(event)) {
                if (event.type == sf::Event::Closed)
                    window.close();
            }


            sf::Vector2i localPosition = sf::Mouse::getPosition(window);
            int x_ = localPosition.x /scale;
            int y_ = (Ny*scale - localPosition.y) /scale;
            if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Left)){
                continueSim = false;
                std::cout << "x: " << x_ << " y: " << y_ << "\n";
                sim.cellState[sim.idxP(x_,y_)] = false;
            } else if (sf::Mouse::isButtonPressed(sf::Mouse::Button::Right)) {
                continueSim = false;
                sim.cellState[sim.idxP(x_,y_)] = true;
            } else {
                continueSim = true;
            }


            if (continueSim){
                sim.step();
            }
            


            
            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    cellIndex = sim.idxP(i,j);
                    double s = sim.dye[cellIndex];

                    color = getSpectrumColor(sim.clamp(s, 0.0, 1.0), sim.cellState[cellIndex]);
                    img.setPixel(i, Ny-1-j, sf::Color(color[0],color[1],color[2]));
                    


                }
            }
            tex.update(img);
            spr.setTexture(tex);

            window.clear(sf::Color::Black);
            window.draw(spr);
            window.display();
        }

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


std::vector<Point> createDiagLinePosY(double startX, double startY, int length){
    std::vector<Point> points(length, Point(0,0));
    for (int i = 0; i < length; i++){
        points[i].x = i + startX;
        points[i].y = i + startY;
    }

    return points;
}


std::vector<Point> createDiagLineNegY(double startX, double startY, int length){
    std::vector<Point> points(length, Point(0,0));
    for (int i = 0; i < length; i++){
        points[i].x = i + startX;
        points[i].y = startY - i;
    }

    return points;
}


void setMultState(FluidSim &_sim, std::vector<Point>_points){
    for (Point point : _points){
        _sim.setCellState(point.x, point.y, false);
    }
}



std::vector<double> getSpectrumColor(double val, bool fluid) {

    if (!fluid){
        return {255.0,255.0,255.0};
    }

    // Clamp input to [0, 1]
    val = std::clamp(val, 0.0, 1.0);

    // Map value to hue: 0 -> red, 1 -> violet (can tweak range)
    double hue = (1.0 - val) * 240.0; // 0° = red, 120° = green, 240° = blue

    double s = 1.0; // full saturation
    double v = 1.0; // full brightness

    double c = v * s;
    double x = c * (1 - fabs(fmod(hue / 60.0, 2) - 1));
    double m = v - c;

    double r, g, b;
    if (hue < 60) { r = c; g = x; b = 0; }
    else if (hue < 120) { r = x; g = c; b = 0; }
    else if (hue < 180) { r = 0; g = c; b = x; }
    else if (hue < 240) { r = 0; g = x; b = c; }
    else if (hue < 300) { r = x; g = 0; b = c; }
    else { r = c; g = 0; b = x; }

    return { (r + m) * 255.0, (g + m) * 255.0, (b + m) * 255.0 };
}
