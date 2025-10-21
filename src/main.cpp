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

std::vector<double> getSpectrumColor(double val);


int main(int argc, char* argv[]) {

    if (argc == 2){
        testMatrixASetup(4, 4, 1.0, 1.0, 1.0);
    } else {

        int Nx = 100, Ny = 100;
        double dx = 1.0 / Nx;

        std::vector<Point> rectangleOne = createRectangle((Nx/2) - 2, (Ny/2) - 2,  2 , 2);
        for (Point point : rectangleOne){
            std::cout << "(" + std::to_string(point.x) + "," + std::to_string(point.y) + ")" << std::endl;
        }

        std::vector<Point> rectangleTwo = createRectangle((Nx/2) - 3, (Ny/2) -4,  2 , 2);

        std::vector<Point> rectangleThree = createRectangle(74, 26, 1, (Nx/2)-1);

        std::vector<Point> diagLineOne = createDiagLinePosY(Nx/2, Ny/2, Nx/4);
        for (Point point : diagLineOne){
            std::cout << "(" + std::to_string(point.x) + "," + std::to_string(point.y) + ")" << std::endl;
        }

        std::vector<Point> diagLineTwo = createDiagLineNegY(Nx/2, Ny/2, Nx/4);
        for (Point point : diagLineTwo){
            std::cout << "(" + std::to_string(point.x) + "," + std::to_string(point.y) + ")" << std::endl;
        }

        FluidSim sim(Nx, Ny, dx);
        sim.setGravity(0);
        // setMultState(sim, rectangleOne);
        // setMultState(sim, rectangleTwo);
        setMultState(sim, rectangleThree);
        setMultState(sim, diagLineOne);
        setMultState(sim, diagLineTwo);

        int scale = 7; // pixel size per cell
        sf::RenderWindow window(sf::VideoMode(Nx*scale, Ny*scale), "2D Fluid Simulation (Projection)");
        window.setFramerateLimit(60);

        sf::Image img;
        img.create(Nx, Ny, sf::Color::Black);
        sf::Texture tex;
        sf::Sprite spr;
        spr.setColor(sf::Color(255, 228, 10));
        tex.create(Nx, Ny);
        spr.setScale((float)scale, (float)scale);

        while (window.isOpen()) {

            std::vector<double> color;

            sf::Event event;
            while (window.pollEvent(event)) {
                if (event.type == sf::Event::Closed)
                    window.close();
            }

            sim.step();

            for (int j = 0; j < Ny; ++j) {
                for (int i = 0; i < Nx; ++i) {
                    double s = sim.dye[sim.idxP(i,j)];

                    color = getSpectrumColor(sim.clamp(s, 0.0, 1.0));

                    img.setPixel(i, Ny-1-j, sf::Color(color[0],color[1],color[2]));

                }
            }
            tex.update(img);
            spr.setTexture(tex);

            window.clear();
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



std::vector<double> getSpectrumColor(double val){
    double r, g, b;
    // Scale value to 0 - 6
    double segment = val * 6;

    // Segment 1: Red -> Yellow (increase Green)
    if (segment > 0 && segment < 1){
        r = 255;
        g = segment * 255;
        b = 0;

    // Segment 2: Yellow -> Green (decrease Red)
    } else if (segment >= 1 && segment < 2){
        r = (2 - segment) * 255;
        g = 255;
        b = 0;

    // Segment 3: Green -> Cyan (increase Blue)
    } else if (segment >= 2 && segment < 3){
        r = 0;
        g = 255;
        b = (segment - 2) * 255;

    // Segment 4: Cyan -> Blue (decrease Green)
    } else if (segment >= 3 && segment < 4){
        r = 0;
        g = (4 - segment) * 255;
        b = 255;

    // Segment 5: Blue -> Magenta (increase Red) - For full spectrum loop
    } else if (segment >= 4 && segment < 5){
        r = (segment - 4) * 255;
        g = 0;
        b = 255;

    // Segment 6: Magenta -> Red (decrease Blue) - For full spectrum loop
    } else if (segment >= 5 && segment <= 6){
        r = 255;
        g = 0;
        b = (6 - segment) * 255;

    // Handle edge cases for t = 1.0 to ensure it's red
    } else if (segment >= 1) {
        r = 255;
        g = 0;
        b = 0;
    } else {
        r = 0;
        g = 0;
        b = 0;
    }

    return {floor(r), floor(g), floor(b)};

}