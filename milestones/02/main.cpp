#include "verlet.h"
#include <iostream>

int main(int argc, char *argv[]) {
    const int frames = 30;
    const double timestep = 1;
    double x = 0;
    double y = 0;
    double z = 0;
    double vx = 0;
    double vy = 0;
    double vz = 0;
    double fx = 2;
    double fy = 2;
    double fz = 2;
    for (int i = 0; i < frames; i++) {
        verlet_step1(x,y,z,vx,vy,vz,fx,fy,fz,timestep);
        verlet_step2(x,y,z,vx,vy,vz,timestep);
        std::cout << "position: " << x << ',' << y << ',' << z << ' '
        << "velocity: " << vx << ',' << vy << ',' << vz << std::endl;
    }
}