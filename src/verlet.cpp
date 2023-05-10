//
// Created by joshua on 4/26/23.
//

#include "verlet.h"
void verlet_step1(double &x, double &y, double &z, double &vx,
                          double &vy, double &vz, double fx, double fy,
                          double fz, double timestep) {
    vx += .5 * fx * timestep;
    vy += .5 * fy * timestep;
    vz += .5 * fz * timestep;
    x += .5 * vx * timestep;
    y += .5 * vy * timestep;
    z += .5 * vz * timestep;
}
void verlet_step2(double &vx, double &vy, double &vz, double fx,
                          double fy, double fz, double timestep) {
    vx += .5 * fx * timestep;
    vy += .5 * fy * timestep;
    vz += .5 * fz * timestep;
}
