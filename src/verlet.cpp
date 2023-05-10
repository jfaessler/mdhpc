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
    x += vx * timestep;
    y += vy * timestep;
    z += vz * timestep;
}
void verlet_step2(double &vx, double &vy, double &vz, double fx,
                          double fy, double fz, double timestep) {
    vx += .5 * fx * timestep;
    vy += .5 * fy * timestep;
    vz += .5 * fz * timestep;
}
void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Velocities_t &forces, Scalar_t timestep) {
    velocities += .5 * forces * timestep;
    positions += velocities * timestep;
}
void verlet_step2(Velocities_t &velocities, const Forces_t &forces,
                  Scalar_t timestep) {
    velocities += .5 * forces * timestep;
}
