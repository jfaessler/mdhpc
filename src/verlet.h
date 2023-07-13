#ifndef __VERLET_H
#define __VERLET_H

#include "types.h"

void verlet_step1(double &x, double &y, double &z, double &vx, double &vy, double &vz,
                  double fx, double fy, double fz, double timestep);
void verlet_step2(double &vx, double &vy, double &vz, double fx, double fy, double fz,
                  double timestep);

void verlet_step1(Positions_t &positions, Velocities_t &velocities,
                  const Velocities_t &forces, Scalar_t timestep, Scalar_t mass = 1.);
void verlet_step2(Velocities_t &velocities, const Forces_t &forces, Scalar_t timestep, Scalar_t mass = 1.);

#endif // __VERLET_H
