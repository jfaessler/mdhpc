#include "verlet.h"
#include <gtest/gtest.h>

TEST(VerletTest, Single) {
    const int steps = 10;
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
    for (int i = 0; i < steps; ++i) {
        verlet_step1(x,y,z,vx,vy,vz,fx,fy,fz,timestep);
//        ... compute forces here ...
        verlet_step2(x,y,z,vx,vy,vz,timestep);
    }
}

TEST(VerletTest, Multiple) {
    const int steps = 10;
    const double timestep = 1;
    int nb_atoms = 10;

    Positions_t positions(3, nb_atoms);
    for (int i = 0; i < nb_atoms; i++) {
        positions(1,i) = i;
    }
    Positions_t r0 (positions);

    Velocities_t velocities(3, nb_atoms);
    velocities.setOnes();
    Velocities_t v0 (velocities);

    Forces_t forces(3, nb_atoms);
    forces(1, Eigen::all).setConstant(1);

    for (int i = 1; i < steps; i++) {
        verlet_step1(positions, velocities, forces, timestep);
        verlet_step2(velocities, forces, timestep);
        int j = 0;
        double t = i * timestep;
        for (const auto& item : positions.reshaped()) {
            EXPECT_NEAR(item, .5 * forces.reshaped()[j] * t * t + v0.reshaped()[j] * t + r0.reshaped()[j], 1e-6);
            j++;
        }
        j = 0;
        for (const auto& item : velocities.reshaped()) {
            EXPECT_NEAR(item, forces.reshaped()[j] * t +v0.reshaped()[j], 1e-6);
            j++;
        }
    }
}