#include "lj_direct_summation.h"
#include "verlet.h"
#include "xyz.h"
#include <iostream>

int main(int argc, char *argv[]) {
    double timestep = 0.010;
    int steps = 100;

    std::ofstream traj("traj.xyz");
    std::ofstream ts(std::to_string(timestep) + ".timestep");

    auto [names, positions, velocities]{read_xyz_with_velocities("lj54.xyz")};
    Atoms atoms(positions, velocities);

    write_xyz(traj, atoms);
    ts << 0 * timestep << ", " << pow(atoms.velocities.sum(), 2) / 2 << std::endl; // mass = 1
    for (int i = 0; i < steps; ++i) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep);
        lj_direct_summation(atoms);
        verlet_step2(atoms.velocities, atoms.forces, timestep);
        write_xyz(traj, atoms);
        ts << (i + 1) * timestep << "," << pow(atoms.velocities.sum(), 2) / 2 << std::endl; // mass = 1
    }

    traj.close();
}