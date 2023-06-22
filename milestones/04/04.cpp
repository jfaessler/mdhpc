#include "lj_direct_summation.h"
#include "verlet.h"
#include "xyz.h"
#include <iostream>

int main(int argc, char *argv[]) {
    double timestep = 0.00001;
    int steps = static_cast<int>(1. / timestep / 10.);

    std::ofstream traj("traj.xyz");
    std::ofstream ts(std::to_string(timestep) + ".timestep");

    auto [names, init_positions,
          init_velocities]{read_xyz_with_velocities("lj54.xyz")};
    Atoms atoms(init_positions, init_velocities);

    write_xyz(traj, atoms);
    ts << 0 * timestep << ", " << pow(atoms.velocities.sum(), 2) / 2
       << std::endl;
    for (int i = 0; i < steps; ++i) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep);
        double pot = lj_direct_summation(atoms) / -1; // epsilon = 1
        verlet_step2(atoms.velocities, atoms.forces, timestep);
        write_xyz(traj, atoms);
        auto kinetic = Eigen::pow(atoms.velocities, 2).sum() / 2;
        ts << (i + 1) * timestep << ","
           << pot + kinetic << std::endl;
    }

    traj.close();
}