#include "lj_direct_summation.h"
#include "verlet.h"
#include "xyz.h"
#include <iostream>

int main(int argc, char *argv[]) {
    double timestep;
    if (argc > 1) {
        timestep = std::stod(argv[1]);
    } else {
        timestep = 0.001;
    }
    int steps = static_cast<int>(1. / timestep / 1.);

    std::ofstream traj("traj.xyz");
    std::ofstream ts(std::to_string(timestep) + ".timestep");
    ts << "#PARAMS:Timestep=" << timestep << std::endl;
    ts << "Time,Total Energy" << std::endl;

    auto [names, init_positions,
          init_velocities]{read_xyz_with_velocities("lj54.xyz")};
    Atoms atoms(init_positions, init_velocities);

    write_xyz(traj, atoms);
    for (int i = 0; i < steps; ++i) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep);
        double pot = lj_direct_summation(atoms);
        verlet_step2(atoms.velocities, atoms.forces, timestep);
        write_xyz(traj, atoms);
        auto kinetic = Eigen::pow(atoms.velocities, 2).sum() / 2;
        ts << (i + 1) * timestep << ","
           << pot + kinetic << std::endl;
    }

    traj.close();
}