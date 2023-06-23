#include "lj_direct_summation.h"
#include "thermostat.h"
#include "verlet.h"
#include "xyz.h"
#include <iostream>

int main(int argc, char *argv[]) {
    const double timestep = 0.001;
    const int steps = 30000;
    const int snapshot_interval = 100;

    const double epsilon = 1.0;
    const double sigma = 1.0;

    const double target_temp = 0.;
    const double relaxation = 1.0;

    const int lattice_size = 5;
    const double spacing = 1. * sigma; // Lattice constant
    const int nb_atoms = 100;
    assert(lattice_size * lattice_size * lattice_size >= nb_atoms);

    Positions_t init_positions;
    init_positions.resize(3, nb_atoms);
    for (int i = 0; i < nb_atoms; i++) {
        double x = i % lattice_size * spacing;
        double y = i / lattice_size % lattice_size * spacing;
        double z = i / lattice_size / lattice_size % lattice_size * spacing;
        init_positions.col(i) << x, y, z;
    }

    Atoms atoms(init_positions);
//    atoms.velocities.setRandom();

    std::ofstream traj("traj.xyz");
    std::cout << "Time,Total Energy,Potential,Kinetic,Temperature" << std::endl;

        write_xyz(traj, atoms);
    for (int i = 0; i < steps; ++i) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep);
        double pot = lj_direct_summation(atoms, epsilon, sigma);
        verlet_step2(atoms.velocities, atoms.forces, timestep);
        berendsen_thermostat(atoms, target_temp, timestep, relaxation);

        if (i % snapshot_interval == 0) {
            write_xyz(traj, atoms);
            auto kinetic = Eigen::pow(atoms.velocities, 2).sum() / 2;
            std::cout << (i + 1) * timestep << "," << pot + kinetic << ","
                      << pot << "," << kinetic << "," << atoms_temperature(atoms) << std::endl;
        }
    }

    traj.close();
}