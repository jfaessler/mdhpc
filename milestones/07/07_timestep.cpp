#include "thermostat.h"
#include "verlet.h"
#include "xyz.h"
#include "neighbors.h"
#include "ducastelle.h"
#include <iostream>

int main(int argc, char *argv[]) {
    double timestep;
    if (argc > 1) {
        timestep = std::stod(argv[1]);
    } else {
        timestep = 1.;
    }
    int steps = static_cast<int>(500. / timestep);
    std::ofstream traj("traj.xyz");
    std::ofstream ts(std::to_string(timestep) + ".timestep");

    auto [names, init_positions]{read_xyz("cluster_923.xyz")};
    double gold_mass =
        20413.15887; // Gold in system's mass units where [m] = 0.009649 g/mol
    Atoms atoms(init_positions, gold_mass);
    atoms.k_b = 8.617333262e-5; // Boltzmann constant in eV/K

    const double cutoff = 5.0;
    NeighborList neighborList;
    neighborList.update(atoms, cutoff);

    write_xyz(traj, atoms);

    std::cout << "Time,Total Energy,Potential,Kinetic,Temperature" << std::endl;
    std::cout.precision(10);
    ts << "Time,Total Energy" << std::endl;
    ts.precision(10);
    for (int i = 0; i < steps; ++i) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, atoms.mass);
        neighborList.update(atoms, cutoff);
        double pot = ducastelle(atoms, neighborList, cutoff);
        verlet_step2(atoms.velocities, atoms.forces, timestep, atoms.mass);

        write_xyz(traj, atoms);
        auto kinetic = Eigen::pow(atoms.velocities, 2).sum() / 2;
        std::cout << (i + 1) * timestep << "," << pot + kinetic << ","
                  << pot << "," << kinetic << "," << temperature(atoms) << std::endl;
        ts << (i + 1) * timestep << ","
           << pot + kinetic << std::endl;
    }

    traj.close();
}