#include "thermostat.h"
#include "verlet.h"
#include "xyz.h"
#include "neighbors.h"
#include "ducastelle.h"
#include <iostream>



int main(int argc, char *argv[]) {
    constexpr double timestep = 5.;
    constexpr int steps = 10000;

    constexpr int snapshot_interval = steps / 100; // 100 total frames
    // TODO set according to time accumulated

    std::ofstream traj("traj.xyz");
    std::ofstream ts(std::to_string(timestep) + ".timestep");

    auto [names, init_positions]{read_xyz("cluster_923.xyz")};
    Atoms atoms(init_positions);
    atoms.mass = 20413.15887; // Mass system's mass units where g/mol = 0.009649
    atoms.k_b = 8.617333262e-5; // Boltzmann constant in eV/K

    const double cutoff = 10.0; // Default cutoff from ducastelle
    NeighborList neighborList;
    neighborList.update(atoms, cutoff);

    write_xyz(traj, atoms);

//    const double target_temp = 0.01;
//    const double relaxation = .1;

    std::cout << "Time,Total Energy,Potential,Kinetic,Temperature" << std::endl;
    std::cout.precision(10);
    ts.precision(10);
    for (int i = 0; i < steps; ++i) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep);
        neighborList.update(atoms, cutoff);
        double pot = ducastelle(atoms, neighborList, cutoff);
        verlet_step2(atoms.velocities, atoms.forces, timestep);
//        berendsen_thermostat(atoms, target_temp, timestep, relaxation);

        if (i % snapshot_interval == 0) {
            write_xyz(traj, atoms);
            auto kinetic = kinetic_energy(atoms);
            std::cout << (i + 1) * timestep << "," << pot + kinetic << ","
                      << pot << "," << kinetic << "," << temperature(atoms) << std::endl;
            ts << (i + 1) * timestep << ","
               << pot + kinetic << std::endl;
        }
    }

    traj.close();
}