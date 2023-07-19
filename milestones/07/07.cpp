#include "thermostat.h"
#include "verlet.h"
#include "xyz.h"
#include "neighbors.h"
#include "ducastelle.h"
#include <iostream>



int main(int argc, char *argv[]) {
    constexpr double timestep = 1.0;
    constexpr int steps = 8000;

    constexpr int snapshot_interval = steps / 100; // 100 total frames
    // TODO set according to time accumulated

    std::ofstream traj("traj.xyz");
    std::ofstream ts(std::to_string(timestep) + ".timestep");

    auto [names, init_positions]{read_xyz("cluster_923.xyz")};
    Atoms atoms(init_positions);
    atoms.mass = 20413.15887; // Mass system's mass units where g/mol = 0.009649
    atoms.k_b = 8.617333262e-5; // Boltzmann constant in eV/K

    const double cutoff = 5.0; // Default cutoff from ducastelle
    NeighborList neighborList;
    neighborList.update(atoms, cutoff);

    double target_temp = 100.0;
    const double relaxation = 1000;

    const int s = 800; // Scale for setting new temp value

    std::cout << "Step,Time,Total Energy,Potential,Kinetic,Temperature" << std::endl;
    std::cout.precision(10);
    ts.precision(10);
    for (int i = 0; i < steps; ++i) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, atoms.mass);
        neighborList.update(atoms, cutoff);
        double pot = ducastelle(atoms, neighborList, cutoff);
        verlet_step2(atoms.velocities, atoms.forces, timestep, atoms.mass);
        if ((i / s) % 2 == 0)
            berendsen_thermostat(atoms, target_temp, timestep, relaxation);

        if (i == s) {
            target_temp = 500;
        }
        if (i == 3 * s) {
            target_temp = 600;
        }

        if (i == 5 * s) {
            target_temp = 650;
        }

        if (i == 7 * s) {
            target_temp = 700;
        }

        if (i == 9 * s) {
            target_temp = 725;
        }

        if (i == 11 * s) {
            target_temp = 750;
        }

        if (i == 12 * s) {
            target_temp = 800;
        }

        if (i == 14 * s) {
            target_temp = 900;
        }
        if (i % snapshot_interval == 0) {
            write_xyz(traj, atoms);
            auto kinetic = kinetic_energy(atoms);
            std::cout << i << "," << (i) * timestep << "," << pot + kinetic << ","
                      << pot << "," << kinetic << "," << temperature(atoms) << std::endl;
            ts << (i + 1) * timestep << ","
               << pot + kinetic << std::endl;
        }
    }

    traj.close();
}