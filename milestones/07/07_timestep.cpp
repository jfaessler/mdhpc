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
    constexpr int eq_steps_thermostat = 4000;
    constexpr int eq_steps_release = 4000;
    constexpr int eq_steps = eq_steps_thermostat + eq_steps_release;
    constexpr double eq_timestep = 15;
    constexpr double eq_temp = 500.0;
    constexpr double eq_relax = 10000.0;
    int steps = static_cast<int>(5000. / timestep);
    int snapshot_interval = steps > 100 ? steps / 100 : 1;
    std::ofstream traj("traj.xyz");
    std::ofstream ts(std::to_string(timestep) + ".timestep");

    auto [names, init_positions]{read_xyz("cluster_eq.xyz")};
    double gold_mass =
        20413.15887; // Gold in system's mass units where [m] = 0.009649 g/mol
    Atoms atoms(init_positions, gold_mass);
    atoms.k_b = 8.617333262e-5; // Boltzmann constant in eV/K

    const double cutoff = 5.0;
    NeighborList neighborList;
    neighborList.update(atoms, cutoff);

    write_xyz(traj, atoms);

    ts << "#PARAMS:Timestep=" << timestep << std::endl;
    ts << "Time,Total Energy" << std::endl;
    ts.precision(10);
    for (int i = 0; i < steps + eq_steps; ++i) {
        if (i < eq_steps) {
            verlet_step1(atoms.positions, atoms.velocities, atoms.forces,
                         eq_timestep, atoms.mass);
            neighborList.update(atoms, cutoff);
            ducastelle(atoms, neighborList, cutoff);
            verlet_step2(atoms.velocities, atoms.forces, eq_timestep, atoms.mass);
            if (i < eq_steps_thermostat)
                berendsen_thermostat(atoms, eq_temp, eq_timestep, eq_relax);
        } else {
            verlet_step1(atoms.positions, atoms.velocities, atoms.forces,
                         timestep, atoms.mass);
            neighborList.update(atoms, cutoff);
            double pot = ducastelle(atoms, neighborList, cutoff);
            verlet_step2(atoms.velocities, atoms.forces, timestep, atoms.mass);

            if (i % snapshot_interval == 0) {
                write_xyz(traj, atoms);
                auto kinetic = Eigen::pow(atoms.velocities, 2).sum() / 2;
                std::cout << (i - eq_steps + 1) * timestep << ","
                          << pot + kinetic << "," << pot << "," << kinetic
                          << "," << temperature(atoms) << std::endl;
                ts << (i - eq_steps + 1) * timestep << "," << pot + kinetic
                   << std::endl;
            }
        }
    }

    traj.close();
}