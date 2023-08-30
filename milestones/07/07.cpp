#include "ducastelle.h"
#include "neighbors.h"
#include "thermostat.h"
#include "verlet.h"
#include "xyz.h"
#include <iostream>

int main(int argc, char *argv[]) {
    constexpr double timestep = 15.0;
    constexpr int steps = 58001;
    constexpr double eq_temp = 500.0;
    constexpr int eq_steps_thermostat = 4000;
    // Let free after heating to finish equilibrating
    constexpr int eq_steps_release = 4000;
    constexpr int eq_steps = eq_steps_thermostat + eq_steps_release;
    constexpr double eq_relax = 10000.0;
    constexpr int tau_relax = 50;
    // Increase in energy per heating cycle
    // Is actually implemented as velocity rescaling by a temperature factor,
    // to allow for different clusters to take roughly the same amount of
    // timesteps to heat up
    constexpr double delta_q = 4.0;

    constexpr int snapshot_interval = steps / 100; // 100 total frames

    std::ofstream traj("traj.xyz");
    std::string in_filename;
    std::string data_filename;
    if (argc >= 3) {
        in_filename = argv[1];
        data_filename = argv[2];
    } else {
        in_filename = "cluster_923.xyz";
        data_filename = "average.csv";
    }
    std::ofstream data(data_filename);

    auto [names, init_positions]{read_xyz(in_filename)};
    double gold_mass =
        20413.15887; // Gold in system's mass units where [m] = 0.009649 g/mol
    Atoms atoms(init_positions, gold_mass);
    atoms.k_b = 8.617333262e-5; // Boltzmann constant in eV/K

    const double cutoff = 5.0;
    NeighborList neighborList;
    neighborList.update(atoms, cutoff);

    std::cout << "Step,Time,Total Energy,Potential,Kinetic,Temperature"
              << std::endl;
    std::cout.precision(10);
    //    ts.precision(10);

    // Data structures for storing averages
    std::vector<double> average_temp;
    std::vector<double> average_kinetic;
    std::vector<double> average_pot;

    // Accumulators for averages
    double temp_totals = 0;
    double kinetic_totals = 0;
    double pot_totals = 0;

    for (int i = 0; i < steps; ++i) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep,
                     atoms.mass);
        neighborList.update(atoms, cutoff);
        double pot = ducastelle(atoms, neighborList, cutoff);
        verlet_step2(atoms.velocities, atoms.forces, timestep, atoms.mass);
        if (i < eq_steps) {
            if (i < eq_steps_thermostat)
                berendsen_thermostat(atoms, eq_temp, timestep, eq_relax);
        } else {
            int cycle = (i - eq_steps) % (tau_relax * 2);
            if (cycle == 0) {
                double current_temp;
                // Deposit energy into the system, take averages
                if (i != eq_steps) {
                    average_temp.emplace_back(temp_totals / tau_relax);
                    current_temp = average_temp.back();
                    average_kinetic.emplace_back(kinetic_totals / tau_relax);
                    average_pot.emplace_back(pot_totals / tau_relax);

                    temp_totals = 0;
                    kinetic_totals = 0;
                    pot_totals = 0;
                } else {
                    current_temp = eq_temp;
                }

                // Deposit kinetic energy into the system via velocity rescaling
                atoms.velocities *=
                    sqrt((current_temp + delta_q) / current_temp);
            } else if (cycle >= tau_relax) {
                // Measure numerator of average temperature,
                // making sure we have relaxed to let the system settle first
                temp_totals += temperature(atoms);
                kinetic_totals += kinetic_energy(atoms);
                pot_totals += pot;
            }
        }

        if (i % snapshot_interval == 0) {
            write_xyz(traj, atoms);
            auto kinetic = kinetic_energy(atoms);
            std::cout << i << "," << (i)*timestep << "," << pot + kinetic << ","
                      << pot << "," << kinetic << "," << temperature(atoms)
                      << std::endl;
        }
    }

    data << "#PARAMS:"
         << "size=" << atoms.nb_atoms() << ",delta_q=" << delta_q
         << ",eq_steps=" << eq_steps << std::endl;
    data << "Step,Cycle,Time,Average Temperature,Kinetic,Potential"
         << std::endl;
    for (size_t i = 0; i < average_temp.size(); i++) {
        size_t step = i * tau_relax * 2 + eq_steps; // full cycle is 2tau
        data << step << ',' << i << ',' << static_cast<double>(i) * timestep
             << ',' << average_temp[i] << ',' << average_kinetic[i] << ','
             << average_pot[i] << std::endl;
    }

    traj.close();
    data.close();
}