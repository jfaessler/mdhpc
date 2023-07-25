#include "thermostat.h"
#include "verlet.h"
#include "xyz.h"
#include "neighbors.h"
#include "ducastelle.h"
#include <iostream>



int main(int argc, char *argv[]) {
    constexpr double timestep = 15.0;
    constexpr int steps = 52001;
    constexpr double eq_temp = 500.0;
    constexpr int eq_steps = 2000;
    constexpr int eq_relax = 10.0; // Set to very low value to equilibrate the system
    constexpr int tau_relax = 50;
    constexpr double delta_q = 4.0; // Increase in temperature per heating cycle
    // TODO heating cycles should not use berendensen thermostat!
    // Rescale for a certain change in kinetic, then let it wait for 1000
    // And take the average for the last like 500 to let it equilibriate


    constexpr int snapshot_interval = steps / 100; // 100 total frames
    // TODO set according to time accumulated

    std::ofstream traj("traj.xyz");
    std::ofstream ts(std::to_string(timestep) + ".timestep");

    auto [names, init_positions]{read_xyz("cluster_923.xyz")};
    double gold_mass = 20413.15887; // Gold in system's mass units where g/mol = 0.009649
    Atoms atoms(init_positions, gold_mass);
    atoms.k_b = 8.617333262e-5; // Boltzmann constant in eV/K

    const double cutoff = 5.0; // Default cutoff from ducastelle
    NeighborList neighborList;
    neighborList.update(atoms, cutoff);

    std::cout << "Step,Time,Total Energy,Potential,Kinetic,Temperature" << std::endl;
    std::cout.precision(10);
    ts.precision(10);

    std::vector<double> average_temp;
    double temp_totals = 0;
    for (int i = 0; i < steps; ++i) {
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, atoms.mass);
        neighborList.update(atoms, cutoff);
        double pot = ducastelle(atoms, neighborList, cutoff);
        verlet_step2(atoms.velocities, atoms.forces, timestep, atoms.mass);
        if (i < eq_steps)
            berendsen_thermostat(atoms, eq_temp, timestep, eq_relax);
        else {
            int cycle = (i - eq_steps) % (tau_relax * 2);
            if (cycle == 0) {
                double current_temp;
                // Deposit energy into the system, take average
                if (i != eq_steps) {
                    average_temp.emplace_back(temp_totals / tau_relax);
                    current_temp = average_temp.back();
                    temp_totals = 0;
                } else {
                    current_temp = eq_temp;
                }

                // Deposit kinetic energy into the system via velocity rescaling
                atoms.velocities *= sqrt((current_temp + delta_q) / current_temp);
            } else if (cycle >= tau_relax) {
                // Measure numerator of average temperature,
                // making sure we have relaxed to let the system settle first
                temp_totals += temperature(atoms);
            }
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

    std::cout << "Average temperature" << std::endl;
    for (const auto& temp : average_temp) {
        std::cout << temp << std::endl;
    }

    traj.close();
}