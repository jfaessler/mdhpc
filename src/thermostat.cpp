#include "thermostat.h"
void berendsen_thermostat(Atoms &atoms, double target, double timestep,
                          double relaxation_time) {
    double current = temperature(atoms);
    double lambda = std::sqrt(1+(target / current - 1) * timestep / relaxation_time);
    atoms.velocities *= lambda;
}

double temperature(const Atoms &atoms) {
    return atoms.mass / atoms.k_b * atoms.velocities.square().sum() / 3. / static_cast<double>(atoms.nb_atoms());
}

double kinetic_energy(const Atoms &atoms) {
    return atoms.mass * atoms.velocities.square().sum() / 2.;
}