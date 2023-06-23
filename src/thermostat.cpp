#include "thermostat.h"
void berendsen_thermostat(Atoms &atoms, double target, double timestep,
                          double relaxation_time) {
    // TODO assume mass and k is 1?
    double current = atoms_temperature(atoms);
    double lambda = std::sqrt(1+(target / current - 1) * timestep / relaxation_time);
    atoms.velocities *= lambda;
}

double atoms_temperature(const Atoms &atoms) {
    return atoms.velocities.square().sum() / 3 / static_cast<double>(atoms.nb_atoms());
}
