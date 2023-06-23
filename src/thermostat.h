#ifndef MD_CODE_THERMOSTAT_H
#define MD_CODE_THERMOSTAT_H

#include "atoms.h"

void berendsen_thermostat(Atoms &atoms, double temperature, double timestep,
                          double relaxation_time);

double atoms_temperature(const Atoms &atoms);

#endif // MD_CODE_THERMOSTAT_H
