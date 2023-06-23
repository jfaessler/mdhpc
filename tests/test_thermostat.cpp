#include "thermostat.h"
#include <gtest/gtest.h>

TEST(ThermostatTest, Zero) {
    Atoms atoms(10);
    atoms.velocities.setOnes();
    Velocities_t starting_velocities(atoms.velocities);
    berendsen_thermostat(atoms, 1., 0.01, 0.1);
    for (int i = 0; i < atoms.velocities.size(); i++) {
        EXPECT_EQ(atoms.velocities.reshaped()[i],
                  starting_velocities.reshaped()[i]);
    }
}

TEST(ThermostatTest, Equilibruim) {
    Atoms atoms(10);
    atoms.velocities.setRandom();
    for (int i = 0; i < 100; i++) {
        berendsen_thermostat(atoms, 3, 0.01, 0.1);
    }
    EXPECT_NEAR(atoms_temperature(atoms), 3, 0.001);
}