#ifndef MD_CODE_ATOMS_H
#define MD_CODE_ATOMS_H

#include "types.h"

class Atoms {
  public:
    Positions_t positions;
    Velocities_t velocities;
    Forces_t forces;

    // Default to LB units
    double mass = 1;
    double k_b = 1;

    explicit Atoms(int nb_atoms);

    explicit Atoms(const Positions_t &p);

    Atoms(const Positions_t &p, const Velocities_t &v);

    [[nodiscard]] size_t nb_atoms() const;
};

#endif // MD_CODE_ATOMS_H
