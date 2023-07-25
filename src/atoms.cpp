#include "atoms.h"
Atoms::Atoms(const Positions_t &p)
    : positions{p},
      velocities{3, p.cols()},
      forces{3, p.cols()},
      masses{p.cols()} {
    velocities.setZero();
    forces.setZero();
    masses.setOnes();
}
Atoms::Atoms(const Positions_t &p, const Velocities_t &v)
    : positions{p},
      velocities{v},
      forces{3, p.cols()},
      masses{p.cols()} {
    assert(p.cols() == v.cols());
    forces.setZero();
    masses.setOnes();
}
Eigen::Index Atoms::nb_atoms() const {
    return positions.cols();
}
Atoms::Atoms(int nb_atoms)
    : positions{3, nb_atoms},
      velocities{3, nb_atoms},
      forces{3, nb_atoms},
      masses{nb_atoms} {
    positions.setZero();
    velocities.setZero();
    forces.setZero();
    masses.setOnes();
}

Atoms::Atoms(const Positions_t &p, double atom_mass)
    : positions{p},
      velocities{3, p.cols()},
      forces{3, p.cols()},
      masses{p.cols()},
      mass(atom_mass) {
    positions.setZero();
    velocities.setZero();
    forces.setZero();
    masses.setConstant(mass);
}

void Atoms::resize(int i) {
    positions.conservativeResize(Eigen::NoChange, i);
    velocities.conservativeResize(Eigen::NoChange, i);
    forces.conservativeResize(Eigen::NoChange, i);
    masses.conservativeResize(i);
}
