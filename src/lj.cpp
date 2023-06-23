#include "lj_direct_summation.h"
#include "neighbors.h"
#include <cmath>
double lj(Atoms &atoms, NeighborList &neighbor_list, double cutoff,
          double epsilon, double sigma) {
    double pairwise_pot = 0;
    // Calculate the amount to subtract from potential due to the cutoff
    double cutoff_energy =
        4 * epsilon *
        (std::pow(sigma / cutoff, 12) - std::pow(sigma / cutoff, 6));
    atoms.forces.setZero();
    neighbor_list.update(atoms, cutoff);
    for (auto [i, j] : neighbor_list) {
        if (i < j) {
            Eigen::Array3d r = (atoms.positions(Eigen::all, i) -
                                atoms.positions(Eigen::all, j));
            double rn = r.matrix().norm();
            // Don't count self and avoid divide by 0
            Eigen::Array3d sdr = sigma / r;
            Eigen::Array3d edr = epsilon / r;
            auto force =
                r / rn * 24 * (epsilon / rn) *
                (2 * std::pow(sigma / rn, 12) - std::pow(sigma / rn, 6));
            pairwise_pot +=
                4 * epsilon *
                    (std::pow(sigma / rn, 12) - std::pow(sigma / rn, 6)) -
                cutoff_energy;
            atoms.forces(Eigen::all, i) += force;
            atoms.forces(Eigen::all, j) -= force;
        }
    }
    return pairwise_pot;
}
