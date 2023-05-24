#include "lj_direct_summation.h"
#include <cmath>
double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    double pairwise_pot = 0;
    for (int i = 0; i < atoms.nb_atoms(); i++) {
        Eigen::Array3d force_sum;
        force_sum.setZero();
        for (int j = i + 1; j < atoms.nb_atoms(); j++) {
            Eigen::Array3d r = (atoms.positions(Eigen::all, i) -
                                atoms.positions(Eigen::all, j));
            //            force_sum += 24 * epsilon / r * (2 * sigma /
            //            Eigen::pow(r,12) - sigma / Eigen::pow(r, 6));
            double rn = r.matrix().norm();
            // Don't count self and avoid divide by 0
            if (rn > 1e-5) {
                force_sum +=
                    24 * epsilon / r *
                    (2 * Eigen::pow(sigma / r, 12) - Eigen::pow(sigma / r, 6));
                //            double r =
                //            (atoms.positions(Eigen::all,i)-atoms.positions(Eigen::all,
                //            j)).matrix().norm(); force_sum += 24 * epsilon / r *
                //            (2 * std::pow(sigma / r, 12) - std::pow(sigma / r,
                //            6));
                pairwise_pot +=
                    4 * epsilon *
                    (std::pow(sigma / rn, 12) - std::pow(sigma / rn, 6));
            }
        }
        atoms.forces(Eigen::all, i) = force_sum;
    }
    return pairwise_pot / 2; // Don't count self
}
