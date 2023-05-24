#include "lj_direct_summation.h"
#include <cmath>
double lj_direct_summation(Atoms &atoms, double epsilon, double sigma) {
    double sum = 0;
    for (int i = 0; i < atoms.nb_atoms(); i++) {
        for (int j = i + 1; j < atoms.nb_atoms(); j++) {
            double r = (atoms.positions(Eigen::all,i)-atoms.positions(Eigen::all, j)).matrix().norm();
            sum += 24 * epsilon / r * (2 * std::pow(sigma / r, 12) - std::pow(sigma / r, 6));
        }
    }
}
