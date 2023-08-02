#ifndef MD_CODE_ANALYSIS_H
#define MD_CODE_ANALYSIS_H

#include "atoms.h"
#include "neighbors.h"
#include "domain.h"

Eigen::Matrix3d stress_parallel(const Atoms &atoms, const NeighborList &neighbor_list,
                       const Domain &domain, double cutoff = 10.0,
                       double A = 0.2061, double xi = 1.790, double p = 10.229,
                       double q = 4.036, double re = 4.079 / sqrt(2));
#endif // MD_CODE_ANALYSIS_H
