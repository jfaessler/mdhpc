#ifndef MD_CODE_LJ_DIRECT_SUMMATION_H
#define MD_CODE_LJ_DIRECT_SUMMATION_H

#include "atoms.h"
#include "neighbors.h"

double lj(Atoms &atoms, NeighborList &neighbor_list, double cutoff = 5.0,
          double epsilon = 1.0, double sigma = 1.0);

#endif // MD_CODE_LJ_DIRECT_SUMMATION_H
