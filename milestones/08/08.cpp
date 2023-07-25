#include "mpi.h"
#include "domain.h"
#include "xyz.h"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    auto [names, init_positions]{read_xyz("cluster_923.xyz")};
    double gold_mass = 20413.15887; // Gold in system's mass units where g/mol = 0.009649
    Atoms atoms(init_positions, gold_mass);
    atoms.k_b = 8.617333262e-5; // Boltzmann constant in eV/K
    Domain domain (MPI_COMM_WORLD, {30, 30, 30}, {2, 2, 1}, {0,0,0});
    domain.enable(atoms);
    domain.exchange_atoms(atoms);
    domain.update_ghosts(atoms, 10.0);
    domain.disable(atoms);
    MPI_Finalize();
}