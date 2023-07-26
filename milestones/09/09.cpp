#include "domain.h"
#include "mpi.h"
#include "neighbors.h"
#include "xyz.h"
#include "verlet.h"
#include "ducastelle.h"
#include "thermostat.h"

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    constexpr double timestep = 5.0;
    constexpr int steps = 50001;
    constexpr int snapshot_interval = steps / 100; // 100 total frames
    constexpr double cutoff = 5.0;
    constexpr int eq_steps = 2000;
    constexpr double eq_temp = 100.0;
    constexpr double eq_relax = 10.0;

    NeighborList neighborList;

    // Perform IO on one rank to avoid duplicates
    std::ofstream *traj = nullptr;
    if (rank == 0) {
        traj = new std::ofstream ("traj.xyz");
        std::cout << "Step,Time,Total Energy,Potential,Kinetic,Temperature" << std::endl;
        std::cout.precision(10);
    }

    auto [names, init_positions]{read_xyz("whisker_small.xyz")};
    double gold_mass =
        20413.15887; // Gold in system's mass units where g/mol = 0.009649
    Atoms atoms(init_positions, gold_mass);
    atoms.k_b = 8.617333262e-5; // Boltzmann constant in eV/K
    Domain domain(MPI_COMM_WORLD, {40.39, 40.8, 144.24978336}, {1, 1, 4}, {0, 0, 1});
    domain.enable(atoms);
    domain.exchange_atoms(atoms);
    domain.update_ghosts(atoms, 2 * cutoff);
    std::cout << domain.nb_local() << " atoms for rank " << rank << std::endl;

    for (int i = 0; i < steps; i++) {
        // TODO masses???
        auto len = domain.domain_length();
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep, atoms.mass);
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, 2 * cutoff);
        neighborList.update(atoms, cutoff);
        double pot = ducastelle(atoms, neighborList, cutoff);
        verlet_step2(atoms.velocities, atoms.forces, timestep, atoms.mass);
        domain.scale(atoms, {len[0], len[1], len[2] + 0.001});

        if (i % snapshot_interval == 0) {
            domain.disable(atoms);
            if (traj != nullptr) {
                assert(rank == 0);
                write_xyz(*traj, atoms);
                auto kinetic = kinetic_energy(atoms);
                std::cout << i << "," << (i) * timestep << "," << pot + kinetic << ","
                          << pot << "," << kinetic << "," << temperature(atoms) << std::endl;
            }
            domain.enable(atoms);
        }
    }
    //    domain.disable(atoms);

    if (rank == 0) {
        delete traj;
    }
    MPI_Finalize();
}