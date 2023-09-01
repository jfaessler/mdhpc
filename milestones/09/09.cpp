#include "analysis.h"
#include "domain.h"
#include "ducastelle.h"
#include "mpi.h"
#include "neighbors.h"
#include "thermostat.h"
#include "verlet.h"
#include "xyz.h"
#include <iostream>

int main(int argc, char *argv[]) {
    std::string filename;
    std::string out_filename;
    int domain_x;
    int domain_y;
    int domain_z;
    if (argc > 5) {
        filename = argv[1];
        out_filename = argv[2];
        domain_x = static_cast<int>(strtol(argv[3], nullptr, 10));
        domain_y = static_cast<int>(strtol(argv[4], nullptr, 10));
        domain_z = static_cast<int>(strtol(argv[5], nullptr, 10));
    }
    else {
        filename = "whisker_small.xyz";
        out_filename = "traj.xyz";
        domain_x = 2;
        domain_y = 2;
        domain_z = 10;
    }
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    constexpr double timestep = 5.0;
    constexpr int steps = 80001;
    constexpr int snapshot_interval = steps / 100; // 100 total frames
    constexpr double cutoff = 10.0;
    constexpr int eq_steps = 10000;
    constexpr double eq_temp = 500.0;
    constexpr double eq_relax = 20.0;
    constexpr double strain_rate = 0.0015; // Strain per frame
    double original_length;

    NeighborList neighborList;

    auto [names, init_positions]{read_xyz(filename)};
    double gold_mass =
        20413.15887; // Gold in system's mass units where g/mol = 0.009649
    Atoms atoms(init_positions, gold_mass);
    atoms.k_b = 8.617333262e-5; // Boltzmann constant in eV/K
    Domain domain(MPI_COMM_WORLD, {40.38993934, 49.799999999, 144.24978336}, {domain_x, domain_y, domain_z},
                  {0, 0, 1});
    domain.enable(atoms);
    domain.exchange_atoms(atoms);
    domain.update_ghosts(atoms, 2 * cutoff);
    std::cout << domain.nb_local() << " atoms for rank " << rank << std::endl;
    original_length = domain.domain_length(2);

    // Perform IO on one rank to avoid duplicates
    std::ofstream *traj = nullptr;
    if (rank == 0) {
        traj = new std::ofstream(out_filename);
        std::cout << "Step,Time,Total "
                     "Energy,Potential,Kinetic,Temperature,Strain,Stress"
                  << std::endl;
        std::cout.precision(10);
    }

    for (int i = 0; i < steps; i++) {
        // TODO masses???
        auto len = domain.domain_length();
        verlet_step1(atoms.positions, atoms.velocities, atoms.forces, timestep,
                     atoms.mass);
        domain.exchange_atoms(atoms);
        domain.update_ghosts(atoms, 2 * cutoff);
        neighborList.update(atoms, cutoff);
        double pot = ducastelle(atoms, neighborList, cutoff);
        verlet_step2(atoms.velocities, atoms.forces, timestep, atoms.mass);
        if (i < eq_steps) {
            if (i < eq_steps / 2 && i % (eq_steps / 20) == 0) {
                domain.disable(atoms);
                    if (rank == 0) {
                        auto t = temperature(atoms);
                        atoms.velocities *= sqrt((t + 150) / t);
                    }
                domain.enable(atoms);
            }
        } else {
            domain.scale(atoms, {len[0], len[1], len[2] + strain_rate});
        }

        if (i % snapshot_interval == 0) {
            domain.exchange_atoms(atoms);
            domain.update_ghosts(atoms, 2 * cutoff);
            neighborList.update(atoms, cutoff);
            Eigen::Matrix3d local_stress =
                stress_parallel(atoms, neighborList, domain, cutoff);
            Eigen::Matrix3Xd recv(3, domain.size() * 3);
            Eigen::ArrayXi recvcount(3 * domain.size());
            recvcount.setConstant(9);
            Eigen::ArrayXi displ(3 * domain.size());
            for (int d = 0; d < domain.size(); d++) {
                displ[d] = 9 * d;
            }
            MPI_Allgatherv(local_stress.data(), 9, MPI_DOUBLE, recv.data(),
                           recvcount.data(), displ.data(), MPI_DOUBLE,
                           MPI_COMM_WORLD);
            Eigen::Matrix3d stress;
            stress.setZero();
            if (domain.rank() == 0) {
                for (int d = 0; d < domain.size(); d++) {
                    stress += recv.block<3, 3>(0, d * 3);
                }
            }
            domain.disable(atoms);
            if (traj != nullptr) {
                assert(rank == 0);
                write_xyz(*traj, atoms);
                auto kinetic = kinetic_energy(atoms);
                std::cout << i << "," << (i)*timestep << "," << pot + kinetic
                          << "," << pot << "," << kinetic << ","
                          << temperature(atoms) << ","
                          << (len[2] - original_length) / original_length << ","
                          << stress(2, 2) << std::endl;
            }
            domain.enable(atoms);
        }
    }

    if (rank == 0) {
        delete traj;
    }
    MPI_Finalize();
}
