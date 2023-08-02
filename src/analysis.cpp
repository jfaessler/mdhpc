#include "analysis.h"
Eigen::Matrix3d stress_parallel(const Atoms &atoms, const NeighborList &neighbor_list,
                       const Domain &domain, double cutoff, double A, double xi,
                       double p, double q, double re) {
    auto cutoff_sq{cutoff * cutoff};
    double xi_sq{xi * xi};

    // Initialize 0 matrix to set value
    Eigen::Matrix3d stress;
    stress.setZero();

    // compute embedding energies
    Eigen::ArrayXd embedding(
        atoms.nb_atoms()); // contains first density, later energy
    embedding.setZero();
    for (auto [i, j] : neighbor_list) {
        if (i < j) {
            Eigen::Vector3d distance_vector{atoms.positions.col(i) -
                                            atoms.positions.col(j)};
            auto distance_sq = distance_vector.squaredNorm();
            if (distance_sq < cutoff_sq) {
                double density_contribution{
                    xi_sq *
                    std::exp(-2 * q * (std::sqrt(distance_sq) / re - 1.0))};
                embedding(i) += density_contribution;
                embedding(j) += density_contribution;
            }
        }
    }

    // compute embedding contribution to the potential energy
    embedding = -embedding.sqrt();

    // per-atom energies
    Eigen::ArrayXd energies{embedding};

    // compute forces
    for (auto [i, j] : neighbor_list) {
        // Only compute within the subdomain
        if (i < j && i < domain.nb_local()) {
            double d_embedding_density_i{0};
            // this is the derivative of sqrt(embedding)
            if (embedding(i) != 0)
                d_embedding_density_i = 1 / (2 * embedding(i));

            Eigen::Vector3d distance_vector{atoms.positions.col(i) -
                                            atoms.positions.col(j)};
            auto distance_sq = distance_vector.squaredNorm();
            if (distance_sq < cutoff_sq) {
                double distance{std::sqrt(distance_sq)};
                double d_embedding_density_j{0};
                // this is the derivative of sqrt(embedding)
                if (embedding(j) != 0)
                    d_embedding_density_j = 1 / (2 * embedding(j));

                // repulsive energy and derivative of it with respect to
                // distance
                double repulsive_energy{2 * A *
                                        std::exp(-p * (distance / re - 1.0))};
                double d_repulsive_energy{-repulsive_energy * p / re};

                // derivative of embedding energy contributions
                double fac{-2 * q / re * xi_sq *
                           std::exp(-2 * q * (distance / re - 1.0))};

                // pair force
                Eigen::Array3d pair_force{
                    (d_repulsive_energy +
                     fac * (d_embedding_density_i + d_embedding_density_j)) *
                    distance_vector.normalized()};

                stress += distance_vector * pair_force.matrix().transpose();
            }
        }
    }

    Eigen::Array3d len = domain.domain_length();
    // Divide stress by volume of entire domain, so that we can add each rank
    stress /= len[0] * len[1] * len[2];
    return stress;
}
