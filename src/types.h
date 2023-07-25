#ifndef MD_CODE_TYPES_H
#define MD_CODE_TYPES_H

#include "Eigen/Dense"

using Positions_t = Eigen::Array3Xd;
using Velocities_t = Eigen::Array3Xd;
using Forces_t = Eigen::Array3Xd;
using Masses_t = Eigen::ArrayXd;
using Scalar_t = double;
using Names_t = std::vector<std::basic_string<char>>;

#endif // MD_CODE_TYPES_H
