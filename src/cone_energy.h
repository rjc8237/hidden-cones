#pragma once

#include <Eigen/Core>
#include "generate_flat_cones.h"


struct ConeEnergyParameters {

};

/// @brief Given a VF mesh, compute a cone energy for each vertex
///
/// @param[in] V: mesh vertices
/// @param[in] F: mesh faces
/// @param[out] cone_energy: cone energy at all vertices
/// @param[in] cone_energy_params: parameters for the cone energy
void compute_cone_energy(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  std::vector<double> &cone_energy,
  const ConeEnergyParameters &cone_energy_params
);
