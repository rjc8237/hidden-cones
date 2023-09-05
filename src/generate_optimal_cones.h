#pragma once

#include <Eigen/Core>
#include "generate_flat_cones.h"


/// @brief Given a VF mesh, generate cone angles that are optimal by some energy
///
/// @param[in] V: mesh vertices
/// @param[in] F: mesh faces
/// @param[out] cone_vertices: vertices with non-flat cones (including boundary)
/// @param[out] cone_angles: angles at all vertices (including flat vertices)
/// @param[in] num_cones_per_index: cones per index (default four for pi/2 angles)
void generate_optimal_cones(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  std::vector<int> &cone_vertices,
  std::vector<double> &cone_angles,
  int num_cones_per_index=4
);
