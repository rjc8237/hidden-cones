#pragma once

#include <Eigen/Core>


/// @brief Given a the faces of a VF mesh, generate a list of the boundary vertices
///
/// @param[in] F: mesh faces
/// @param[in] B: mesh boundary vertices
void compute_boundary_vertices(
  const Eigen::MatrixXi &F,
  Eigen::VectorXi &B
);


/// @brief Given a VF mesh, generate trivial flat cone angles.
///
/// Note that this method will not generally produce valid cone angles except.
///
/// @param[in] V: mesh vertices
/// @param[in] F: mesh faces
/// @param[out] cone_angles: flat angles at vertices
void generate_flat_cones(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  std::vector<double> &cone_angles
);
