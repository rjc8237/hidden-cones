#include "cone_energy.h"
#include <iostream>
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
#include <igl/gaussian_curvature.h>


void compute_cone_energy(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  std::vector<double> &cone_energy,
  const ConeEnergyParameters &cone_energy_params
) {
  size_t num_vertices = V.rows();
  spdlog::debug("Computing cone energy");

  // Compute the Guassian curvature at each vertex
  Eigen::VectorXd K;
  igl::gaussian_curvature(V, F, K);

  // Build the cone energy
  cone_energy.resize(num_vertices);
  for (size_t vi = 0; vi < num_vertices; ++vi)
  {
    // Strong weight for negative curvature
    if (K[vi] < 0)
    {
      cone_energy[vi] = 10.0 * std::abs(K[vi]);
    } else {
      cone_energy[vi] = std::abs(K[vi]);
    }
  }

  spdlog::debug("Computing cone energy");
}