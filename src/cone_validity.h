#pragma once

#include <Eigen/Core>
#include "generate_flat_cones.h"



/// Check if a vertex is k adjacent to a cone or boundary
///
/// @param[in] adjacency_list: list of lists of vertex adjacencies
/// @param[in] is_boundary_vertex: boolean mask of boundary vertices
/// @param[in] is_cone: boolean mask of current cone indices
/// @param[in] vertex: vertex index to check
/// @param[in] k: distance to use
bool is_k_adjacent_to_cone_or_boundary(
  const std::vector<std::vector<size_t>> &adjacency_list,
  const std::vector<bool> &is_boundary_vertex,
  const std::vector<bool> &is_cone,
  size_t vertex,
  size_t k
);


// Generate range 0,...,n-1
inline Eigen::VectorXi arange(size_t n)
{
  Eigen::VectorXi v(n);
  for (size_t i = 0; i < n; ++i)
  {
    v[i] = i;
  }

  return v;
}
