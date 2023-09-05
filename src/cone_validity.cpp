#include "cone_validity.h"
#include <iostream>
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
#include <igl/gaussian_curvature.h>


// Check if vertex is an interior cone
bool is_interior_cone(
  const std::vector<bool> &is_boundary_vertex,
  const std::vector<bool> &is_cone,
  size_t vertex
) {
  return ((!is_boundary_vertex[vertex]) && (is_cone[vertex]));
}


// Check if this is a straight boundary vertex
bool is_boundary_cone(
  const std::vector<bool> &is_boundary_vertex,
  const std::vector<bool> &is_cone,
  size_t vertex
) {
  return ((is_boundary_vertex[vertex]) && (is_cone[vertex]));
}


// Check if vertex is k edges away from a cone
bool is_k_adjacent_to_cone_or_boundary(
  const std::vector<std::vector<size_t>> &adjacency_list,
  const std::vector<bool> &is_boundary_vertex,
  const std::vector<bool> &is_cone,
  size_t vertex,
  size_t k
) {
  // In base case 0, just check if the vertex is a cone or boundary
  if (k == 0)
  {
    return ((is_cone[vertex]) || (is_boundary_vertex[vertex]));
  }

  // Vertex is k > 0 adjacent to a cone if and only if an adjacent vertex is k - 1 adjacent or...
  for (size_t i = 0; i < adjacency_list[vertex].size(); ++i)
  {
    size_t vi = adjacency_list[vertex][i];
    if (is_k_adjacent_to_cone_or_boundary(adjacency_list, is_boundary_vertex, is_cone, vi, k - 1))
    {
      return true;
    }
  }

  // ...if the vertex is a cone or boundary itself
  return ((is_cone[vertex]) || (is_boundary_vertex[vertex]));
}

