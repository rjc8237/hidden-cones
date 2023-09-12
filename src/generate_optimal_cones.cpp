#include <iostream>
#include <random>
#include <queue>
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
#include <igl/adjacency_list.h>
#include <igl/partition.h>
#include <igl/boundary_loop.h>
#include <igl/is_border_vertex.h>
#include <igl/edges.h>
#include "cone_validity.h"
#include "cone_energy.h"


/**
 * Helper function for compute the topological info (#bd, #genus) of the input mesh
 * @param V dim #v*3 matrix, each row corresponds to mesh vertex coordinates
 * @param F dim #f*3 matrix, each row corresponds to three vertex ids of each facet
 * @return n_genus, int, genus of the input mesh
 * @return n_boundary, int, number of boundary loops of the input mesh
 */ 
static
std::pair<int,int> count_genus_and_boundary(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F){
  int nv = V.rows();
  int nf = F.rows();
  Eigen::MatrixXi E;
  igl::edges(F, E);

  std::vector<std::vector<int>> bds;
  igl::boundary_loop(F, bds);

  int n_bd = bds.size();
  int ne = E.rows();
  int eu = nv - ne + nf + n_bd;
  int n_genus = (2-eu)/2;
  return std::make_pair(n_genus, n_bd);
}

bool add_optimal_cone(
  const Eigen::VectorXi &candidate_cones,
  const std::vector<std::vector<size_t>> &adjacency_list,
  const std::vector<bool> &is_boundary_vertex,
  std::vector<bool> &is_cone,
  std::vector<double> &cone_energy,
  int &cone_deficit
) {
  // Build a priority queue of cones with energy
  typedef std::pair<double, int> cone_pair;
  std::priority_queue<
    cone_pair,
    std::vector<cone_pair>,
    std::greater<cone_pair>
  > optimal_candidate_cones;
  for (size_t i = 0; i < candidate_cones.size(); ++ i)
  {
    int vi = candidate_cones[i];
    optimal_candidate_cones.push(std::make_pair(cone_energy[vi], vi));
  }

  // Try all possible cones
  while (!optimal_candidate_cones.empty())
  {
    // Get a potential cone update
    int cone_vertex = optimal_candidate_cones.top().second;
    spdlog::debug(
      "Attempting to add cone at {}",
      cone_vertex
    );

    // Check if this is a valid update
    bool invalid = is_k_adjacent_to_cone_or_boundary(
      adjacency_list,
      is_boundary_vertex,
      is_cone,
      cone_vertex,
      2
    );

    // Continue if it is not valid or perform the update if it is
    if (invalid)
    {
      spdlog::debug("Rejecting cone");
      optimal_candidate_cones.pop();
      continue;
    } else {
      // Update the cone if it is valid
      is_cone[cone_vertex] = true;
      cone_deficit--;
      spdlog::debug(
        "Added cone at {}. Deficit of {} remaining",
        cone_vertex,
        cone_deficit
      );
      return true;
    }
  }

  // Failure if cannot place any cones
  return false;
}


void generate_boundary_cones(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const Eigen::VectorXi &boundary_vertices,
  int cone_index_deficit,
  std::vector<int> &cone_vertices,
  std::vector<double> &cone_angles
){
  size_t num_vertices = V.rows();
  size_t n_faces = F.rows();
  cone_vertices.clear();
  cone_angles.clear();

  // Initialize all angles to 2 pi
  cone_angles = std::vector<double>(num_vertices, 2 * M_PI);

  int n_bd_vertices = boundary_vertices.size();
  double cone_angle_deficit = 2.0 * M_PI * cone_index_deficit;
  double change_in_cone_angle = cone_angle_deficit / n_bd_vertices;
  spdlog::info(
    "Distributing cone angle {} deficit evenly by {} to boundary",
    cone_angle_deficit,
    change_in_cone_angle
  );
  for (size_t i = 0; i < n_bd_vertices; ++i)
  {
    int bvi = boundary_vertices[i];
    cone_angles[bvi] = M_PI + change_in_cone_angle;
    cone_vertices.push_back(bvi);
  }
}


//void generate_global_optimal_cones(
//  const Eigen::MatrixXd &V,
//  const Eigen::MatrixXi &F,
//  const Eigen::VectorXi &boundary_vertices,
//  int cone_index_deficit,
//  std::vector<int> &cone_vertices,
//  std::vector<double> &cone_angles
//){
//  size_t num_vertices = V.rows();
//  size_t num_faces = F.rows();
//  cone_vertices.clear();
//  cone_angles.clear();
//
//  // Get vertex one rings
//  std::vector<std::vector<size_t>> adjacency_list;
//  igl::adjacency_list(F, adjacency_list);
//
//  // Compute index of cones to place
//  int change_in_cone_index = (cone_index_deficit > 0) ? -1 : 1;
//  spdlog::info(
//    "Must place {} signed cones of index {}",
//    std::abs(cone_index_deficit),
//    change_in_cone_index
//  );
//
//  // First build cone index array with 0 for interior and 2 for boundary
//  std::vector<int> cone_indices(num_vertices, 4);
//  std::vector<bool> is_boundary_vertex(num_vertices, false);
//  for (size_t i = 0; i < boundary_vertices.size(); ++i)
//  {
//    int bvi = boundary_vertices[i];
//    cone_indices[bvi] = 2;
//    is_boundary_vertex[bvi] = true;
//  }
//
//  // Generate cone energy
//  std::vector<double> cone_energy;
//  ConeEnergyParameters cone_energy_params;
//  compute_cone_energy(
//    V,
//    F,
//    cone_indices,
//    cone_energy,
//    cone_energy_params
//  );
//
//  // Compute globally optimal cones
//  spdlog::info("Placing optimal cones");
//  Eigen::VectorXi all_vertices = arange(num_vertices);
//  while (cone_index_deficit != 0)
//  {
//    // Add a random cone
//    bool success = add_optimal_cone(
//      all_vertices,
//      adjacency_list,
//      is_boundary_vertex,
//      cone_indices,
//      cone_energy,
//      cone_index_deficit,
//      change_in_cone_index
//    );
//    if (!success) {
//      spdlog::error(
//        "Failed to find a valid cone assignment"
//      );
//      cone_vertices.clear();
//      cone_angles.clear();
//      return;
//    }
//  }
//
//  // Update cone angles according to the indices
//  cone_angles.resize(num_vertices);
//  for (size_t vi = 0; vi < num_vertices; ++vi)
//  {
//    // Change cone angles
//    cone_angles[vi] = (M_PI / 2.0) * cone_indices[vi];
//
//    // Record vertices with cones
//    if (cone_indices[vi] != 0)
//    {
//      cone_vertices.push_back(vi);
//    }
//  }
//}


void generate_partitioned_optimal_cones(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const Eigen::VectorXi &boundary_vertices,
  int cone_index_deficit,
  std::vector<int> &cone_vertices,
  std::vector<double> &cone_angles,
  int num_cones_per_index=4
) {
  size_t num_vertices = V.rows();
  size_t num_faces = F.rows();
  cone_vertices.clear();
  cone_angles.clear();
  spdlog::info("Using partitioned optimal cones");

  int cone_deficit = num_cones_per_index * std::abs<int>(cone_index_deficit);
  double cone_angle_change = -(2 * M_PI * cone_index_deficit) / cone_deficit;

  // Get vertex one rings
  std::vector<std::vector<size_t>> adjacency_list;
  igl::adjacency_list(F, adjacency_list);

  // Compute index of cones to place
  spdlog::info(
    "Must place {} cones",
    cone_deficit
  );

  // Nothing to do if 0 cones to place
  if (cone_deficit <= 0) return;

  // Partition the mesh into groups by vertex distances
  Eigen::VectorXi group_indices;
  Eigen::VectorXi seed_vertices;
  Eigen::VectorXd distances_squared;
  igl::partition(
    V,
    cone_deficit,
    group_indices,
    seed_vertices,
    distances_squared
  );

  // Build the groups
  std::vector<std::vector<int>> vec_groups(cone_deficit);
  std::vector<Eigen::VectorXi> groups(cone_deficit);
  for (size_t i = 0; i < cone_deficit; ++i)
  {
    vec_groups[i].reserve(2 * num_vertices / cone_deficit);
  }
  for (size_t i = 0; i < num_vertices; ++i)
  {
    vec_groups[group_indices[i]].push_back(i);
  }
  for (size_t i = 0; i < cone_deficit; ++i)
  {
    groups[i].resize(vec_groups[i].size());
    for (size_t j = 0; j < vec_groups[i].size(); ++j)
    {
      groups[i][j] = vec_groups[i][j];
    }
  }

  // First build cone index array with 0 for interior and 2 for boundary
  std::vector<bool> is_boundary_vertex(num_vertices, false);
  std::vector<bool> is_cone(num_vertices, false);
  for (size_t i = 0; i < boundary_vertices.size(); ++i)
  {
    int bvi = boundary_vertices[i];
    is_boundary_vertex[bvi] = true;
  }

  // Generate cone energy
  std::vector<double> cone_energy;
  ConeEnergyParameters cone_energy_params;
  compute_cone_energy(
    V,
    F,
    cone_energy,
    cone_energy_params
  );

  // Add a cone for each group
  int num_cones = cone_deficit;
  for (size_t i = 0; i < num_cones; ++i)
  {
    bool success = add_optimal_cone(
      groups[i],
      adjacency_list,
      is_boundary_vertex,
      is_cone,
      cone_energy,
      cone_deficit
    );
  }

  // Update cone angles according to the indices
  cone_angles.resize(num_vertices);
  for (size_t vi = 0; vi < num_vertices; ++vi)
  {
    // Change cone angles
    if (is_cone[vi])
    {
      cone_angles[vi] =  2.0 * M_PI + cone_angle_change;
      cone_vertices.push_back(vi);
    }
    else if (is_boundary_vertex[vi])
    {
      cone_angles[vi] = M_PI;
    }
    else
    {
      cone_angles[vi] = 2.0 * M_PI;
    }
  }
}


void generate_optimal_cones(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  std::vector<int> &cone_vertices,
  std::vector<double> &cone_angles,
  int num_cones_per_index
) {
  size_t num_vertices = V.rows();
  size_t n_faces = F.rows();
  cone_vertices.clear();
  cone_angles.clear();

  // Get genus and boundary
  auto gb = count_genus_and_boundary(V, F);
  int n_genus = gb.first;
  int n_bd = gb.second;
  spdlog::info(
    "Building cones for genus {} surface with {} boundary components",
    n_genus,
    n_bd
  );
  if (n_genus < 0)
  {
    spdlog::error("Impossible negative genus encountered");
    return;
  }
  if (n_bd < 0)
  {
    spdlog::error("Impossible negative boundary count encountered");
    return;
  }

  // Get cone index deficit
  int cone_index_deficit = (2 - 2 * n_genus - n_bd);

  // Get boundary vertices
  Eigen::VectorXi boundary_vertices;
  compute_boundary_vertices(F, boundary_vertices);
  size_t n_bd_vertices = boundary_vertices.size();
  spdlog::info("{} boundary vertices found", n_bd_vertices);

  // Optionally just handle boundary vertices with uniform angles
  //if (n_bd > 0)
  //{
  //  generate_boundary_cones(
  //    V,
  //    F,
  //    boundary_vertices,
  //    cone_index_deficit,
  //    cone_vertices,
  //    cone_angles
  //  );
  //  return;
  //}

  // Use partitioned cone selection otherwise
  generate_partitioned_optimal_cones(
    V,
    F,
    boundary_vertices,
    cone_index_deficit,
    cone_vertices,
    cone_angles,
    num_cones_per_index
  );

  spdlog::info("Found valid cone assignment");
}
