#include "generate_flat_cones.h"
#include <iostream>
#include "spdlog/spdlog.h"
#include "spdlog/fmt/ostr.h"
#include <igl/boundary_facets.h>
#include <igl/facet_components.h>
#include <igl/unique.h>
#include <igl/edges.h>


void compute_boundary_vertices(
  const Eigen::MatrixXi &F,
  Eigen::VectorXi &B
) {
  Eigen::MatrixXi E;
  Eigen::VectorXi J, K, E_flat, IA, IC;
  igl::boundary_facets(F, E, J, K);
  spdlog::debug("Redundant boundary vertices: {}", E);

  E_flat.resize(E.size());
  E_flat.head(E.rows()) = E.col(0);
  E_flat.tail(E.rows()) = E.col(1);
  spdlog::debug("Flattened redundant boundary vertices: {}", E_flat);


  igl::unique(E_flat, B, IA, IC);
  spdlog::debug("Unique boundary vertices: {}", B);
}

void generate_flat_cones(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  std::vector<double> &cone_angles
) {
  cone_angles.clear();

  Eigen::MatrixXi E;
  igl::edges(F, E);
  int num_vertices = V.rows();
  int num_edges = E.rows();
  int num_faces = F.rows();
  spdlog::info(
    "Surface has {} vertices, {} edges, and {} faces",
    num_vertices,
    num_edges,
    num_faces
  );

  // Get number of connected components
  Eigen::VectorXi components;
  igl::facet_components(F, components);
  if (components.size() == 0)
  {
    spdlog::error("No components");
    return;
  }
  int num_components = components.maxCoeff() - components.minCoeff();
  spdlog::info("Surface has {} connected components", num_components);

  // Initialize all angles to 2 pi
  cone_angles = std::vector<double>(num_vertices, 2 * M_PI);

  // Set boundary vertices to pi
  Eigen::VectorXi boundary_vertices;
  compute_boundary_vertices(F, boundary_vertices);
  for (size_t i = 0; i < boundary_vertices.size(); ++i)
  {
    size_t vi = boundary_vertices[i];
    cone_angles[vi] = M_PI;
  }
}

