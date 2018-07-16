#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Variational_shape_approximation.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type Vertex_point_map;

typedef CGAL::Variational_shape_approximation<Polyhedron, Vertex_point_map> L21_approx;
typedef L21_approx::Error_metric L21_metric;

bool check_strict_ordering(const std::vector<FT> &error)
{
  if (error.empty()) {
    std::cout << "Empty error sequence." << std::endl;
    return false;
  }
  FT pre = error.front();
  for (std::vector<FT>::const_iterator itr = error.begin(); itr != error.end(); ++itr)
    if (pre < *itr)
      return false;

  return true;
}

/**
 * This file tests the decrease of the relaxing error on a sphere shape.
 */
int main()
{
  Polyhedron mesh;
  std::ifstream input("./data/sphere_iso.off");
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // algorithm instance
  L21_metric error_metric(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)));
  L21_approx approx(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)),
    error_metric);

  approx.initialize_seeds(CGAL::VSA::Random, 100);
  std::vector<FT> error;
  for (std::size_t i = 0; i < 30; ++i) {
    approx.run();
    error.push_back(approx.compute_total_error());
  }

  if (check_strict_ordering(error)) {
    std::cout << "Pass the decrease test." << std::endl;
    return EXIT_SUCCESS;
  }
  else {
    std::cout << "Failed the decrease test." << std::endl;
    return EXIT_FAILURE;
  }
}
