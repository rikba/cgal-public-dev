#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/approximate_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef typename boost::graph_traits<Polyhedron>::face_descriptor face_descriptor;
typedef boost::unordered_map<face_descriptor, std::size_t> Face_index_map;
typedef boost::associative_property_map<Face_index_map> Face_proxy_pmap;

/**
 * This file tests the free function CGAL::VSA::approximate_mesh.
 */
int main()
{
  Polyhedron input;
  std::ifstream file("data/sphere_iso.off");
  if (!file || !(file >> input) || input.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  Face_index_map fidx_map;
  BOOST_FOREACH(face_descriptor f, faces(input))
    fidx_map[f] = 0;
  Face_proxy_pmap fpxmap(fidx_map);

  std::vector<Kernel::Vector_3> proxies;

  // free function interface with named parameters
  CGAL::VSA::approximate_mesh(input,
    CGAL::VSA::parameters::seeding_method(CGAL::VSA::Hierarchical). // hierarchical seeding
    max_nb_of_proxies(200). // both maximum number of proxies stop criterion,
    min_error_drop(0.05). // and minimum error drop stop criterion are specified
    nb_of_iterations(30). // number of clustering iterations after seeding
    nb_of_relaxations(5). // number of relaxations in seeding
    face_proxy_map(fpxmap). // output indexed face set
    proxies(std::back_inserter(proxies))); // number of iterations after seeding

  return EXIT_SUCCESS;
}
