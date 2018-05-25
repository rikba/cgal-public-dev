#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/approximate_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef boost::graph_traits<Polyhedron>::face_descriptor Face_descriptor;
typedef boost::unordered_map<Face_descriptor, std::size_t> Facet_index_map;
typedef boost::associative_property_map<Facet_index_map> Facet_proxy_pmap;

int main()
{
  // reads input mesh 
  Polyhedron input;
  std::ifstream file("data/mask.off");
  file >> input;

  Facet_index_map fidx_map;
  BOOST_FOREACH(Face_descriptor f, faces(input))
    fidx_map[f] = 0;

  // facet proxy index property map
  Facet_proxy_pmap fpxmap(fidx_map);

  // free function interface with named parameters
  CGAL::approximate_mesh(input,
  CGAL::Surface_mesh_approximation::parameters::max_nb_proxies(200). // first stop criterion
    min_error_drop(0.05). // second stop criterion
    nb_of_iterations(30). // number of relaxation iterations after seeding
    facet_proxy_map(fpxmap)); // output facet-proxy map

  // TODO: iterates over segments and outputs to console

  return EXIT_SUCCESS;
}
