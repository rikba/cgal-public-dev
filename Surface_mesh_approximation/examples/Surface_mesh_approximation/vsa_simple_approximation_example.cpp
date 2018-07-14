#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/approximate_mesh.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Polyhedron_3<Kernel> Polyhedron;

int main()
{
  // read input surface triangle mesh 
  Polyhedron input;
  std::ifstream file("data/mask.off");
  file >> input;

  // The output will be an indexed triangle mesh
  std::vector<Kernel::Point_3> anchors;
  std::vector<CGAL::cpp11::array<std::size_t, 3> > triangles;

  // free function interface with named parameters
  CGAL::approximate_mesh(input,
    CGAL::VSA::parameters::verbose_level(CGAL::Main_steps).
    max_nb_proxies(200).
    anchors(std::back_inserter(anchors)). // anchor points
    triangles(std::back_inserter(triangles))); // indexed triangles

  std::cout << "#anchor points: " << anchors.size() << std::endl;
  std::cout << "#triangles: " << triangles.size() << std::endl;

  return EXIT_SUCCESS;
}
