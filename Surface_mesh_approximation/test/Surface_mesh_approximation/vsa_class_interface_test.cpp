#include <iostream>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/IO/Polyhedron_iostream.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/Variational_shape_approximation.h>
#include <CGAL/L2_metric_plane_proxy.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron;
typedef Polyhedron::Facet_handle Facet_handle;
typedef boost::associative_property_map<std::map<Facet_handle, std::size_t> > Facet_proxy_map;
typedef boost::property_map<Polyhedron, boost::vertex_point_t>::type Vertex_point_map;

typedef CGAL::VSA::L2_metric_plane_proxy<Polyhedron> L2_metric_plane_proxy;
typedef CGAL::Variational_shape_approximation<Polyhedron, Vertex_point_map, L2_metric_plane_proxy> L2_approx;
typedef L2_approx::Proxy Plane_proxy;

/**
 * This file tests the main class API and the L2 metric.
 * It should cover all the APIs.
 */
int main()
{
  Polyhedron mesh;
  std::ifstream input("./data/sphere_iso.off");
  if (!input || !(input >> mesh) || mesh.empty()) {
    std::cerr << "Invalid off file." << std::endl;
    return EXIT_FAILURE;
  }

  // facet area map
  std::map<Facet_handle, std::size_t> facet_index;
  for (Polyhedron::Facet_iterator fitr = mesh.facets_begin();
    fitr != mesh.facets_end(); ++fitr)
    facet_index.insert(std::pair<Facet_handle, std::size_t>(fitr, 0));
  Facet_proxy_map proxy_pmap(facet_index);

  // create L2_approx L2 metric approximation algorithm instance
  std::cout << "setup algorithm instance" << std::endl;
  L2_approx approx(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)));
  L2_metric_plane_proxy error_metric(mesh,
    get(boost::vertex_point, const_cast<Polyhedron &>(mesh)));
  approx.set_metric(error_metric);

  // random seeding and run
  std::cout << "random seeding and run" << std::endl;
  approx.seeding(CGAL::VSA::Random, 10);
  approx.run(10);
  if (approx.proxies_size() != 10)
    return EXIT_FAILURE;

  // incremental add and run to convergence
  std::cout << "incremental add and run to convergence" << std::endl;
  approx.add_to_furthest_proxies(3, 5);
  if (approx.run_to_convergence(0.1))
    std::cout << "Converged." << std::endl;
  if (approx.proxies_size() != 13)
    return EXIT_FAILURE;

  std::cout << "hierarchical add and run" << std::endl;
  approx.add_proxies_error_diffusion(3);
  approx.run(10);
  if (approx.proxies_size() != 16)
    return EXIT_FAILURE;

  // merge and teleport the proxies from local minimal
  std::cout << "teleport" << std::endl;
  approx.teleport_proxies(3);
  approx.run(10);
  if (approx.proxies_size() != 16)
    return EXIT_FAILURE;

  // split proxy 0 into 2 proxies
  // precondition: proxy 0 should have more than 2 facets
  std::cout << "spliting" << std::endl;
  if (!approx.split(0, 2, 10))
    return EXIT_FAILURE;
  if (approx.proxies_size() != 17)
    return EXIT_FAILURE;

  // extract the approximation polyhedron
  std::cout << "meshing" << std::endl;
  if (approx.extract_mesh(CGAL::VSA::parameters::subdivision_ratio(1.0)))
    std::cout << "manifold." << std::endl;
  else
    std::cout << "non-manifold" << std::endl;

  // get outputs
  std::cout << "get outputs" << std::endl;
  approx.proxy_map(proxy_pmap);

  for (std::size_t i = 0; i < approx.proxies_size(); ++i) {
    std::list<Facet_handle> patch;
    approx.proxy_region(i, std::back_inserter(patch));
  }

  std::vector<Plane_proxy> proxies;
  approx.proxies(std::back_inserter(proxies));

  std::vector<Point> anchor_pos;
  approx.anchor_points(std::back_inserter(anchor_pos));

  std::vector<Polyhedron::Vertex_handle> anchor_vtx;
  approx.anchor_vertices(std::back_inserter(anchor_vtx));

  std::vector<CGAL::cpp11::array<std::size_t, 3> > tris;
  approx.indexed_triangles(std::back_inserter(tris));

  std::vector<std::vector<std::size_t> > boundary;
  approx.indexed_boundary_polygons(std::back_inserter(boundary));

  const FT drop(0.001);
  const std::size_t iterations = 5;
  std::cout << "rebuild and hierarchical seeding" << std::endl;
  approx.rebuild();
  if (approx.proxies_size() != 0)
    return EXIT_FAILURE;
  approx.seeding(CGAL::VSA::Hierarchical, boost::none, drop, iterations);
  approx.run(10);
  std::cout << "#proxies " << approx.proxies_size() << std::endl;

  std::cout << "rebuild and incremental seeding" << std::endl;
  approx.rebuild();
  if (approx.proxies_size() != 0)
    return EXIT_FAILURE;
  approx.seeding(CGAL::VSA::Incremental, boost::none, drop, iterations);
  approx.run(10);
  std::cout << "#proxies " << approx.proxies_size() << std::endl;

  return EXIT_SUCCESS;
}
