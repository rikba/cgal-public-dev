#ifndef SCENE_H
#define SCENE_H

#include <QtOpenGL/qgl.h>
#include <iostream>
#include <cmath>
#include <cstdlib>

#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>

#include <CGAL/boost/graph/graph_traits_Polyhedron_3.h>
#include <CGAL/property_map.h>
#include "VSA_approximation_wrapper.h"

typedef CGAL::Simple_cartesian<double> Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;

typedef CGAL::Polyhedron_3<Kernel> Polyhedron_3;
typedef Polyhedron_3::Halfedge_handle Halfedge_handle;
typedef Polyhedron_3::Edge_iterator Edge_iterator;
typedef Polyhedron_3::Facet_handle Facet_handle;
typedef Polyhedron_3::Facet_iterator Facet_iterator;
typedef CGAL::Bbox_3 Bbox_3;

typedef VSA_approximation_wrapper<Polyhedron_3, Kernel> Approximation_wrapper;
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
typedef Approximation_wrapper::L21_proxy_wrapper L21_proxy_wrapper;
#endif

typedef Approximation_wrapper::Indexed_triangle Indexed_triangle;

class Scene
{
public:
  Scene() :
    m_pmesh(NULL),
    m_fidx_pmap(m_fidx_map),
    m_view_polyhedron(false),
    m_view_wireframe(false),
    m_view_boundary(false),
    m_view_proxies(false),
    m_view_anchors(false),
    m_view_approximation(false) {}

  ~Scene() {
    if (m_pmesh)
      delete m_pmesh;
  }

  void update_bbox();
  Bbox_3 bbox() { return m_bbox; }

  // file menu
  int open(QString filename);
  void save_approximation(const std::string &filename);

  // algorithms
  void set_metric(const int init);
  void seeding(const CGAL::Approximation_seeding_tag method,
    const boost::optional<std::size_t> num_proxies,
    const boost::optional<FT> min_error_drop,
    const std::size_t nb_relaxations,
    const std::size_t nb_iterations);
  void extract_mesh(const double chord_error,
    const bool is_relative_to_chord,
    const bool with_dihedral_angle,
    const bool if_optimize_anchor_location,
    const bool pca_plane);
  void run_one_step();
  void add_one_proxy();
  void teleport_one_proxy();
  void split(const std::size_t px_idx, const std::size_t n, const std::size_t nb_relaxations);

  // toggle view options
  void toggle_view_polyhedron() {
    m_view_polyhedron = !m_view_polyhedron;
  }

  void toggle_view_wireframe() {
    m_view_wireframe = !m_view_wireframe;
  }

  void toggle_view_boundary() {
    m_view_boundary = !m_view_boundary;
  }

  void toggle_view_proxies() {
    m_view_proxies = !m_view_proxies;
  }

  void toggle_view_anchors() {
    m_view_anchors = !m_view_anchors;
  }

  void toggle_view_approximation() {
    m_view_approximation = !m_view_approximation;
  }

  void draw();

private:
  Vector_3 normalize(const Vector_3& v) {
    return v / std::sqrt(v * v);
  }

  // pseudorandom number for proxy color mapping
  std::size_t rand_0_255() {
    return static_cast<std::size_t>(std::rand() % 255);
  }

  // rendering
  void render_polyhedron();
  void render_wireframe();
  void render_boundary();
  void render_anchors();
  void render_borders();
  void render_approximation();
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  void render_proxies();
#endif

private:
  // member data
  Bbox_3 m_bbox;
  Polyhedron_3 *m_pmesh;

  // property-map for segment-idx
  std::map<Facet_handle, std::size_t> m_fidx_map;
  boost::associative_property_map<std::map<Facet_handle, std::size_t> > m_fidx_pmap;

  // algorithm instance
  Approximation_wrapper m_approx;

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  std::vector<L21_proxy_wrapper> m_proxies;
#endif
  std::vector<std::size_t> m_px_color;
  std::vector<Point_3> m_anchor_pos;
  std::vector<Polyhedron_3::Vertex_handle> m_anchor_vtx;
  std::vector<std::vector<std::size_t> > m_bdrs; // anchor borders
  std::vector<Indexed_triangle> m_tris;

  // view options
  bool m_view_polyhedron;
  bool m_view_wireframe;
  bool m_view_boundary;
  bool m_view_proxies;
  bool m_view_anchors;
  bool m_view_approximation;
}; // end class Scene

#endif // SCENE_H
