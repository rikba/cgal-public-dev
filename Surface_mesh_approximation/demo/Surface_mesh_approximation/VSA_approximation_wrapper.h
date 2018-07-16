#ifndef VSA_APPROXIMAITON_WRAPPER_H
#define VSA_APPROXIMAITON_WRAPPER_H

#include <CGAL/Variational_shape_approximation.h>
#include <CGAL/L2_metric_plane_proxy.h>
#include <CGAL/property_map.h>

template <typename TriangleMesh, typename GeomTraits>
class VSA_approximation_wrapper {
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Vector_3 Vector_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type Vertex_point_map;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > Face_area_map;
  typedef boost::associative_property_map<std::map<face_descriptor, Point_3> > Face_center_map;

#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Variational_shape_approximation<TriangleMesh, Vertex_point_map,
    CGAL::Default, GeomTraits, CGAL::Parallel_tag> L21_approx;
#else
  typedef CGAL::Variational_shape_approximation<TriangleMesh, Vertex_point_map,
    CGAL::Default, GeomTraits> L21_approx;
#endif
  typedef typename L21_approx::Error_metric L21_metric;

  typedef CGAL::VSA::L2_metric_plane_proxy<TriangleMesh> L2_metric;
#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Variational_shape_approximation<TriangleMesh, Vertex_point_map,
    L2_metric, GeomTraits, CGAL::Parallel_tag> L2_approx;
#else
  typedef CGAL::Variational_shape_approximation<TriangleMesh, Vertex_point_map,
    L2_metric, GeomTraits> L2_approx;
#endif

  // user defined point-wise compact metric
  struct Compact_metric_point_proxy {
    typedef Point_3 Proxy;

    Compact_metric_point_proxy(const Face_center_map &_center_pmap,
      const Face_area_map &_area_pmap)
      : center_pmap(_center_pmap), area_pmap(_area_pmap) {}

    FT compute_error(const TriangleMesh &tm, const face_descriptor f, const Proxy &px) const {
      (void)(tm);
      return FT(std::sqrt(CGAL::to_double(
        CGAL::squared_distance(center_pmap[f], px))));
    }

    template <typename FaceRange>
    Proxy fit_proxy(const FaceRange &faces, const TriangleMesh &tm) const {
      (void)(tm);
      CGAL_assertion(!faces.empty());

      // fitting center
      Vector_3 center = CGAL::NULL_VECTOR;
      FT area(0.0);
      BOOST_FOREACH(const face_descriptor f, faces) {
        center = center + (center_pmap[f] - CGAL::ORIGIN) * area_pmap[f];
        area += area_pmap[f];
      }
      center = center / area;
      return CGAL::ORIGIN + center;
    }

    const Face_center_map center_pmap;
    const Face_area_map area_pmap;
  };
  typedef Compact_metric_point_proxy Compact_metric;

#ifdef CGAL_LINKED_WITH_TBB
  typedef CGAL::Variational_shape_approximation<TriangleMesh, Vertex_point_map,
    Compact_metric, GeomTraits, CGAL::Parallel_tag> Compact_approx;
#else
  typedef CGAL::Variational_shape_approximation<TriangleMesh, Vertex_point_map,
    Compact_metric, GeomTraits> Compact_approx;
#endif

public:
  enum Metric { L21, L2, Compact };

  typedef CGAL::cpp11::array<std::size_t, 3> Indexed_triangle;

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  typedef typename L21_approx::Proxy_wrapper L21_proxy_wrapper;
#endif

  VSA_approximation_wrapper()
    : m_metric(L21),
    m_center_pmap(m_face_centers),
    m_area_pmap(m_face_areas),
    m_pl21_metric(NULL),
    m_l21_approx(NULL),
    m_pl2_metric(NULL),
    m_l2_approx(NULL),
    m_pcompact_metric(NULL),
    m_iso_approx(NULL) {}

  ~VSA_approximation_wrapper() {
    if (m_l21_approx)
      delete m_l21_approx;
    if (m_pl21_metric)
      delete m_pl21_metric;
    if (m_l2_approx)
      delete m_l2_approx;
    if (m_pl2_metric)
      delete m_pl2_metric;
    if (m_iso_approx)
      delete m_iso_approx;
    if (m_pcompact_metric)
      delete m_pcompact_metric;
  }

  void set_mesh(const TriangleMesh &mesh) {
    Vertex_point_map vpm = get(boost::vertex_point, const_cast<TriangleMesh &>(mesh));

    m_face_centers.clear();
    m_face_areas.clear();
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      const halfedge_descriptor he = halfedge(f, mesh);
      const Point_3 &p0 = vpm[source(he, mesh)];
      const Point_3 &p1 = vpm[target(he, mesh)];
      const Point_3 &p2 = vpm[target(next(he, mesh), mesh)];

      m_face_centers.insert(std::pair<face_descriptor, Point_3>(
        f, CGAL::centroid(p0, p1, p2)));

      FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
      m_face_areas.insert(std::pair<face_descriptor, FT>(f, area));
    }

    if (m_l21_approx)
      delete m_l21_approx;
    if (m_pl21_metric)
      delete m_pl21_metric;
    if (m_l2_approx)
      delete m_l2_approx;
    if (m_pl2_metric)
      delete m_pl2_metric;
    if (m_iso_approx)
      delete m_iso_approx;
    if (m_pcompact_metric)
      delete m_pcompact_metric;

    m_pl21_metric = new L21_metric(mesh, vpm);
    m_l21_approx = new L21_approx(mesh, vpm, *m_pl21_metric);

    m_pl2_metric = new L2_metric(mesh, vpm);
    m_l2_approx = new L2_approx(mesh, vpm, *m_pl2_metric);

    m_pcompact_metric = new Compact_metric(m_center_pmap, m_area_pmap);
    m_iso_approx = new Compact_approx(mesh, vpm, *m_pcompact_metric);
  }

  void set_metric(const Metric &m) { m_metric = m; }

  std::size_t initialize_seeds(const CGAL::VSA::Seeding_method method,
    const boost::optional<std::size_t> max_nb_of_proxies,
    const boost::optional<FT> min_error_drop,
    const std::size_t nb_relaxations) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->initialize_seeds(
          method, max_nb_of_proxies, min_error_drop, nb_relaxations);
      case L2:
        return m_l2_approx->initialize_seeds(
          method, max_nb_of_proxies, min_error_drop, nb_relaxations);
      case Compact:
        return m_iso_approx->initialize_seeds(
          method, max_nb_of_proxies, min_error_drop, nb_relaxations);
    }
    return 0;
  }

  void run(const std::size_t nb_iterations) {
    FT err(0.0);
    switch (m_metric) {
      case L21:
        m_l21_approx->run(nb_iterations);
        break;
      case L2:
        m_l2_approx->run(nb_iterations);
        break;
      case Compact:
        m_iso_approx->run(nb_iterations);
        break;
    }
  }

  std::size_t add_one_proxy() {
    switch (m_metric) {
      case L21:
        return m_l21_approx->add_to_furthest_proxies(1, 0);
      case L2:
        return m_l2_approx->add_to_furthest_proxies(1, 0);
      case Compact:
        return m_iso_approx->add_to_furthest_proxies(1, 0);
    }
    return 0;
  }

  std::size_t teleport_one_proxy() {
    switch (m_metric) {
      case L21:
        return m_l21_approx->teleport_proxies(1, 0, true);
      case L2:
        return m_l2_approx->teleport_proxies(1, 0, true);
      case Compact:
        return m_iso_approx->teleport_proxies(1, 0, true);
    }
    return 0;
  }

  bool split(const std::size_t px_idx, const std::size_t n, const std::size_t nb_relaxations) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->split(px_idx, n, nb_relaxations);
      case L2:
        return m_l2_approx->split(px_idx, n, nb_relaxations);
      case Compact:
        return m_iso_approx->split(px_idx, n, nb_relaxations);
    }
    return false;
  }

  bool extract_mesh(const FT subdivision_ratio,
    const bool relative_to_chord,
    const bool with_dihedral_angle,
    const bool optimize_anchor_location,
    const bool pca_plane) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->extract_mesh(
          CGAL::VSA::parameters::subdivision_ratio(subdivision_ratio).
          relative_to_chord(relative_to_chord).
          with_dihedral_angle(with_dihedral_angle).
          optimize_anchor_location(optimize_anchor_location).
          pca_plane(pca_plane));
      case L2:
        return m_l2_approx->extract_mesh(
          CGAL::VSA::parameters::subdivision_ratio(subdivision_ratio).
          relative_to_chord(relative_to_chord).
          with_dihedral_angle(with_dihedral_angle).
          optimize_anchor_location(optimize_anchor_location).
          pca_plane(pca_plane));
      case Compact:
        return m_iso_approx->extract_mesh(
          CGAL::VSA::parameters::subdivision_ratio(subdivision_ratio).
          relative_to_chord(relative_to_chord).
          with_dihedral_angle(with_dihedral_angle).
          optimize_anchor_location(optimize_anchor_location).
          pca_plane(pca_plane));
    }
    return false;
  }

  std::size_t proxies_size() {
    switch (m_metric) {
      case L21:
        return m_l21_approx->proxies_size();
      case L2:
        return m_l2_approx->proxies_size();
      case Compact:
        return m_iso_approx->proxies_size();
    }
    return 0;
  }

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  template <typename OutputIterator>
  void get_l21_proxies(OutputIterator outitr) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->wrapped_proxies(outitr);
      default:
        return;
    }
  }
#endif

  template <typename FaceProxyMap>
  void proxy_map(FaceProxyMap &fpmap) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->proxy_map(fpmap);
      case L2:
        return m_l2_approx->proxy_map(fpmap);
      case Compact:
        return m_iso_approx->proxy_map(fpmap);
    }
  }

  template <typename OutputIterator>
  void indexed_triangles(OutputIterator outitr) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->indexed_triangles(outitr);
      case L2:
        return m_l2_approx->indexed_triangles(outitr);
      case Compact:
        return m_iso_approx->indexed_triangles(outitr);
    }
  }

  template <typename OutputIterator>
  void anchor_points(OutputIterator outitr) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->anchor_points(outitr);
      case L2:
        return m_l2_approx->anchor_points(outitr);
      case Compact:
        return m_iso_approx->anchor_points(outitr);
    }
  }

  template <typename OutputIterator>
  void anchor_vertices(OutputIterator outitr) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->anchor_vertices(outitr);
      case L2:
        return m_l2_approx->anchor_vertices(outitr);
      case Compact:
        return m_iso_approx->anchor_vertices(outitr);
    }
  }

  template <typename OutputIterator>
  void indexed_boundary_polygons(OutputIterator outitr) {
    switch (m_metric) {
      case L21:
        return m_l21_approx->indexed_boundary_polygons(outitr);
      case L2:
        return m_l2_approx->indexed_boundary_polygons(outitr);
      case Compact:
        return m_iso_approx->indexed_boundary_polygons(outitr);
    }
  }

private:
  Metric m_metric; // current metric

  // face property maps
  std::map<face_descriptor, Point_3> m_face_centers;
  Face_center_map m_center_pmap;
  std::map<face_descriptor, FT> m_face_areas;
  Face_area_map m_area_pmap;

  L21_metric *m_pl21_metric;
  L21_approx *m_l21_approx;

  L2_metric *m_pl2_metric;
  L2_approx *m_l2_approx;

  Compact_metric *m_pcompact_metric;
  Compact_approx *m_iso_approx;
};

#endif // VSA_APPROXIMAITON_WRAPPER_H
