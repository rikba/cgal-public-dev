#include <CGAL/vsa_approximation.h>
#include <CGAL/property_map.h>

template <typename TriangleMesh, typename GeomTraits>
class VSAWrapper {
  typedef typename GeomTraits::FT FT;
  typedef typename GeomTraits::Point_3 Point_3;
  typedef typename GeomTraits::Vector_3 Vector_3;

  typedef typename boost::graph_traits<TriangleMesh>::face_descriptor face_descriptor;
  typedef typename boost::graph_traits<TriangleMesh>::halfedge_descriptor halfedge_descriptor;

  typedef typename boost::property_map<TriangleMesh, boost::vertex_point_t>::type VertexPointMap;
  typedef boost::associative_property_map<std::map<face_descriptor, FT> > FacetAreaMap;
  typedef boost::associative_property_map<std::map<face_descriptor, Point_3> > FacetCenterMap;

  typedef typename CGAL::VSA_seeding VSA_seeding;

  typedef CGAL::VSA_approximation<TriangleMesh, VertexPointMap,
    CGAL::Default, CGAL::Default, GeomTraits> L21VSA;
  typedef typename L21VSA::ErrorMetric L21Metric;
  typedef typename L21VSA::ProxyFitting L21ProxyFitting;

  typedef CGAL::L2Metric<TriangleMesh> L2Metric;
  typedef CGAL::L2ProxyFitting<TriangleMesh> L2ProxyFitting;
  typedef CGAL::VSA_approximation<TriangleMesh, VertexPointMap,
    L2Metric, L2ProxyFitting, GeomTraits> L2VSA;

  // user defined point-wise compact metric
  struct CompactMetric {
    typedef Point_3 Proxy;

    CompactMetric(const FacetCenterMap &_center_pmap)
      : center_pmap(_center_pmap) {}

    FT operator()(const face_descriptor &f, const Proxy &px) const {
      return FT(std::sqrt(CGAL::to_double(
        CGAL::squared_distance(center_pmap[f], px))));
    }

    const FacetCenterMap center_pmap;
  };

  struct PointProxyFitting {
    typedef Point_3 Proxy;

    PointProxyFitting(const FacetCenterMap &_center_pmap,
      const FacetAreaMap &_area_pmap)
      : center_pmap(_center_pmap),
      area_pmap(_area_pmap) {}

    template <typename FacetIterator>
    Proxy operator()(const FacetIterator beg, const FacetIterator end) const {
      CGAL_assertion(beg != end);

      // fitting center
      Vector_3 center = CGAL::NULL_VECTOR;
      FT area(0);
      for (FacetIterator fitr = beg; fitr != end; ++fitr) {
        center = center + (center_pmap[*fitr] - CGAL::ORIGIN) * area_pmap[*fitr];
        area += area_pmap[*fitr];
      }
      center = center / area;
      return CGAL::ORIGIN + center;
    }

    const FacetCenterMap center_pmap;
    const FacetAreaMap area_pmap;
  };
  typedef CGAL::VSA_approximation<TriangleMesh, VertexPointMap,
    CompactMetric, PointProxyFitting, GeomTraits> CompactVSA;

public:
  enum Metric { L21, L2, Compact };

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  typedef typename L21VSA::ProxyWrapper L21Proxy;
#endif

  VSAWrapper()
    : m_metric(L21),
    m_center_pmap(m_facet_centers),
    m_area_pmap(m_facet_areas),
    m_pl21_metric(NULL),
    m_pl21_proxy_fitting(NULL),
    m_pl2_metric(NULL),
    m_pl2_proxy_fitting(NULL),
    m_pcompact_metric(NULL),
    m_pcompact_proxy_fitting(NULL) {}

  ~VSAWrapper() {
    if (m_pl21_metric)
      delete m_pl21_metric;
    if (m_pl21_proxy_fitting)
      delete m_pl21_proxy_fitting;
    if (m_pl2_metric)
      delete m_pl2_metric;
    if (m_pl2_proxy_fitting)
      delete m_pl2_proxy_fitting;
    if (m_pcompact_metric)
      delete m_pcompact_metric;
    if (m_pcompact_proxy_fitting)
      delete m_pcompact_proxy_fitting;
  }

  void set_mesh(const TriangleMesh &mesh) {
    VertexPointMap vpm = get(boost::vertex_point, const_cast<TriangleMesh &>(mesh));

    m_facet_centers.clear();
    m_facet_areas.clear();
    BOOST_FOREACH(face_descriptor f, faces(mesh)) {
      const halfedge_descriptor he = halfedge(f, mesh);
      const Point_3 &p0 = vpm[source(he, mesh)];
      const Point_3 &p1 = vpm[target(he, mesh)];
      const Point_3 &p2 = vpm[target(next(he, mesh), mesh)];

      m_facet_centers.insert(std::pair<face_descriptor, Point_3>(
        f, CGAL::centroid(p0, p1, p2)));

      FT area(std::sqrt(CGAL::to_double(CGAL::squared_area(p0, p1, p2))));
      m_facet_areas.insert(std::pair<face_descriptor, FT>(f, area));
    }

    if (m_pl21_metric)
      delete m_pl21_metric;
    if (m_pl21_proxy_fitting)
      delete m_pl21_proxy_fitting;
    if (m_pl2_metric)
      delete m_pl2_metric;
    if (m_pl2_proxy_fitting)
      delete m_pl2_proxy_fitting;
    if (m_pcompact_metric)
      delete m_pcompact_metric;
    if (m_pcompact_proxy_fitting)
      delete m_pcompact_proxy_fitting;

    m_pl21_metric = new L21Metric(mesh);
    m_pl21_proxy_fitting = new L21ProxyFitting(mesh);
    m_pl2_metric = new L2Metric(mesh);
    m_pl2_proxy_fitting = new L2ProxyFitting(mesh);
    m_pcompact_metric = new CompactMetric(m_center_pmap);
    m_pcompact_proxy_fitting = new PointProxyFitting(m_center_pmap, m_area_pmap);

    m_vsa_l21.set_mesh(mesh, vpm);
    m_vsa_l21.set_metric(*m_pl21_metric, *m_pl21_proxy_fitting);

    m_vsa_l2.set_mesh(mesh, vpm);
    m_vsa_l2.set_metric(*m_pl2_metric, *m_pl2_proxy_fitting);

    m_vsa_compact.set_mesh(mesh, vpm);
    m_vsa_compact.set_metric(*m_pcompact_metric, *m_pcompact_proxy_fitting);
  }

  void set_metric(const Metric m) {
    m_metric = m;
    switch (m_metric) {
      case L21: return m_vsa_l21.rebuild();
      case L2: return m_vsa_l2.rebuild();
      case Compact: return m_vsa_compact.rebuild();
    }
  }

  std::size_t init_by_number(const int init, const std::size_t num_seed, const std::size_t iterations) {
    switch (m_metric) {
      case L21:
        return m_vsa_l21.init_by_number(static_cast<VSA_seeding>(init), num_seed, iterations);
      case L2:
        return m_vsa_l2.init_by_number(static_cast<VSA_seeding>(init), num_seed, iterations);
      case Compact:
        return m_vsa_compact.init_by_number(static_cast<VSA_seeding>(init), num_seed, iterations);
    }
    return 0;
  }

  std::size_t init_by_error(const int init, const FT drop, const std::size_t iterations) {
    switch (m_metric) {
      case L21:
        return m_vsa_l21.init_by_error(static_cast<VSA_seeding>(init), drop, iterations);
      case L2:
        return m_vsa_l2.init_by_error(static_cast<VSA_seeding>(init), drop, iterations);
      case Compact:
        return m_vsa_compact.init_by_error(static_cast<VSA_seeding>(init), drop, iterations);
    }
    return 0;
  }

  std::size_t run_one_step() {
    switch (m_metric) {
      case L21:
        return m_vsa_l21.run_one_step();
      case L2:
        return m_vsa_l2.run_one_step();
      case Compact:
        return m_vsa_compact.run_one_step();
    }
    return 0;
  }

  std::size_t add_one_proxy() {
    switch (m_metric) {
      case L21:
        return m_vsa_l21.add_proxies_furthest(1, 0);
      case L2:
        return m_vsa_l2.add_proxies_furthest(1, 0);
      case Compact:
        return m_vsa_compact.add_proxies_furthest(1, 0);
    }
    return 0;
  }

  std::size_t teleport_one_proxy() {
    switch (m_metric) {
      case L21:
        return m_vsa_l21.teleport_proxies(1, 0, true);
      case L2:
        return m_vsa_l2.teleport_proxies(1, 0, true);
      case Compact:
        return m_vsa_compact.teleport_proxies(1, 0, true);
    }
    return 0;
  }

  template <typename PolyhedronSurface>
  bool meshing(PolyhedronSurface &mesh_out, const FT split = FT(0.2), bool pca_plane = false) {
    switch (m_metric) {
      case L21:
        return m_vsa_l21.meshing(mesh_out, split, pca_plane);
      case L2:
        return m_vsa_l2.meshing(mesh_out, split, pca_plane);
      case Compact:
        return m_vsa_compact.meshing(mesh_out, split, pca_plane);
    }
    return false;
  }

  std::size_t get_proxies_size() {
    switch (m_metric) {
      case L21:
        return m_vsa_l21.get_proxies_size();
      case L2:
        return m_vsa_l2.get_proxies_size();
      case Compact:
        return m_vsa_compact.get_proxies_size();
    }
    return 0;
  }

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  template <typename OutputIterator>
  void get_l21_proxies(OutputIterator outitr) {
    switch (m_metric) {
      case L21:
        return m_vsa_l21.get_wrapped_proxies(outitr);
      default:
        return;
    }
  }
#endif

  template <typename FacetProxyMap>
  void get_proxy_map(FacetProxyMap &fpmap) {
    switch (m_metric) {
      case L21:
        return m_vsa_l21.get_proxy_map(fpmap);
      case L2:
        return m_vsa_l2.get_proxy_map(fpmap);
      case Compact:
        return m_vsa_compact.get_proxy_map(fpmap);
    }
  }

  template <typename OutputIterator>
  void get_indexed_triangles(OutputIterator outitr) {
    switch (m_metric) {
      case L21:
        return m_vsa_l21.get_indexed_triangles(outitr);
      case L2:
        return m_vsa_l2.get_indexed_triangles(outitr);
      case Compact:
        return m_vsa_compact.get_indexed_triangles(outitr);
    }
  }

  template <typename OutputIterator>
  void get_anchor_points(OutputIterator outitr) {
    switch (m_metric) {
      case L21:
        return m_vsa_l21.get_anchor_points(outitr);
      case L2:
        return m_vsa_l2.get_anchor_points(outitr);
      case Compact:
        return m_vsa_compact.get_anchor_points(outitr);
    }
  }

  template <typename OutputIterator>
  void get_anchor_vertices(OutputIterator outitr) {
    switch (m_metric) {
      case L21:
        return m_vsa_l21.get_anchor_vertices(outitr);
      case L2:
        return m_vsa_l2.get_anchor_vertices(outitr);
      case Compact:
        return m_vsa_compact.get_anchor_vertices(outitr);
    }
  }

  template <typename OutputIterator>
  void get_indexed_boundary_polygons(OutputIterator outitr) {
    switch (m_metric) {
      case L21:
        return m_vsa_l21.get_indexed_boundary_polygons(outitr);
      case L2:
        return m_vsa_l2.get_indexed_boundary_polygons(outitr);
      case Compact:
        return m_vsa_compact.get_indexed_boundary_polygons(outitr);
    }
  }

private:
  Metric m_metric; // current metric

  // facet property maps
  std::map<face_descriptor, Point_3> m_facet_centers;
  FacetCenterMap m_center_pmap;
  std::map<face_descriptor, FT> m_facet_areas;
  FacetAreaMap m_area_pmap;

  L21Metric *m_pl21_metric;
  L21ProxyFitting *m_pl21_proxy_fitting;
  L21VSA m_vsa_l21;

  L2Metric *m_pl2_metric;
  L2ProxyFitting *m_pl2_proxy_fitting;
  L2VSA m_vsa_l2;

  CompactMetric *m_pcompact_metric;
  PointProxyFitting *m_pcompact_proxy_fitting;
  CompactVSA m_vsa_compact;
};
