#include "Scene.h"

#include <iostream>
#include <fstream>

#include <QApplication>
#include <QString>
#include <QTextStream>
#include <QFileInfo>
#include <QInputDialog>

#include <CGAL/Timer.h>
#include <CGAL/IO/Polyhedron_iostream.h>
#include <CGAL/centroid.h>

#include "ColorCheatSheet.h"

void Scene::update_bbox()
{
  if (m_pmesh == NULL) {
    std::cout << "failed (no polyhedron)." << std::endl;
    return;
  }
  
  std::cout << "Compute bbox...";

  m_bbox = CGAL::bbox_3(m_pmesh->points_begin(), m_pmesh->points_end());
  
  std::cout << "done (" << m_pmesh->size_of_facets()
    << " facets)" << std::endl;
}

int Scene::open(QString filename)
{
  QTextStream cerr(stderr);
  cerr << QString("Opening file \"%1\"\n").arg(filename);
  QApplication::setOverrideCursor(QCursor(::Qt::WaitCursor));

  QFileInfo fileinfo(filename);
  std::ifstream in(filename.toUtf8());

  if (!in || !fileinfo.isFile() || ! fileinfo.isReadable()) {
    std::cerr << "unable to open file" << std::endl;
    QApplication::restoreOverrideCursor();
    return -1;
  }

  if (m_pmesh != NULL)
    delete m_pmesh;

  // allocate new polyhedron
  m_pmesh = new Polyhedron_3;
  in >> *m_pmesh;
  if (!in) {
    std::cerr << "invalid OFF file" << std::endl;
    QApplication::restoreOverrideCursor();

    delete m_pmesh;
    m_pmesh = NULL;
    return -1;
  }

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  // degeneracy check
  std::cerr << "Degeneracy check." << std::endl;
  for (Facet_iterator fitr = m_pmesh->facets_begin(); fitr != m_pmesh->facets_end(); ++fitr) {
    Halfedge_handle he = fitr->halfedge();
    if (CGAL::collinear(
      he->opposite()->vertex()->point(),
      he->vertex()->point(),
      he->next()->vertex()->point()))
      std::cerr << "Warning: degenerate facet" << std::endl;
  }
  std::cerr << "Done." << std::endl;
#endif

  m_fidx_map.clear();
  for (Facet_iterator fitr = m_pmesh->facets_begin(); fitr != m_pmesh->facets_end(); ++fitr)
    m_fidx_map.insert(std::pair<Facet_handle, std::size_t>(fitr, 0));

  m_approx.set_mesh(*m_pmesh);
  m_approx.set_metric(Approximation_wrapper::L21);

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
#endif
  m_px_color.clear();
  m_anchor_pos.clear();
  m_anchor_vtx.clear();
  m_bdrs.clear();
  m_tris.clear();

  m_view_polyhedron = true;
  m_view_wireframe = false;
  m_view_boundary = false;
  m_view_proxies = false;
  m_view_anchors = false;
  m_view_approximation = false;

  QApplication::restoreOverrideCursor();
  return 0;
}

void Scene::save_approximation(const std::string &filename)
{
  if (m_tris.empty())
    return;

  std::ofstream ofs(filename);
  if (!ofs.is_open()) {
    std::cerr << "Error: open " << filename << " failed." << std::endl;
    return;
  }

  ofs << "OFF\n" << m_anchor_pos.size() << ' ' << m_tris.size() << ' ' << "0\n";
  BOOST_FOREACH(const Point_3 &pt, m_anchor_pos)
    ofs << pt.x() << ' ' << pt.y() << ' ' << pt.z() << ' ' << '\n';
  BOOST_FOREACH(const std::vector<std::size_t> &t, m_tris)
    ofs << 3 << ' ' << t[0] << ' ' << t[1] << ' ' << t[2] << '\n';
  ofs.flush();
  ofs.close();
}

void Scene::set_metric(const int m) {
  if (m < 0 || m > 2) {
    std::cout << "Error: unknown metric." << std::endl;
    return;
  }

  for (Facet_iterator fitr = m_pmesh->facets_begin(); fitr != m_pmesh->facets_end(); ++fitr)
    m_fidx_map[fitr] = 0;

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
#endif
  m_px_color.clear();
  m_anchor_pos.clear();
  m_anchor_vtx.clear();
  m_bdrs.clear();
  m_tris.clear();

  m_view_polyhedron = true;
  m_view_wireframe = false;
  m_view_boundary = false;
  m_view_proxies = false;
  m_view_anchors = false;
  m_view_approximation = false;

  switch(m) {
    case 0: return m_approx.set_metric(Approximation_wrapper::L21);
    case 1: return m_approx.set_metric(Approximation_wrapper::L2);
    case 2: return m_approx.set_metric(Approximation_wrapper::Compact);
  }
}

void Scene::seeding(
  const CGAL::Approximation_seeding_tag method,
  const boost::optional<std::size_t> max_nb_proxies,
  const boost::optional<FT> min_error_drop,
  const std::size_t nb_relaxations,
  const std::size_t nb_iterations)
{
  if (!m_pmesh)
    return;

  m_approx.seeding(method, max_nb_proxies, min_error_drop, nb_relaxations);
  m_approx.run(nb_iterations);
  m_approx.proxy_map(m_fidx_pmap);
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
  m_approx.get_l21_proxies(std::back_inserter(m_proxies));
#endif

  // generate proxy color map
  m_px_color.clear();
  for (std::size_t i = 0; i < m_approx.proxies_size(); i++)
    m_px_color.push_back(rand_0_255());

  // update display options
  m_view_boundary = true;
}

void Scene::run_one_step()
{
  m_approx.run(1);
  m_approx.proxy_map(m_fidx_pmap);
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
  m_approx.get_l21_proxies(std::back_inserter(m_proxies));
#endif
}

void Scene::add_one_proxy()
{
  if (m_approx.add_one_proxy() == 0)
    return;

  m_approx.proxy_map(m_fidx_pmap);
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
  m_approx.get_l21_proxies(std::back_inserter(m_proxies));
#endif

  // add one proxy color
  m_px_color.push_back(rand_0_255());
}

void Scene::teleport_one_proxy()
{
  m_approx.teleport_one_proxy();
  m_approx.proxy_map(m_fidx_pmap);
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
  m_approx.get_l21_proxies(std::back_inserter(m_proxies));
#endif
}

void Scene::split(const std::size_t px_idx, const std::size_t n, const std::size_t nb_relaxations)
{
  if (m_approx.split(px_idx, n, nb_relaxations)) {
    std::cerr << "split succeeded" << std::endl;
    // add colors
    for (std::size_t i = 1; i < n; ++i)
      m_px_color.push_back(rand_0_255());
  }
  else
    std::cerr << "split failed" << std::endl;
  m_approx.proxy_map(m_fidx_pmap);
#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  m_proxies.clear();
  m_approx.get_l21_proxies(std::back_inserter(m_proxies));
#endif
}

void Scene::extract_mesh(const double chord_error,
  const bool is_relative_to_chord,
  const bool with_dihedral_angle,
  const bool if_optimize_anchor_location,
  const bool pca_plane)
{
  m_tris.clear();
  m_anchor_pos.clear();
  m_anchor_vtx.clear();
  m_bdrs.clear();

  m_approx.extract_mesh(chord_error,
    is_relative_to_chord,
    with_dihedral_angle,
    if_optimize_anchor_location,
    pca_plane);
  m_approx.indexed_triangles(std::back_inserter(m_tris));
  m_approx.anchor_points(std::back_inserter(m_anchor_pos));
  m_approx.anchor_vertices(std::back_inserter(m_anchor_vtx));
  m_approx.indexed_boundary_polygons(std::back_inserter(m_bdrs));

  // update display options
  m_view_anchors = true;
  m_view_approximation = true;
}

void Scene::draw()
{
  if (m_view_polyhedron) {
    if (m_view_wireframe || m_view_boundary) {
      ::glEnable(GL_POLYGON_OFFSET_FILL);
      ::glPolygonOffset(3.0f, 1.0f);
    }
    render_polyhedron();
  }

  if (m_view_wireframe)
    render_wireframe();
  
  if (m_view_boundary)
    render_boundary();

  if (m_view_anchors) {
    render_anchors();
    render_borders();
  }

  if (m_view_approximation)
    render_approximation();

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
  if (m_view_proxies)
    render_proxies();
#endif
}

void Scene::render_polyhedron()
{
  if (!m_pmesh)
    return;

  const std::size_t px_num = m_approx.proxies_size();
  ::glEnable(GL_LIGHTING);
  ::glColor3ub(200, 200, 200);
  ::glBegin(GL_TRIANGLES);
  for (Facet_iterator fitr = m_pmesh->facets_begin(); fitr != m_pmesh->facets_end(); ++fitr) {
    Halfedge_handle he = fitr->halfedge();
    const Point_3 &a = he->opposite()->vertex()->point();
    const Point_3 &b = he->vertex()->point();
    const Point_3 &c = he->next()->vertex()->point();

    if (px_num) {
      std::size_t cidx = m_px_color[m_fidx_pmap[fitr]];
      ::glColor3ub(ColorCheatSheet::r(cidx), ColorCheatSheet::g(cidx), ColorCheatSheet::b(cidx));
    }

    Vector_3 norm = CGAL::unit_normal(a, b, c);
    ::glNormal3d(norm.x(), norm.y(), norm.z());
    ::glVertex3d(a.x(), a.y(), a.z());
    ::glVertex3d(b.x(), b.y(), b.z());
    ::glVertex3d(c.x(), c.y(), c.z());
  }
  ::glEnd();
}

void Scene::render_wireframe()
{
  if (!m_pmesh)
    return;
  
  // draw black edges
  ::glDisable(GL_LIGHTING);
  ::glColor3ub(0, 0, 0);
  ::glLineWidth(1.0f);
  ::glBegin(GL_LINES);
  for (Edge_iterator he = m_pmesh->edges_begin();
    he != m_pmesh->edges_end(); he++) {
    const Point_3& a = he->vertex()->point();
    const Point_3& b = he->opposite()->vertex()->point();
    ::glVertex3d(a.x(),a.y(),a.z());
    ::glVertex3d(b.x(),b.y(),b.z());
  }
  ::glEnd();
}

void Scene::render_boundary()
{
  if (!m_pmesh || !m_approx.proxies_size())
    return;

  ::glDisable(GL_LIGHTING);
  ::glColor3ub(0, 0, 0);
  ::glLineWidth(1.0);
  ::glBegin(GL_LINES);
  for (Edge_iterator eitr = m_pmesh->edges_begin(); eitr != m_pmesh->edges_end(); ++eitr) {
    std::size_t segid0 = std::numeric_limits<std::size_t>::max();
    if (!eitr->is_border())
      segid0 = m_fidx_pmap[eitr->facet()];
    std::size_t segid1 = std::numeric_limits<std::size_t>::max();
    if (!eitr->opposite()->is_border())
      segid1 = m_fidx_pmap[eitr->opposite()->facet()];

    if (segid0 != segid1) {
      const Point_3 &p0 = eitr->vertex()->point();
      const Point_3 &p1 = eitr->opposite()->vertex()->point();
      ::glVertex3d(p0.x(), p0.y(), p0.z());
      ::glVertex3d(p1.x(), p1.y(), p1.z());
    }
  }
  ::glEnd();
}

void Scene::render_anchors()
{
  ::glDisable(GL_LIGHTING);
  ::glColor3ub(0, 0, 0);
  ::glPointSize(5.0f);
  ::glBegin(GL_POINTS);
  BOOST_FOREACH(const Point_3 &pt, m_anchor_pos) {
    ::glVertex3d(pt.x(), pt.y(), pt.z());
  }
  ::glEnd();

  ::glColor3ub(255, 255, 255);
  ::glPointSize(5.0f);
  ::glBegin(GL_POINTS);
  BOOST_FOREACH(const Polyhedron_3::Vertex_handle &vtx, m_anchor_vtx) {
    const Point_3 &pt = vtx->point();
    ::glVertex3d(pt.x(), pt.y(), pt.z());
  }
  ::glEnd();

  ::glLineWidth(1.0f);
  ::glColor3ub(0, 0, 255);
  ::glBegin(GL_LINES);
  for (std::size_t i = 0; i < m_anchor_pos.size(); ++i) {
    const Point_3 &ps = m_anchor_vtx[i]->point();
    ::glVertex3d(ps.x(), ps.y(), ps.z());
    const Point_3 &pt = m_anchor_pos[i];
    ::glVertex3d(pt.x(), pt.y(), pt.z());
  }
  ::glEnd();
}

void Scene::render_borders()
{
  ::glDisable(GL_LIGHTING);
  ::glLineWidth(3.0f);
  ::glColor3ub(255, 0, 0);
  BOOST_FOREACH(const std::vector<std::size_t> &b, m_bdrs) {
    ::glBegin(GL_LINE_LOOP);
    BOOST_FOREACH(const std::size_t &a, b) {
      const Point_3 &pt = m_anchor_pos[a];
      ::glVertex3d(pt.x(), pt.y(), pt.z());
    }
    ::glEnd();
  }
}

void Scene::render_approximation()
{
  ::glEnable(GL_LIGHTING);
  ::glPolygonOffset(3.0, 1.0);
  ::glLineWidth(1.0f);
  ::glColor3ub(0, 0, 255);
  BOOST_FOREACH(const std::vector<std::size_t> &t, m_tris) {
    ::glBegin(GL_LINE_LOOP);
    const Point_3 &p0 = m_anchor_pos[t[0]];
    ::glVertex3d(p0.x(), p0.y(), p0.z());
    const Point_3 &p1 = m_anchor_pos[t[1]];
    ::glVertex3d(p1.x(), p1.y(), p1.z());
    const Point_3 &p2 = m_anchor_pos[t[2]];
    ::glVertex3d(p2.x(), p2.y(), p2.z());
    ::glEnd();
  }

  ::glColor3ub(200, 200, 200);
  // ::glPolygonMode(GL_FRONT, GL_FILL);
  ::glBegin(GL_TRIANGLES);
  BOOST_FOREACH(const std::vector<std::size_t> &t, m_tris) {
    const Point_3 &p0 = m_anchor_pos[t[0]];
    const Point_3 &p1 = m_anchor_pos[t[1]];
    const Point_3 &p2 = m_anchor_pos[t[2]];
    Vector_3 n = CGAL::unit_normal(p0, p1, p2);
    ::glNormal3d(n.x(), n.y(), n.z());
    ::glVertex3d(p0.x(), p0.y(), p0.z());
    ::glVertex3d(p1.x(), p1.y(), p1.z());
    ::glVertex3d(p2.x(), p2.y(), p2.z());
  }
  ::glEnd();
}

#ifdef CGAL_SURFACE_MESH_APPROXIMATION_DEBUG
void Scene::render_proxies()
{
  ::glDisable(GL_LIGHTING);
  ::glColor3ub(255, 0, 0);
  ::glLineWidth(3.0f);
  ::glBegin(GL_LINES);
  BOOST_FOREACH(const L21_proxy_wrapper &pxw, m_proxies) {
    const Halfedge_handle he = pxw.seed->halfedge();
    const Vector_3 norm = pxw.px;
    const Point_3 cen = CGAL::centroid(he->opposite()->vertex()->point(),
      he->vertex()->point(),
      he->next()->vertex()->point());
    const Point_3 end = cen + norm;
    ::glVertex3d(cen.x(), cen.y(), cen.z());
    ::glVertex3d(end.x(), end.y(), end.z());
  }
  ::glEnd();
}
#endif
