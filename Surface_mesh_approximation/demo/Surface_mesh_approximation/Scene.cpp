#include "Scene.h"

#include <iostream>
#include <fstream>

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
  if(m_pmesh == NULL) {
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

  if(!in || !fileinfo.isFile() || ! fileinfo.isReadable()) {
    std::cerr << "unable to open file" << std::endl;
    QApplication::restoreOverrideCursor();
    return -1;
  }

  if(m_pmesh != NULL)
    delete m_pmesh;

  // allocate new polyhedron
  m_pmesh = new Polyhedron_3;
  in >> *m_pmesh;
  if(!in) {
    std::cerr << "invalid OFF file" << std::endl;
    QApplication::restoreOverrideCursor();

    delete m_pmesh;
    m_pmesh = NULL;
    return -1;
  }

  m_fidx_map.clear();
  for(Facet_iterator fitr = m_pmesh->facets_begin(); fitr != m_pmesh->facets_end(); ++fitr)
    m_fidx_map.insert(std::pair<Facet_handle, std::size_t>(fitr, 0));

  m_vsa.set_metric(VSA::L21);
  m_vsa.set_mesh(*m_pmesh);

  m_tris.clear();
  m_anchor_pos.clear();
  m_anchor_vtx.clear();
  m_bdrs.clear();

  m_view_polyhedron = true;
  m_view_wireframe = false;
  m_view_boundary = false;
  m_view_anchors = false;
  m_view_approximation = false;

  QApplication::restoreOverrideCursor();
  return 0;
}

void Scene::save_approximation(const std::string &filename)
{
  if(m_tris.empty())
    return;

  std::ofstream ofs(filename);
  if(!ofs.is_open()) {
    std::cerr << "Error: open " << filename << " failed." << std::endl;
    return;
  }

  ofs << "OFF\n" << m_anchor_pos.size() << ' ' << m_tris.size() / 3 << ' ' << "0\n";
  BOOST_FOREACH(const Point_3 &pt, m_anchor_pos)
    ofs << pt.x() << ' ' << pt.y() << ' ' << pt.z() << ' ' << '\n';
  for(std::vector<int>::iterator titr = m_tris.begin(); titr != m_tris.end(); titr += 3)
    ofs << 3 << ' ' << *titr << ' ' << *(titr + 1) << ' ' << *(titr + 2) << '\n';
  ofs.flush();
  ofs.close();
}

void Scene::set_metric(const int m) {
  if (m < 0 || m > 2) {
    std::cout << "Error: unknown metric." << std::endl;
    return;
  }

  switch(m) {
    case 0: return m_vsa.set_metric(VSA::L21);
    case 1: return m_vsa.set_metric(VSA::L2);
    case 2: return m_vsa.set_metric(VSA::Compact);
  }
}

void Scene::seeding_by_number(
  const int init,
  const std::size_t num_proxies,
  const std::size_t num_iterations)
{
  if(!m_pmesh)
    return;
  m_px_num = m_vsa.seeding_by_number(init, num_proxies, 5);
  for (std::size_t i = 0; i < num_iterations; ++i)
    m_vsa.run_one_step();
  m_vsa.get_proxy_map(m_fidx_pmap);
  m_view_boundary = true;
}

void Scene::run_one_step()
{
  m_vsa.run_one_step();
  m_vsa.get_proxy_map(m_fidx_pmap);
}

void Scene::add_one_proxy()
{
  m_vsa.add_one_proxy();
  m_vsa.get_proxy_map(m_fidx_pmap);
}

void Scene::teleport_one_proxy()
{
  m_vsa.teleport_one_proxy();
  m_vsa.get_proxy_map(m_fidx_pmap);
}

void Scene::meshing()
{
  Polyhedron_3 out_mesh;
  m_tris.clear();
  m_anchor_pos.clear();
  m_anchor_vtx.clear();

  m_vsa.meshing(out_mesh);
  m_vsa.get_indexed_triangles(std::back_inserter(m_tris));
  m_vsa.get_anchor_points(std::back_inserter(m_anchor_pos));
  m_vsa.get_anchor_vertices(std::back_inserter(m_anchor_vtx));
  m_bdrs = m_vsa.get_indexed_boundary_polygons();
}

void Scene::draw()
{
  if (m_view_polyhedron) {
    if(m_view_wireframe || m_view_boundary) {
      ::glEnable(GL_POLYGON_OFFSET_FILL);
      ::glPolygonOffset(3.0f, 1.0f);
    }
    ::glEnable(GL_LIGHTING);
    render_polyhedron();
  }

  if(m_view_wireframe)
    render_wireframe();
  
  if(m_view_boundary)
    render_boundary();

  if (m_view_anchors) {
    render_anchors();
    render_borders();
  }

  if (m_view_approximation)
    render_approximation();
}

void Scene::render_polyhedron()
{
  if(!m_pmesh)
    return;

  ::glColor3ub(200, 200, 200);
  ::glBegin(GL_TRIANGLES);
  for(Facet_iterator fitr = m_pmesh->facets_begin();
    fitr != m_pmesh->facets_end(); ++fitr) {
    Halfedge_around_facet_circulator he = fitr->facet_begin();
    const Point_3 &a = he->vertex()->point();
    const Point_3 &b = he->next()->vertex()->point();
    const Point_3 &c = he->prev()->vertex()->point();

    //Vector_3 norm = CGAL::normal(a, b, c);
    Vector_3 norm = CGAL::unit_normal(a, b, c);
    ::glNormal3d(norm.x(), norm.y(), norm.z());

    if(m_px_num) {
      std::size_t cidx = std::floor(static_cast<double>(m_fidx_pmap[fitr]) / static_cast<double>(m_px_num) * 256.0);
      ::glColor3ub(ColorCheatSheet::r(cidx), ColorCheatSheet::g(cidx), ColorCheatSheet::b(cidx));
    }

    ::glVertex3d(a.x(), a.y(), a.z());
    ::glVertex3d(b.x(), b.y(), b.z());
    ::glVertex3d(c.x(), c.y(), c.z());
  }
  ::glEnd();
}

void Scene::render_wireframe()
{
  if(!m_pmesh)
    return;
  
  // draw black edges
  ::glDisable(GL_LIGHTING);
  ::glColor3ub(0, 0, 0);
  ::glLineWidth(1.0f);
  ::glBegin(GL_LINES);
  for(Edge_iterator he = m_pmesh->edges_begin();
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
  if(!m_pmesh || !m_px_num)
    return;

  ::glDisable(GL_LIGHTING);
  ::glColor3ub(0, 0, 0);
  ::glLineWidth(1.0);
  ::glBegin(GL_LINES);
  for(Edge_iterator eitr = m_pmesh->edges_begin();
    eitr != m_pmesh->edges_end(); ++eitr) {
    std::size_t segid0 = std::numeric_limits<std::size_t>::max();
    if(!eitr->is_border())
      segid0 = m_fidx_pmap[eitr->facet()];
    std::size_t segid1 = std::numeric_limits<std::size_t>::max();
    if(!eitr->opposite()->is_border())
      segid1 = m_fidx_pmap[eitr->opposite()->facet()];

    if(segid0 != segid1) {
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
  for (std::vector<std::vector<std::size_t> >::iterator bitr = m_bdrs.begin(); bitr != m_bdrs.end(); ++bitr) {
    ::glBegin(GL_LINE_LOOP);
    for (std::vector<std::size_t>::iterator aitr = bitr->begin(); aitr != bitr->end(); ++aitr) {
      const Point_3 &pt = m_anchor_pos[*aitr];
      ::glVertex3d(pt.x(), pt.y(), pt.z());
    }
    ::glEnd();
  }
}

void Scene::render_approximation()
{
  ::glEnable(GL_LIGHTING);
  // ::glDisable(GL_LIGHTING);
  ::glPolygonOffset(3.0, 1.0);
  ::glLineWidth(1.0f);
  ::glColor3ub(0, 0, 255);
  for (std::vector<int>::iterator vitr = m_tris.begin(); vitr != m_tris.end(); vitr += 3) {
    ::glBegin(GL_LINE_LOOP);
    const Point_3 &p0 = m_anchor_pos[*vitr];
    ::glVertex3d(p0.x(), p0.y(), p0.z());
    const Point_3 &p1 = m_anchor_pos[*(vitr + 1)];
    ::glVertex3d(p1.x(), p1.y(), p1.z());
    const Point_3 &p2 = m_anchor_pos[*(vitr + 2)];
    ::glVertex3d(p2.x(), p2.y(), p2.z());
    ::glEnd();
  }

  ::glColor3ub(200, 200, 200);
  // ::glPolygonMode(GL_FRONT, GL_FILL);
  ::glBegin(GL_TRIANGLES);
  for (std::vector<int>::iterator vitr = m_tris.begin(); vitr != m_tris.end(); vitr += 3) {
    const Point_3 &p0 = m_anchor_pos[*vitr];
    const Point_3 &p1 = m_anchor_pos[*(vitr + 1)];
    const Point_3 &p2 = m_anchor_pos[*(vitr + 2)];
    Vector_3 n = CGAL::unit_normal(p0, p1, p2);
    ::glNormal3d(n.x(), n.y(), n.z());
    ::glVertex3d(p0.x(), p0.y(), p0.z());
    ::glVertex3d(p1.x(), p1.y(), p1.z());
    ::glVertex3d(p2.x(), p2.y(), p2.z());
  }
  ::glEnd();
}
