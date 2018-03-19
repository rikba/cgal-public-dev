// Copyright (c) 2010-2016  INRIA Sophia Antipolis, INRIA Nancy (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: 
// $Id: 
// 
//
// Author(s)     : Mikhail Bogdanov
//                 Monique Teillaud <Monique.Teillaud@inria.fr>

#ifndef CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H
#define CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H

#include <CGAL/Hyperbolic_triangulation_face_base_2.h>
#include <CGAL/Delaunay_triangulation_2.h>

#include <stack>
#include <set>

namespace CGAL {
 
template < class Gt, 
  class Tds = Triangulation_data_structure_2 <
                           Triangulation_vertex_base_2<Gt>, 
                           Hyperbolic_triangulation_face_base_2<Gt> > >
  class Hyperbolic_Delaunay_triangulation_2 
  : public Delaunay_triangulation_2<Gt,Tds>
{
public:
  typedef Hyperbolic_Delaunay_triangulation_2<Gt, Tds> Self;
  typedef Delaunay_triangulation_2<Gt,Tds> Base;
  
  typedef typename Base::size_type             size_type;
  
  typedef typename Base::Vertex_handle Vertex_handle;
  typedef typename Base::Face_handle   Face_handle;
  typedef typename Base::Edge          Edge;
  
#ifndef CGAL_CFG_USING_BASE_MEMBER_BUG_2  
  using Base::cw;
  using Base::ccw;
  using Base::geom_traits;
#endif
  
  typedef typename Base::Edge_circulator       Edge_circulator;
  typedef typename Base::Face_circulator       Face_circulator;
  typedef typename Base::Vertex_circulator     Vertex_circulator;
  
  typedef typename Base::All_vertices_iterator    All_vertices_iterator;
  typedef typename Base::All_edges_iterator       All_edges_iterator;
  typedef typename Base::All_faces_iterator       All_faces_iterator;
 
  typedef Gt Geom_traits;
  typedef typename Geom_traits::FT                    FT;
  typedef typename Geom_traits::Point_2               Point;
  typedef typename Geom_traits::Voronoi_point_2       Voronoi_point;
  typedef typename Geom_traits::Hyperbolic_segment_2  Hyperbolic_segment;
  
  typedef typename Base::Locate_type                  Locate_type;

  typedef typename Geom_traits::Side_of_hyperbolic_triangle_2 Side_of_hyperbolic_triangle;
  typedef typename Geom_traits::Is_hyperbolic         Is_hyperbolic;

  Hyperbolic_Delaunay_triangulation_2(const Gt& gt = Gt())
  : Delaunay_triangulation_2<Gt,Tds>(gt) {}
  
  Hyperbolic_Delaunay_triangulation_2(
	     const Hyperbolic_Delaunay_triangulation_2<Gt,Tds> &tr)
       : Delaunay_triangulation_2<Gt,Tds>(tr)
  {   CGAL_triangulation_postcondition( this->is_valid() );}

  
  void mark_star(Vertex_handle v) const
  {
    if(!is_star_bounded(v)) {
      mark_star_faces(v);
    }
  }
  
  Vertex_handle insert(const Point  &p, 
                       Face_handle start = Face_handle() )
  {
    Vertex_handle v = Base::insert(p, start);
    mark_star(v);
    
    return v;
  }
  
  Vertex_handle insert(const Point& p,
                       typename Base::Locate_type lt,
                       Face_handle loc, int li )
  {
    Vertex_handle v = Base::insert(p, lt, loc, li);
    mark_star(v);    
    
    return v;
  }
    
#ifndef CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO
  template < class InputIterator >
  std::ptrdiff_t
  insert( InputIterator first, InputIterator last,
         typename boost::enable_if<
         boost::is_base_of<
         Point,
         typename std::iterator_traits<InputIterator>::value_type
         >
         >::type* = NULL
         )
#else
  template < class InputIterator >
  std::ptrdiff_t
  insert(InputIterator first, InputIterator last)
#endif //CGAL_TRIANGULATION_2_DONT_INSERT_RANGE_OF_POINTS_WITH_INFO 
  {
    size_type n = Base::insert(first, last);
    
    mark_finite_non_hyperbolic_faces();
    
    return n;
  }
  
  bool is_infinite(Vertex_handle v) const
  {
    return Base::is_infinite(v);
  }

  bool is_non_hyperbolic(Face_handle f) const
  {
    return has_infinite_vertex(f) || is_finite_non_hyperbolic(f);
  }
  
  bool is_non_hyperbolic(Face_handle f, int i) const 
  {
    return has_infinite_vertex(f, i) || is_finite_non_hyperbolic(f, i);
  }
  
  bool is_non_hyperbolic(const Edge& e) const 
  {
    return is_non_hyperbolic(e.first, e.second);
  }
  
  bool is_non_hyperbolic(const Edge_circulator& ec) const 
  {
    return is_non_hyperbolic(*ec);
  }
  
  bool is_non_hyperbolic(const All_edges_iterator& ei) const 
  {
    return is_non_hyperbolic(*ei);
  }
  
  // is_infinite functions are kept in order to reuse Triangulation_2 demo :
  //              apply_to_range is called by Qt/TriangulationGraphicsItem.h  
  // TODO: document that is_infinite functions are not inherited from Triangulation_2
  bool is_infinite(Face_handle f) const { return is_non_hyperbolic(f); }
  bool is_infinite(Face_handle f, int i) const { return is_non_hyperbolic(f,i); }
  bool is_infinite(const Edge e) const { return is_non_hyperbolic(e); }
  bool is_infinite(const Edge_circulator& ec) const { return is_non_hyperbolic(ec); }
  bool is_infinite(const All_edges_iterator& ei) const { return is_non_hyperbolic(ei); }

private:
  
  bool has_infinite_vertex(Face_handle f) const
  {
    return Base::is_infinite(f);
  }
  
  bool has_infinite_vertex(Face_handle f, int i) const
  {
    return Base::is_infinite(f, i);
  }
  
  bool has_infinite_vertex(const Edge& e) const
  {
    return Base::is_infinite(e);
  }
  
  int get_finite_non_hyperbolic_edge(Face_handle f) const
  {
    assert(is_finite_non_hyperbolic(f));
    
    return f->get_non_hyperbolic_edge(); 
  }
  
  bool is_finite_non_hyperbolic(Face_handle f) const
  {
    return f->is_finite_non_hyperbolic();
  }
  
  bool is_finite_non_hyperbolic(Face_handle f, int i) const
  {
    if(this->dimension() <= 1) {
      return false;
    }
    
    if(is_finite_non_hyperbolic(f) && get_finite_non_hyperbolic_edge(f) == i) {
      return true;
    }
    
    // another incident face and corresponding index
    Face_handle f2 = f->neighbor(i);
    int i2 =  f2->index(f);
    
    if(is_finite_non_hyperbolic(f2) && get_finite_non_hyperbolic_edge(f2) == i2) {
      return true;
    }
    
    return false;
  }
  
  bool is_finite_non_hyperbolic(const Edge& e) const
  {
    return is_finite_non_hyperbolic(e.first, e.second);
  }
  
  // Depth-first search (dfs) and marking the finite non_hyperbolic faces.
  void mark_finite_non_hyperbolic_faces() const
  {
    if(this->dimension() <= 1) return;
      
    std::set<Face_handle> visited_faces;
    
    // maintain a stack to be able to backtrack
    // to the most recent faces which neighbors are not visited
    std::stack<Face_handle> backtrack;
    
    // start from a face with infinite vertex
    Face_handle current = Base::infinite_face();
    
    // mark it as visited
    visited_faces.insert(current);
    
    // put the element whose neighbors we are going to explore.
    backtrack.push(current);
    
    // test whether a face is finite non_hyperbolic or not
    Mark_face test(*this);
    
    Face_handle next;
    
    while(!backtrack.empty()) {
      // take a face
      current = backtrack.top();
      
      // start visiting the neighbors
      int i = 0;
      for(; i < 3; i++) {
        next = current->neighbor(i);
        
        // if a neighbor is already visited, then stop going deeper
        if(visited_faces.find(next) != visited_faces.end()) {
          continue;
        }
        
        visited_faces.insert(next);
        mark_face(next, test);
        
        // go deeper if the neighbor is non_hyperbolic
        if(is_non_hyperbolic(next)) {
          backtrack.push(next);
          break;
        }
      }
      
      // if all the neighbors are already visited, then remove "current" face.
      if(i == 3) {
        backtrack.pop();
      }
    }
    
  }
  
  // check if a star is bounded by finite faces
  // TODO: rename this function name
  bool is_star_bounded(Vertex_handle v) const
  {
    if(this->dimension() <= 1) {
      return true;
    }
    
    Face_handle f = v->face();
    Face_handle next;
    int i;
    Face_handle start(f);
    
    Face_handle opposite_face;
    
    do {
      i = f->index(v);
      next = f->neighbor(ccw(i));  // turn ccw around v
      
      opposite_face = f->neighbor(i);
      if(this->is_non_hyperbolic(opposite_face)) {
        return false;
      }
      
      f = next;
    } while(next != start);
    
    return true;
  }
  
  //use the function: insert_and_give_new_faces?
  
  void mark_star_faces(Vertex_handle v) const
  {
    // TODO: think of it
    if(this->dimension() <= 1) return;
    
    Mark_face test(*this);
    
    Face_handle f = v->face();
    Face_handle start(f), next;
    int i;
    do {
      i = f->index(v);
      next = f->neighbor(ccw(i));  // turn ccw around v
      
      mark_face(f, test);
      
      f = next;
    } while(next != start);
    return;
  }
  
  template<class Mark_face_test>
  void mark_face(const Face_handle& f, const Mark_face_test& test) const
  {
    f->set_finite_non_hyperbolic(test(f));
  }
  
  void mark_face(const Face_handle& f) const
  {
    Mark_face test(*this);
    mark_face(f, test);
  }
    
  class Mark_face
  {
  public:
    Mark_face(const Self& tr) :
      _tr(tr)
    {}
    
    bool operator ()(const Face_handle& f) const
    {
      typedef typename Gt::Is_hyperbolic Is_hyperbolic;
      
      if(_tr.has_infinite_vertex(f)) {
	return false; 
      }
      
      Point p0 = f->vertex(0)->point();
      Point p1 = f->vertex(1)->point();
      Point p2 = f->vertex(2)->point();
      int ind = 0;
      
      Is_hyperbolic is_hyperbolic = _tr.geom_traits().is_hyperbolic_object();
      if(is_hyperbolic(p0, p1, p2, ind) == false) {
	f->set_finite_non_hyperbolic(true); // MT should not be necessary, return true should be enough (?)
	f->set_non_hyperbolic_edge(ind);
	return true; 
      }
      
      // the face is finite and hyperbolic
      return false; 
    }
    
  private:
  
    Mark_face(const Mark_face&);
    Mark_face& operator= (const Mark_face&);
    
    const Self& _tr;
  }; 
  
public:
  // This class is used to generate the Finite_*_iterators.
  class Non_hyperbolic_tester
  {
    const Self *t;
  public:
    Non_hyperbolic_tester() {}
    Non_hyperbolic_tester(const Self *tr)	  : t(tr) {}
    
    bool operator()(const All_vertices_iterator & vit) const  {
      return t->is_infinite(vit);
    }
    bool operator()(const All_faces_iterator & fit) const {
      return t->is_non_hyperbolic(fit);
    }
    bool operator()(const All_edges_iterator & eit ) const {
      return t->is_non_hyperbolic(eit);
    }
  };
  
  Non_hyperbolic_tester
  non_hyperbolic_tester() const
  {
    return Non_hyperbolic_tester(this);
  }
  
  class Hyperbolic_faces_iterator
  : public Filter_iterator<All_faces_iterator, Non_hyperbolic_tester> 
  {
    typedef Filter_iterator<All_faces_iterator, Non_hyperbolic_tester> Base;
    typedef Hyperbolic_faces_iterator                           Self;
  public:
    Hyperbolic_faces_iterator() : Base() {}
    Hyperbolic_faces_iterator(const Base &b) : Base(b) {}
    Self & operator++() { Base::operator++(); return *this; }
    Self & operator--() { Base::operator--(); return *this; }
    Self operator++(int) { Self tmp(*this); ++(*this); return tmp; }
    Self operator--(int) { Self tmp(*this); --(*this); return tmp; }
    operator const Face_handle() const { return Base::base(); }
  };

  Hyperbolic_faces_iterator
  hyperbolic_faces_begin() const
  {
    if ( this->dimension() < 2 )
      return hyperbolic_faces_end();
    return CGAL::filter_iterator(this->all_faces_end(),
                                 Non_hyperbolic_tester(this),
                                 this->all_faces_begin() );
  } 

  Hyperbolic_faces_iterator
  hyperbolic_faces_end() const
  {
    return CGAL::filter_iterator(this->all_faces_end(),
                                 Non_hyperbolic_tester(this)   );
  }

  typedef Filter_iterator<All_edges_iterator, Non_hyperbolic_tester> Hyperbolic_edges_iterator;
  
  Hyperbolic_edges_iterator
  hyperbolic_edges_begin() const
  {
    if ( this->dimension() < 1 )
      return hyperbolic_edges_end();
    return CGAL::filter_iterator(this->all_edges_end(),
                                 Non_hyperbolic_tester(this),
                                 this->all_edges_begin());
  }
  
  Hyperbolic_edges_iterator
  hyperbolic_edges_end() const
  {
    return CGAL::filter_iterator(this->all_edges_end(),
                                 Non_hyperbolic_tester(this) );
  }

  size_type number_of_hyperbolic_faces() const
  {
    return std::distance(hyperbolic_faces_begin(), hyperbolic_faces_end());
  }

  size_type number_of_hyperbolic_edges() const
  {
    return std::distance(hyperbolic_edges_begin(), hyperbolic_edges_end());
  }

  // Finite faces/edges iterators kept for the demo in order to reuse Triangulation_2 demo (see above)
  // TODO: document that they are not inherited from Triangulation_2
  typedef Hyperbolic_faces_iterator Finite_faces_iterator;
  Finite_faces_iterator finite_faces_begin() const { return hyperbolic_faces_begin(); }
  Finite_faces_iterator finite_faces_end() const { return hyperbolic_faces_end(); }
  typedef Hyperbolic_edges_iterator Finite_edges_iterator;
  Finite_edges_iterator finite_edges_begin() const { return hyperbolic_edges_begin(); }
  Finite_edges_iterator finite_edges_end() const { return hyperbolic_edges_end(); }
  
  Voronoi_point
  dual(Face_handle f) const
  {
    CGAL_triangulation_precondition (!this->is_non_hyperbolic(f));
    
    return this->geom_traits().construct_hyperbolic_circumcenter_2_object()
         ( f->vertex(0)->point(), f->vertex(1)->point(), f->vertex(2)->point());
  }

  Hyperbolic_segment
  dual(const Edge& e) const
  { 
    return dual(e.first, e.second);
  }

  Hyperbolic_segment
  dual(Face_handle f, int i) const
  {
    CGAL_triangulation_precondition (!this->is_non_hyperbolic(f,i));
    
    if(this->dimension() == 1) {
      Point p = f->vertex(cw(i))->point();
      Point q = f->vertex(ccw(i))->point();
      
      // hyperbolic line
      Hyperbolic_segment line = this->geom_traits().construct_hyperbolic_bisector_2_object()(p,q);
      return line;
    }
    
    Face_handle n = f->neighbor(i);
    int in = n->index(f);
    //TODO MT store values of bools to avoid recomputing is-hyperbolic several times

    // boths faces are non_hyperbolic, but the incident edge is hyperbolic
    if( is_non_hyperbolic(f) && is_non_hyperbolic(n) ){
      const Point& p = f->vertex(ccw(i))->point();
      const Point& q = f->vertex(cw(i))->point();
      
      // hyperbolic line
      Hyperbolic_segment line = 
          this->geom_traits().construct_hyperbolic_bisector_2_object()(p,q);
      return line;
    }
    
    // both faces are hyperbolic
    if( !is_non_hyperbolic(f) && !is_non_hyperbolic(n) ) {
      const Point& p = f->vertex(ccw(i))->point();
      const Point& q = f->vertex(cw(i))->point();
   
      Hyperbolic_segment s = 
      this->geom_traits().construct_hyperbolic_bisector_2_object()
      (p,q,f->vertex(i)->point(),n->vertex(in)->point());
      //TODO MT cut edge at dual points !!!!
      return s;
    }
    
    // one of the incident faces is non_hyperbolic
    Face_handle hyp_face = f;
    
    if(is_non_hyperbolic(f)) {
      hyp_face = n;
      i = in;
    }
    
    const Point& p = hyp_face->vertex(ccw(i))->point();
    const Point& q = hyp_face->vertex(cw(i))->point();
    
    // ToDo: Line or Segment?
    // hyperbolic line and ray
    Hyperbolic_segment ray = this->geom_traits().construct_hyperbolic_bisector_2_object()(p,q,hyp_face->vertex(i)->point());
    // TODO MT cut edge at dual point !!!
    //    Segment ray = this->geom_traits().construct_ray_2_object()(dual(finite_face), line);
    return ray;
  }

public:
  Face_handle locate(const Point& p, const Face_handle hint = Face_handle()) const {
    Locate_type lt;
    int li;
    return locate(p, lt, li, hint);
  }

  Face_handle locate(const Point& query, Locate_type& lt, int &li, Face_handle hint = Face_handle()) const {
    
    // Perform an Euclidean location first and get close to the hyperbolic face containing the query point
    Face_handle fh = Base::locate(query, lt, li, hint);
    
    if (lt == Base::OUTSIDE_CONVEX_HULL ||
        lt == Base::OUTSIDE_AFFINE_HULL ||
        lt == Base::VERTEX) {
      return Face_handle();
    }
    
    // This case corresponds to when the point is located on an Euclidean edge.
    if (lt == Base::EDGE) {
      Point p = fh->vertex(0)->point();
      Point q = fh->vertex(1)->point();
      Point r = fh->vertex(2)->point();
      if (Is_hyperbolic()(p, q, r)) {
        Bounded_side side = Side_of_hyperbolic_triangle()(p, q, r, query, li);
        if (side == ON_BOUNDARY) {
          lt = Base::EDGE;
          return fh;
        } else {
          if (side == ON_BOUNDED_SIDE) {
            lt = Base::FACE;
            return fh;
          } else {
            // do nothing -- we still have to check the neighboring face
          }
        }
      }

      p = fh->vertex(ccw(li))->point();
      q = fh->mirror_vertex(li)->point();
      r = fh->vertex(cw(li))->point();
      if (Is_hyperbolic()(p, q, r)) {
        Bounded_side side = Side_of_hyperbolic_triangle()(p, q, r, query, li);
        if (side == ON_BOUNDARY) {
          lt = Base::EDGE;
          return fh;
        } else {
          if (side == ON_BOUNDED_SIDE) {
            lt = Base::FACE;
            return fh;
          } else {
            // There is nothing to be done now -- the point is outside the convex hull of the triangulation
            lt = Base::OUTSIDE_CONVEX_HULL;
            return Face_handle();
          }
        }
      }
    }

    // Here, the face has been located in the Euclidean face lh
    Point p = fh->vertex(0)->point();
    Point q = fh->vertex(1)->point();
    Point r = fh->vertex(2)->point();
    if (!Is_hyperbolic()(p, q, r)) {
      lt = Base::OUTSIDE_CONVEX_HULL;
      return Face_handle();
    }

    Bounded_side side = Side_of_hyperbolic_triangle()(p, q, r, query, li);
    if (side == ON_BOUNDED_SIDE) {
      lt = Base::FACE;
      return fh;
    } else {
      if (side == ON_BOUNDARY) {
        lt = Base::EDGE;
        return fh;
      } else {
        // Here, the point lies in a face that is a neighbor to fh
        for (int i = 0; i < 3; i++) {
          Face_handle nfh = fh->neighbor(i);
          if (Is_hyperbolic()(nfh->vertex(0)->point(),nfh->vertex(1)->point(),nfh->vertex(2)->point())) {
            Bounded_side nside = Side_of_hyperbolic_triangle()(nfh->vertex(0)->point(),nfh->vertex(1)->point(),nfh->vertex(2)->point(), query, li);
            if (nside == ON_BOUNDED_SIDE) {
              lt = Base::FACE;
              return nfh;
            } else if (nside == ON_BOUNDARY) {
              lt = Base::EDGE;
              return nfh;
            }
          }
        }

        // At this point, the point lies outside of the convex hull of the triangulation,
        // since it has not been found in any of the hyperbolic faces adjacent to fh.
        lt = Base::OUTSIDE_CONVEX_HULL;
        return Face_handle();
      }
    }

    // We never reach this point, but we have to make the compiler happy
    lt = Base::OUTSIDE_CONVEX_HULL;
    return Face_handle();
  }

};
  
} //namespace CGAL

#endif // CGAL_HYPERBOLIC_DELAUNAY_TRIANGULATION_2_H
