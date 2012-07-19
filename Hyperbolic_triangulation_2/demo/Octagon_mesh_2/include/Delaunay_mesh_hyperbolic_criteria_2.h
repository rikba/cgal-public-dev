// Copyright (c) 2003-2006  INRIA Sophia-Antipolis (France).
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
// $URL: svn+ssh://scm.gforge.inria.fr/svn/cgal/branches/next/Mesh_2/include/CGAL/Delaunay_mesh_criteria_2.h $
// $Id: Delaunay_mesh_criteria_2.h 67117 2012-01-13 18:14:48Z lrineau $
// 
//
// Author(s)     : Laurent RINEAU

#ifndef CGAL_DELAUNAY_MESH_HYPERBOLIC_CRITERIA_2_H
#define CGAL_DELAUNAY_MESH_HYPERBOLIC_CRITERIA_2_H

#include <CGAL/Mesh_2/Face_badness.h>

namespace CGAL {

template <class Tr>
class Delaunay_mesh_hyperbolic_criteria_2
{
  double B;

protected:
  typedef typename Tr::Geom_traits Geom_traits;
  Geom_traits traits;

public:
  typedef typename Tr::Face_handle Face_handle;

  Delaunay_mesh_hyperbolic_criteria_2(const double bound = 0.125,
                           const Geom_traits& traits = Geom_traits())
    : B(bound), traits(traits) {}

  typedef double Quality;

  inline
  double bound() const { return B; }

  inline 
  void set_bound(const double bound) { B = bound; }

  class Is_bad
  {
  protected:
    const double B;
    const Geom_traits& traits;
  public:
    typedef typename Tr::Point Point_2;
      
    Is_bad(const double bound, const Geom_traits& traits)
      : B(bound), traits(traits) {}
      
    Mesh_2::Face_badness operator()(const Quality q) const
    {
      if( q > B )
	return Mesh_2::BAD;
      else
	return Mesh_2::NOT_BAD;
    }

    Mesh_2::Face_badness operator()(const Face_handle& fh,
				    Quality& q) const
    {
     //return the hyperbolic area of the triangle
		typedef typename Tr::Geom_traits Geom_traits;
		
		typedef typename Geom_traits::Construct_weighted_circumcenter_2 Construct_weighted_circumcenter_2;
		typedef typename Geom_traits::Weighted_point_2 Weighted_point_2;
		typedef typename Geom_traits::Bare_point Bare_point;
		typedef typename Geom_traits::Point_2 Point_2;
		typedef typename Geom_traits::Triangle_2 Triangle_2;
		typedef typename Geom_traits::Compute_area_2 Compute_area_2;
		typedef typename Geom_traits::Compute_squared_distance_2 Compute_squared_distance_2;
		typedef	typename Geom_traits::Collinear_2 Collinear_2;
		
		Geom_traits traits;
		
		Compute_area_2 area_2 = traits.compute_area_2_object();
		Compute_squared_distance_2 squared_distance = 
		traits.compute_squared_distance_2_object();
		
		double area = 0.;
		
		for (int i=0; i<3; i++) 
		{
			
			
			Point_2& pa = fh->vertex(i)->point();
			Point_2& pb = fh->vertex((i+1)%3)->point();
			Point_2& pc = fh->vertex((i+2)%3)->point();
			const Point_2&  po = Point_2(0,0);
			
			if(!Collinear_2()(pa,pb,po))
			{
				Weighted_point_2 wa(pa);
				Weighted_point_2 wb(pb);
				Weighted_point_2 wo(po, 10000);
				Bare_point center = Construct_weighted_circumcenter_2()(wa, wb, wo);
				
				double radius = Compute_squared_distance_2()(pa, center);
				Triangle_2 t = traits.construct_triangle_2_object()(pa,pb,center);
				double temp_area = CGAL::to_double(area_2(t));
				temp_area = temp_area*temp_area;
				double sin_2 = 4*temp_area/(radius*radius);
				double alpha = 0.25*acos(1.-2.*sin_2);
				double dist_c = Compute_squared_distance_2()(pc,center);
				if (dist_c > radius) {
					area = area + alpha;
				}
				else {
					area = area - alpha;
				}
			}

		}
		q=area;
		

      return operator()(q);
    }
  };

  Is_bad is_bad_object() const
  { return Is_bad(B, traits); }
};

} // end namespace CGAL

#endif
