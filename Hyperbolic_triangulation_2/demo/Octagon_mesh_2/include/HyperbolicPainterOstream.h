// Copyright (c) 2011   INRIA Sophia-Antipolis (France).
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org); you may redistribute it under
// the terms of the Q Public License version 1.0.
// See the file LICENSE.QPL distributed with CGAL.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: svn+ssh://mbogdanov@scm.gforge.inria.fr/svn/cgal/trunk/Triangulation_2/include/CGAL/Mathieu_traits.h $
// $Id: Mathieu_traits.h 57323 2010-07-05 10:07:39Z sloriot $
// 
//
// Author(s)     : Mikhail Bogdanov

#ifndef CGAL_HYPERBOLIC_PAINTER_OSTREAM_H
#define CGAL_HYPERBOLIC_PAINTER_OSTREAM_H

#include <CGAL/Qt/PainterOstream.h>

#include "Mathieu_traits.h"

namespace CGAL{
	
	namespace Qt {
		
		template <typename K>
		class PainterOstream<Mathieu_traits<K> > : public PainterOstream<K> {
		public:
			typedef PainterOstream<K> Base;
			typedef PainterOstream<Mathieu_traits<K> > Self;
			
			typedef Mathieu_traits<K> Gt;
			
		private:
			typedef typename Gt::Segment_2      Segment_2;
			typedef typename Gt::Line_Segment_2 Line_segment_2;
			typedef typename Gt::Arc_2          Arc_2;
			//typedef typename Gt::Line_2         Line_2;
			
			typedef typename K::Point_2    Point_2;
			typedef typename K::Circle_2   Circle_2;
			
		public:
			PainterOstream(QPainter* p, QRectF rect = QRectF())
			: Base(p, rect), qp(p), convert(rect) {}
			
			using Base::operator <<;
			
			PainterOstream& operator << (const Segment_2& s)
			{
				
				const Circle_2& circle = s.getcircle();
				const Point_2& center = circle.center();
				const Point_2& source = s.getarc(1);
				const Point_2& target = s.getarc(2);
				double radius = circle.squared_radius();
				
				if (radius > 100000) {
					const Line_segment_2 seg(source,target);
					operator << (seg);
					return *this;
				}
				else {
				
				double asource = std::atan2( -to_double(source.y() - center.y()),
											to_double(source.x() - center.x())); 
				double atarget = std::atan2( -to_double(target.y() - center.y()),
												to_double(target.x() - center.x()));
					
				std::swap(asource, atarget);
				double aspan = atarget - asource;
					
				if(aspan < 0.)
					aspan += 2 * CGAL_PI;
					
				const double coeff = 180*16/CGAL_PI;
				this->qp->drawArc(this->convert(circle.bbox()), 
									(int)(asource * coeff), 
									(int)(aspan * coeff));

				return *this;
				}
			}
			
		private:
			// ToDo: These objects must be deleted
			// Copies of these objects are in the base class.
			// We need access to the copies in the base class.
			QPainter* qp;
			Converter<K> convert;      
		};
		
	} //namespace Qt
	
} //namespace CGAL

#endif // CGAL_HYPERBOLIC_PAINTER_OSTREAM_H

