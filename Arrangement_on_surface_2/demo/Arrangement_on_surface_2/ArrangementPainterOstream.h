// Copyright (c) 2012  Tel-Aviv University (Israel).
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
// Author(s)     : Alex Tsui <alextsui05@gmail.com>

#ifndef CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H
#define CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H

#include <QRectF>
#include <vector>

// TODO: should be included in PainterOstream.h
#include <CGAL/Kernel/global_functions.h>
#include <CGAL/Qt/PainterOstream.h>
#include <CGAL/Qt/Converter.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arr_polyline_traits_2.h>
#include <CGAL/Arr_conic_traits_2.h>
#include <CGAL/Arr_linear_traits_2.h>
#include <CGAL/Arr_algebraic_segment_traits_2.h>

#include "Utils.h"

#include <CGAL/Curved_kernel_via_analysis_2/Curve_renderer_facade.h>

class QPainter;

namespace CGAL {
namespace Qt {

template < typename ArrTraits >
class ArrangementPainterOstreamBase : public QGraphicsSceneMixin
{
public:
  // typedefs
  typedef ArrTraits Traits;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Kernel::Point_2                      Point_2;
  typedef typename Kernel::Segment_2                    Segment_2;
  typedef typename Kernel::Ray_2                        Ray_2;
  typedef typename Kernel::Line_2                       Line_2;
  typedef typename Kernel::Triangle_2                   Triangle_2;
  typedef typename Kernel::Iso_rectangle_2              Iso_rectangle_2;
  typedef typename Kernel::Circle_2                     Circle_2;

public:
  /*! Constructor */
  ArrangementPainterOstreamBase( QPainter* p,
                                 QRectF clippingRectangle = QRectF( ) ) :
    painterOstream( p, clippingRectangle ),
    qp( p ),
    convert( clippingRectangle ),
    // scene( NULL ),
    clippingRect( QRectF( ) ), // null rectangle
    scale( 1.0 )
  {
    if ( p != 0 )
    {
      this->scale = p->worldTransform( ).m11( );
    }
  }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstreamBase() {}

  // methods
  template < typename T >
  ArrangementPainterOstreamBase& operator<<( const T& t )
  {
    this->painterOstream << t;
    return *this;
  }

  void setScene( QGraphicsScene* scene_ )
  {
    this->scene = scene_;

    // set the clipping rectangle
    if ( scene_ == NULL )
    {
      return;
    }

    this->clippingRect = this->viewportRect( );
    this->convert = Converter< Kernel >( this->clippingRect );
  }

#if 0
  void setScene( QGraphicsScene* scene_ )
  {
    this->scene = scene_;

    // set the clipping rectangle
    if ( scene_ == NULL )
    {
      return;
    }
    this->clippingRect = this->getViewportRect( );
  }
#endif

protected: // methods
#if 0
  QRectF getViewportRect( ) const
  {
    // assumes scene is not null and attached to exactly one view
    QGraphicsView* view = this->scene->views( ).first( );
    QPointF p1 = view->mapToScene( 0, 0 );
    QPointF p2 = view->mapToScene( view->width( ), view->height( ) );
    QRectF clipRect = QRectF( p1, p2 );

    return clipRect;
  }
#endif

protected:
  // fields
  PainterOstream< Kernel > painterOstream;
  QPainter* qp;
  Converter< Kernel > convert;
  // QGraphicsScene* scene;
  QRectF clippingRect;
  double scale;

}; // class ArrangementPainterOstreamBase

template < typename ArrTraits >
class ArrangementPainterOstream :
    public ArrangementPainterOstreamBase< ArrTraits >
{
public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()):
    ArrangementPainterOstreamBase< ArrTraits >( p, clippingRectangle )
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}
};

template < typename Kernel_ >
class ArrangementPainterOstream<CGAL::Arr_segment_traits_2< Kernel_> >:
  public ArrangementPainterOstreamBase<CGAL::Arr_segment_traits_2<Kernel_> >
{
public: // typedefs
  typedef Kernel_ Kernel;
  typedef CGAL::Arr_segment_traits_2< Kernel > Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()) :
    Superclass(p, clippingRectangle)
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
  {
    std::cout<<"In operator<< Arr_segment_traits_2 X_monotone_curve_2"<<std::endl;
    const Point_2& p1 = curve.source( );
    const Point_2& p2 = curve.target( );
    Segment_2 seg( p1, p2 );

    // skip segments outside our view
    QRectF seg_bb = this->convert( seg.bbox( ) );
    if ( this->clippingRect.isValid( ) &&
         ! this->clippingRect.intersects( seg_bb ) )
    {
      return *this;
    }

    this->painterOstream << seg;
    return *this;
  }

  ArrangementPainterOstream& operator<<( const Point_2& p )
  {
    std::cout<<"In operator<< Arr_segment_traits_2 Point_2"<<std::endl;

    QPointF qpt = this->convert( p );
    // clip the point if possible
    if ( this->clippingRect.isValid( ) &&
         ! this->clippingRect.contains( qpt ) )
    {
      return *this;
    }

    QPen savePen = this->qp->pen( );
    this->qp->setBrush( QBrush( savePen.color( ) ) );
    double radius = savePen.width( ) / 2.0;
    radius /= this->scale;

    this->qp->drawEllipse( qpt, radius, radius );

    this->qp->setBrush( QBrush( ) );
    this->qp->setPen( savePen );
    return *this;
  }

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }
};

template < typename SegmentTraits >
class ArrangementPainterOstream<CGAL::Arr_polyline_traits_2<SegmentTraits> > :
  public ArrangementPainterOstreamBase<CGAL::Arr_polyline_traits_2<
                                         SegmentTraits> >
{
public: // typedefs
  typedef ArrangementPainterOstreamBase<CGAL::Arr_polyline_traits_2<
                                          SegmentTraits> > Superclass;
  typedef typename Superclass::Traits                   Traits;
  typedef typename Superclass::Kernel                   Kernel;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()):
    Superclass( p, clippingRectangle )
  { }

    /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
  {
    std::cout<<"In operator<< Arr_polyline_traits_2 X_monotone_curve_2"<<std::endl;

    int cnt = 0;
    for (typename X_monotone_curve_2::Subcurve_const_iterator it =
           curve.subcurves_begin();
         it != curve.subcurves_end(); ++it)
    {
        cnt++;
        this->painterOstream << *it;
    }

    std::cout<<"cnt: "<< cnt<<std::endl;
    // TODO: implement polyline painting
#if 0
    const Point_2& p1 = curve.source( );
    const Point_2& p2 = curve.target( );
    Segment_2 seg( p1, p2 );
    this->painterOstream << seg;
#endif
    return *this;
  }

  // cloned from segtraits painter
  ArrangementPainterOstream& operator<<( const Point_2& p )
  {
    std::cout<<"In operator<< Arr_polyline_traits_2 Point_2"<<std::endl;

    QPointF qpt = this->convert( p );
    QPen savePen = this->qp->pen( );
    this->qp->setBrush( QBrush( savePen.color( ) ) );
    double radius = savePen.width( ) / 2.0;
    radius /= this->scale;

    // Draw a circle as a blue dot
    this->qp->drawEllipse( qpt, radius, radius );

    this->qp->setBrush( QBrush( ) );
    this->qp->setPen( savePen );
    return *this;
  }

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }
};

template < typename RatKernel, class AlgKernel, class NtTraits >
class ArrangementPainterOstream<CGAL::Arr_conic_traits_2<RatKernel, AlgKernel,
                                                         NtTraits > >:
  public ArrangementPainterOstreamBase< CGAL::Arr_conic_traits_2<RatKernel,
                                                                 AlgKernel,
                                                                 NtTraits> >
{
public: // typedefs
  typedef CGAL::Arr_conic_traits_2< RatKernel, AlgKernel, NtTraits > Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Construct_x_monotone_curve_2
    Construct_x_monotone_curve_2;
  typedef typename Traits::Point_2                      Intersection_point_2;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef typename Traits::Multiplicity                 Multiplicity;
  typedef typename ArrTraitsAdaptor< Traits >::Kernel   Kernel;
  typedef typename Kernel::FT                           FT;

public: // inner classes
  // utility class to use with std::sort on an Intersect_2 result set.
  class Compare_intersection_point_result
  {
  public:
    typedef std::pair< Intersection_point_2, Multiplicity > Result;
    // returns whether the point1 < point2, using x-coord to compare
    bool operator()( const Result& o1, const Result& o2 )
    {
      Point_2 p1 = o1.first;
      Point_2 p2 = o2.first;
      return ( p1.x( ) < p2.x( ) );
    }
  };

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()) :
    Superclass( p, clippingRectangle ),
    //intersect_2( this->traits.intersect_2_object( ) ),
    // Why doesn't this work?
    construct_x_monotone_curve_2(this->
                                 traits.construct_x_monotone_curve_2_object())
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
  {
    // std::cout<< "In ArrangementPainterOstream& operator curve"<<std::endl;
    std::cout<<"In operator<< Arr_conic_traits_2 X_monotone_curve_2"<<std::endl;

    CGAL::Bbox_2 bb = curve.bbox( );
    QRectF qbb = this->convert( bb );
    // quick cull
    if ( this->clippingRect.isValid( ) &&
         ! this->clippingRect.intersects( qbb ) )
    {
      //std::cout << "quick culled curve" << std::endl;
      return *this;
    }

#if 0
    std::cout << "bottom: ("
              << this->clippingRect.bottomLeft( ).x( )
              << " "
              << this->clippingRect.bottomLeft( ).y( )
              << " "
              << this->clippingRect.bottomRight( ).x( )
              << " "
              << this->clippingRect.bottomRight( ).y( )
              << ")"
              << std::endl;
#endif

    if ( this->clippingRect.isValid( ) )
    {
      std::vector< X_monotone_curve_2 > visibleParts;
      if ( this->clippingRect.contains( qbb ) )
      {
        visibleParts.push_back( curve );
      }
      else
      {
        visibleParts = this->visibleParts( curve );
      }
      for ( unsigned int i = 0; i < visibleParts.size( ); ++i )
      {
        X_monotone_curve_2 subcurve = visibleParts[ i ];
        int n;
        if ( this->scene == NULL )
        {
          n = 100; // TODO: get an adaptive approximation
        }
        else
        {
          QGraphicsView* view = this->scene->views( ).first( );
          int xmin, xmax;
          xmin = view->mapFromScene( bb.xmin( ), bb.ymin( ) ).x( );
          xmax = view->mapFromScene( bb.xmax( ), bb.ymin( ) ).x( );
          n = xmax - xmin;
        }
        if ( n == 0 )
        {
          return *this;
        }

        std::pair<double, double>* app_pts =
          new std::pair<double, double>[n + 1];
        std::pair<double, double>* end_pts =
          subcurve.polyline_approximation(n, app_pts);
        std::pair<double, double>* p_curr = app_pts;
        std::pair<double, double>* p_next = p_curr + 1;
        int count = 0;
        do
        {
          QPointF p1( p_curr->first, p_curr->second );
          QPointF p2( p_next->first, p_next->second );
#if 0
          Segment_2 seg( p1, p2 );
          this->painterOstream << seg;
#endif
          this->qp->drawLine( p1, p2 );
          p_curr++;
          p_next++;
          ++count;
        }
        while ( p_next != end_pts );
      }
    }
    else
    { // draw the whole curve
      int n;
      if ( this->scene == NULL )
      {
        n = 100; // TODO: get an adaptive approximation
      }
      else
      {
        QGraphicsView* view = this->scene->views( ).first( );
        int xmin, xmax;
        xmin = view->mapFromScene( bb.xmin( ), bb.ymin( ) ).x( );
        xmax = view->mapFromScene( bb.xmax( ), bb.ymin( ) ).x( );
        n = xmax - xmin;

        std::cout<<"xmin\txmax"<<std::endl;
        std::cout<<xmin<<"\t"<<xmax<<std::endl;
      }
      if ( n == 0 )
      {
        return *this;
      }

      std::pair<double, double>* app_pts = new std::pair<double, double>[n + 1];
      std::pair<double, double>* end_pts =
        curve.polyline_approximation(n, app_pts);
      std::pair<double, double>* p_curr = app_pts;
      std::pair<double, double>* p_next = p_curr + 1;
      int count = 0;
      do
      {
        QPointF p1( p_curr->first, p_curr->second );
        QPointF p2( p_next->first, p_next->second );
#if 0
        Segment_2 seg( p1, p2 );
        this->painterOstream << seg;
#endif
        this->qp->drawLine( p1, p2 );
        p_curr++;
        p_next++;
        ++count;
      }
      while ( p_next != end_pts );
      //std::cout << count << " approximation points" << std::endl;
    }

    return *this;
  }

  // cloned from segtraits painter
  ArrangementPainterOstream& operator<<( const Point_2& p )
  {
    // std::cout<< "In ArrangementPainterOstream& operator Point_2"<<std::endl;
    std::cout<<"In operator<< Arr_conic_traits_2 Point_2"<<std::endl;

    QPointF qpt = this->convert( p );
    QPen savePen = this->qp->pen( );
    this->qp->setBrush( QBrush( savePen.color( ) ) );
    double radius = savePen.width( ) / 2.0;
    radius /= this->scale;

    this->qp->drawEllipse( qpt, radius, radius );

    this->qp->setBrush( QBrush( ) );
    this->qp->setPen( savePen );
    return *this;
  }

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    // std::cout<< "In ArrangementPainterOstream& operator T"<<std::endl;
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }

protected: // methods
  // Returns subcurves of curve that are actually visible in the view.
  // Assumes that clippingRect is valid.
  std::vector< X_monotone_curve_2 > visibleParts( X_monotone_curve_2 curve )
  {
    // see if we intersect the bottom edge of the viewport
    Intersect_2 intersect_2 = this->traits.intersect_2_object( );
    Point_2 bottomLeft = this->convert( this->clippingRect.bottomLeft( ) );
    Point_2 bottomRight = this->convert( this->clippingRect.bottomRight( ) );
    Point_2 topLeft = this->convert( this->clippingRect.topLeft( ) );
    Point_2 topRight = this->convert( this->clippingRect.topRight( ) );
    X_monotone_curve_2 bottom =
      this->construct_x_monotone_curve_2( bottomLeft, bottomRight );
    X_monotone_curve_2 left =
      this->construct_x_monotone_curve_2( bottomLeft, topLeft );
    X_monotone_curve_2 top =
      this->construct_x_monotone_curve_2( topLeft, topRight );
    X_monotone_curve_2 right =
      this->construct_x_monotone_curve_2( topRight, bottomRight );

    std::vector< CGAL::Object > bottomIntersections;
    std::vector< CGAL::Object > leftIntersections;
    std::vector< CGAL::Object > topIntersections;
    std::vector< CGAL::Object > rightIntersections;
    std::vector< CGAL::Object > intersections;

    intersect_2( bottom, curve, std::back_inserter( bottomIntersections ) );
    intersect_2( left, curve, std::back_inserter( leftIntersections ) );
    intersect_2( top, curve, std::back_inserter( topIntersections ) );
    intersect_2( right, curve, std::back_inserter( rightIntersections ) );
    // int total = bottomIntersections.size( )
    //   + leftIntersections.size( )
    //   + topIntersections.size( )
    //   + rightIntersections.size( );

    intersect_2( bottom, curve, std::back_inserter( intersections ) );
    intersect_2( left, curve, std::back_inserter( intersections ) );
    intersect_2( top, curve, std::back_inserter( intersections ) );
    intersect_2( right, curve, std::back_inserter( intersections ) );

    this->filterIntersectionPoints( intersections );
    //std::cout << "total intersections: " << intersections.size( )
    //          << std::endl;
    //this->printIntersectResult( intersections );

    Point_2 leftEndpt = curve.source( );
    Point_2 rightEndpt = curve.target( );

    if ( leftEndpt.x( ) > rightEndpt.x( ) )
    {
      std::swap( leftEndpt, rightEndpt );
    }

    QPointF qendpt1 = this->convert( leftEndpt );
    QPointF qendpt2 = this->convert( rightEndpt );

    std::list< Point_2 > pointList;
    for ( unsigned int i = 0; i < intersections.size( ); ++i )
    {
      CGAL::Object o = intersections[ i ];
      std::pair< Intersection_point_2, Multiplicity > pair;
      if ( CGAL::assign( pair, o ) )
      {
        Point_2 pt = pair.first;
        pointList.push_back( pt );
      }
    }

    bool includeLeftEndpoint = this->clippingRect.contains( qendpt1 );
    bool includeRightEndpoint = this->clippingRect.contains( qendpt2 );
    if ( includeLeftEndpoint )
    {
      pointList.push_front( leftEndpt );
    }

    if ( includeRightEndpoint )
    {
      pointList.push_back( rightEndpt );
    }

    Construct_x_monotone_subcurve_2< Traits > construct_x_monotone_subcurve_2;
    std::vector< X_monotone_curve_2 > clippings;
    typename std::list< Point_2 >::iterator pointListItr = pointList.begin( );
    for ( unsigned int i = 0; i < pointList.size( ); i += 2 )
    {
      Point_2 p1 = *pointListItr++;
      Point_2 p2 = *pointListItr++;
      X_monotone_curve_2 subcurve =
        construct_x_monotone_subcurve_2( curve, p1, p2 );
      clippings.push_back( subcurve );
    }

#if 0
    // std::cout << "pointList size: " << pointList.size( ) << std::endl;
    // if ( intersections.size( ) % 2 == 0 )
    // {
    //   // either both curve endpoints are in view or both are out
    //     if ( this->clippingRect.contains( qendpt1 ) )
    //     {
    //         if ( this->clippingRect.contains( qendpt2 ) )
    //         {
    //             std::cout << "both endpoints are in view" << std::endl;
    //         }
    //     }
    //     else if ( !this->clippingRect.contains( qendpt2 ) )
    //     {
    //         std::cout << "both endpoints are out of view" << std::endl;
    //     }
    // }
    // else
    // { // one curve endpoint is in view
    //     if ( this->clippingRect.contains( qendpt1 ) )
    //     {
    //         std::cout << "left endpoint is in view" << std::endl;
    //     }
    //     else if ( this->clippingRect.contains( qendpt2 ) )
    //     {
    //         std::cout << "right endpoint is in view" << std::endl;
    //     }
    // }

    std::vector< X_monotone_curve_2 > res;
    res.push_back( curve );
    return res;
#endif
    return clippings;
  }

  // keep only the intersection points ie. throw out overlapping curve segments
  void filterIntersectionPoints( std::vector< CGAL::Object >& res )
  {
    std::vector< std::pair< Intersection_point_2, Multiplicity > > tmp;

    // filter out the non-intersection point results
    for ( unsigned int i = 0; i < res.size( ); ++i )
    {
      CGAL::Object obj = res[ i ];
      std::pair< Intersection_point_2, Multiplicity > pair;
      if ( CGAL::assign( pair, obj ) )
      {
        tmp.push_back( pair );
      }
    }
    res.clear( );

    // sort the intersection points by x-coord
    Compare_intersection_point_result compare_intersection_point_result;
    std::sort( tmp.begin( ), tmp.end( ), compare_intersection_point_result );

    // box up the sorted elements
    for ( unsigned int i = 0; i < tmp.size( ); ++i )
    {
      std::pair< Intersection_point_2, Multiplicity > pair = tmp[ i ];
      CGAL::Object o = CGAL::make_object( pair );
      res.push_back( o );
    }
  }

  void printIntersectResult( const std::vector< CGAL::Object >& res )
  {
    for ( std::vector< CGAL::Object >::const_iterator it = res.begin( );
          it != res.end( ); ++it )
    {
      CGAL::Object obj = *it;
      std::pair< Intersection_point_2, Multiplicity > pair;
      if ( CGAL::assign( pair, obj ) )
      {
        Point_2 pt = pair.first;
        /* QPointF qpt = */ this->convert( pt );
        // std::cout << "(" << pt.x( ) << " " << pt.y( ) < ")" << std::endl;
      }
    }
  }

protected: // members
  Traits traits;
  //Intersect_2 intersect_2;
  Construct_x_monotone_curve_2 construct_x_monotone_curve_2;
};

template < typename Kernel_ >
class ArrangementPainterOstream< CGAL::Arr_linear_traits_2< Kernel_ > >:
  public ArrangementPainterOstreamBase< CGAL::Arr_linear_traits_2< Kernel_ > >
{
public: // typedefs
  typedef Kernel_                                       Kernel;
  typedef CGAL::Arr_linear_traits_2< Kernel >           Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()) :
    Superclass( p, clippingRectangle )
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
  {
    std::cout<<"In operator<< Arr_linear_traits_2 X_monotone_curve_2"<<std::endl;

    if ( curve.is_segment( ) )
    {
      Segment_2 seg = curve.segment( );

      // skip segments outside our view
      QRectF seg_bb = this->convert( seg.bbox( ) );
      if ( this->clippingRect.isValid( ) &&
           ! this->clippingRect.intersects( seg_bb ) )
      {
        return *this;
      }

      this->painterOstream << seg;
    }
    else if ( curve.is_ray( ) )
    {
      Ray_2 ray = curve.ray( );
      QLineF qseg = this->convert( ray );
      if ( qseg.isNull( ) )
      { // it's out of view
        return *this;
      }
      Segment_2 seg = this->convert( qseg );
      this-> painterOstream << seg;
    }
    else // curve.is_line( )
    {
      Line_2 line = curve.line( );
      QLineF qseg = this->convert( line );
      if ( qseg.isNull( ) )
      { // it's out of view
        return *this;
      }
      Segment_2 seg = this->convert( qseg );
      this-> painterOstream << seg;
    }
    return *this;
  }

  ArrangementPainterOstream& operator<<( const Point_2& p )
  {

    std::cout<<"In operator<< Arr_linear_traits_2 Point_2"<<std::endl;

    QPointF qpt = this->convert( p );
    // clip the point if possible
    if ( this->clippingRect.isValid( ) &&
         ! this->clippingRect.contains( qpt ) )
    {
      return *this;
    }

    QPen savePen = this->qp->pen( );
    this->qp->setBrush( QBrush( savePen.color( ) ) );
    double radius = savePen.width( ) / 2.0;
    radius /= this->scale;

    this->qp->drawEllipse( qpt, radius, radius );

    this->qp->setBrush( QBrush( ) );
    this->qp->setPen( savePen );
    return *this;
  }

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }
};

template < typename CircularKernel >
class ArrangementPainterOstream< CGAL::Arr_circular_arc_traits_2<
                                   CircularKernel > >:
  public ArrangementPainterOstreamBase< CGAL::Arr_circular_arc_traits_2<
                                          CircularKernel > >
{
public:
  typedef CircularKernel                                Kernel;
  typedef CGAL::Arr_circular_arc_traits_2< Kernel >     Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Superclass::Segment_2                Segment_2;
  typedef typename Superclass::Ray_2                    Ray_2;
  typedef typename Superclass::Line_2                   Line_2;
  typedef typename Superclass::Triangle_2               Triangle_2;
  typedef typename Superclass::Iso_rectangle_2          Iso_rectangle_2;
  typedef typename Superclass::Circle_2                 Circle_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()):
    Superclass( p, clippingRectangle )
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
  {
    std::cout<<"In operator<< Arr_circular_arc_traits_2 X_monotone_curve_2"<<std::endl;

    this->painterOstream << curve;
    return *this;
  }

  ArrangementPainterOstream& operator<<( const Point_2& p )
  {
    std::cout<<"In operator<< Arr_circular_arc_traits_2 Point_2"<<std::endl;
    QPointF qpt = this->convert( p );
    // clip the point if possible
    if ( this->clippingRect.isValid( ) &&
         ! this->clippingRect.contains( qpt ) )
    {
      return *this;
    }

    QPen savePen = this->qp->pen( );
    this->qp->setBrush( QBrush( savePen.color( ) ) );
    double radius = savePen.width( ) / 2.0;
    radius /= this->scale;

    this->qp->drawEllipse( qpt, radius, radius );

    this->qp->setBrush( QBrush( ) );
    this->qp->setPen( savePen );
    return *this;
  }

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }
};


template < typename Coefficient_ >
class ArrangementPainterOstream< CGAL::Arr_algebraic_segment_traits_2<
                                   Coefficient_ > >:
  public ArrangementPainterOstreamBase< CGAL::Arr_algebraic_segment_traits_2<
                                          Coefficient_ > >
{
public:
  typedef Coefficient_                                  Coefficient;
  typedef typename CGAL::Arr_algebraic_segment_traits_2< Coefficient >
                                                        Traits;
  typedef ArrangementPainterOstreamBase< Traits >       Superclass;
  typedef typename Superclass::Point_2                  Point_2;
  typedef typename Traits::X_monotone_curve_2           X_monotone_curve_2;
  typedef typename Traits::Curve_2                      Curve_2;
  typedef typename Traits::Construct_x_monotone_segment_2
    Construct_x_monotone_segment_2;
  typedef typename Traits::Point_2                      Intersection_point_2;
  typedef typename Traits::Intersect_2                  Intersect_2;
  typedef typename Traits::Multiplicity                 Multiplicity;

public:
  /*! Constructor */
  ArrangementPainterOstream(QPainter* p, QRectF clippingRectangle = QRectF()):
    Superclass( p, clippingRectangle ),
    clipRect(clippingRectangle),
    construct_x_monotone_segment_2(this->
                                 traits.construct_x_monotone_segment_2_object())
  { }

  /*! Destructor (virtual) */
  virtual ~ArrangementPainterOstream() {}

  class Compare_intersection_point_result
  {
  public:
    typedef std::pair< Intersection_point_2, Multiplicity > Result;
    // returns whether the point1 < point2, using x-coord to compare
    bool operator()( const Result& o1, const Result& o2 )
    {
      Intersection_point_2 p1 = o1.first;
      Intersection_point_2 p2 = o2.first;
      return ( p1.x( ) < p2.x( ) );
    }
  };

protected:
  std::vector< Intersection_point_2 > visibleParts( X_monotone_curve_2 curve )
  {
    std::cout<<"In visibleParts\n";
    // see if we intersect the bottom edge of the viewport
    typename Traits::Construct_x_monotone_segment_2 constructSegment =
      traits.construct_x_monotone_segment_2_object( );
    typename Traits::Construct_point_2 constructPoint =
      traits.construct_point_2_object( );

    Intersect_2 intersect_2 = this->traits.intersect_2_object( );

    Intersection_point_2 bottomLeft = constructPoint(this->clipRect.bottom(), this->clipRect.left());
    Intersection_point_2 bottomRight = constructPoint(this->clipRect.bottom(), this->clipRect.right());
    Intersection_point_2 topLeft = constructPoint(this->clipRect.top(), this->clipRect.left());
    Intersection_point_2 topRight = constructPoint(this->clipRect.top(), this->clipRect.right());

    std::cout<<this->clipRect.bottom()<<'\t'<<this->clipRect.left()<<std::endl;
    std::cout<<this->clipRect.bottom()<<'\t'<<this->clipRect.right()<<std::endl;
    std::cout<<this->clipRect.top()<<'\t'<<this->clipRect.left()<<std::endl;
    std::cout<<this->clipRect.top()<<'\t'<<this->clipRect.right()<<std::endl;

    std::vector< X_monotone_curve_2 > vec;

    constructSegment( bottomLeft, bottomRight, std::back_inserter(vec) );
    X_monotone_curve_2 bottom = vec[0];
    vec.clear();

    constructSegment( bottomLeft, topLeft, std::back_inserter(vec) );
    X_monotone_curve_2 left = vec[0];
    vec.clear();

    constructSegment( topLeft, topRight, std::back_inserter(vec) );
    X_monotone_curve_2 top = vec[0];
    vec.clear();

    constructSegment( topRight, bottomRight, std::back_inserter(vec) );
    X_monotone_curve_2 right = vec[0];
    vec.clear();

    std::vector< CGAL::Object > intersections;

    intersect_2( bottom, curve, std::back_inserter( intersections ) );
    intersect_2( left, curve, std::back_inserter( intersections ) );
    intersect_2( top, curve, std::back_inserter( intersections ) );
    intersect_2( right, curve, std::back_inserter( intersections ) );

    this->filterIntersectionPoints( intersections );

    std::cout << "total intersections: " << intersections.size( )
             << std::endl;
    //this->printIntersectResult( intersections );



    std::list< Intersection_point_2 > pointList;
    for ( unsigned int i = 0; i < intersections.size( ); ++i )
    {
      CGAL::Object o = intersections[ i ];
      std::pair< Intersection_point_2, Multiplicity > pair;
      if ( CGAL::assign( pair, o ) )
      {
        Intersection_point_2 pt = pair.first;
        pointList.push_back( pt );
      }
    }


    // std::cout<<"leftEndpt\n";
    // std::cout<<leftEndpt.x()<<'\t'<<leftEndpt.y()<<std::endl;

    // std::cout<<"rightEndpt\n";
    // std::cout<<rightEndpt.x()<<'\t'<<rightEndpt.y()<<std::endl;
    
    Arr_compute_y_at_x_2<Traits> arr_compute_y_at_x_2;
    CGAL::Bbox_2 bb = curve.bbox( );

    std::cout<<bb.xmin()<<"\t";
    std::cout<<bb.xmax()<<"\t";
    std::cout<<bb.ymin()<<"\t";
    std::cout<<bb.ymax()<<"\t"<<std::endl;

    if ( !std::isinf(bb.xmin()) )
    {
      double left_endpoint_x = bb.xmin();
      double left_endpoint_y = arr_compute_y_at_x_2.approx(curve, left_endpoint_x);

      std::cout<<"In !std::isinf(bb.xmin():\t"<<left_endpoint_y<<std::endl;
      Intersection_point_2 left_endpoint = constructPoint(left_endpoint_x, left_endpoint_y);
      pointList.push_front( left_endpoint );
    }

    if ( !std::isinf(bb.xmax()) )
    {
      double right_endpoint_x = bb.xmax();
      double right_endpoint_y = arr_compute_y_at_x_2.approx(curve, right_endpoint_x);

      std::cout<<"In !std::isinf(bb.xmax():\t"<<right_endpoint_y<<std::endl;

      Intersection_point_2 right_endpoint = constructPoint(right_endpoint_x, right_endpoint_y);
      pointList.push_back( right_endpoint );
    }
    // if ( leftEndpt.x( ) > rightEndpt.x( ) )
    // {
    //   std::swap( leftEndpt, rightEndpt );
    // }

    // QPointF qendpt1 = this->convert( leftEndpt );
    // QPointF qendpt2 = this->convert( rightEndpt );

    // bool includeLeftEndpoint = this->clippingRect.contains( qendpt1 );
    // bool includeRightEndpoint = this->clippingRect.contains( qendpt2 );
    // if ( includeLeftEndpoint )
    // {
    //   pointList.push_front( leftEndpt );
    // }

    // if ( includeRightEndpoint )
    // {
    //   pointList.push_back( rightEndpt );
    // }
    std::cout << "total points in pointList: " << pointList.size( )
             << std::endl;

    // std::vector< X_monotone_curve_2 > clippings;
    std::vector<Intersection_point_2> seperation_points;

    typename std::list< Intersection_point_2 >::iterator pointListItr = pointList.begin( );
    for ( ; pointListItr != pointList.end(); pointListItr++ )
    {
      std::cout<<"before construct_x_monotone_segment_2:\n";
      seperation_points.push_back(*pointListItr);
      // Intersection_point_2 p1 = *pointListItr++;
      // Intersection_point_2 p2 = *pointListItr++;
      // X_monotone_curve_2 subcurve =
      //   construct_x_monotone_subcurve_2( curve, p1, p2 );
      // clippings.push_back( subcurve );
      
      // std::cout<<CGAL::to_double(p1.x())<<"\t"<<CGAL::to_double(p1.y())<<"\n";
      // std::cout<<CGAL::to_double(p2.x())<<"\t"<<CGAL::to_double(p2.y())<<"\n";

      // this->construct_x_monotone_segment_2(curve.curve(), p1, p2, std::back_inserter(clippings) );
      std::cout<<"after construct_x_monotone_segment_2:\n"<<std::endl;
    }

    std::cout<<"Leaving visibleParts\n";
    return seperation_points;
  }

  void filterIntersectionPoints( std::vector< CGAL::Object >& res )
  {
    std::vector< std::pair< Intersection_point_2, Multiplicity > > tmp;

    // filter out the non-intersection point results
    for ( unsigned int i = 0; i < res.size( ); ++i )
    {
      CGAL::Object obj = res[ i ];
      std::pair< Intersection_point_2, Multiplicity > pair;
      if ( CGAL::assign( pair, obj ) )
      {
        tmp.push_back( pair );
      }
    }
    res.clear( );

    // sort the intersection points by x-coord
    Compare_intersection_point_result compare_intersection_point_result;
    std::sort( tmp.begin( ), tmp.end( ), compare_intersection_point_result );

    // box up the sorted elements
    for ( unsigned int i = 0; i < tmp.size( ); ++i )
    {
      std::pair< Intersection_point_2, Multiplicity > pair = tmp[ i ];
      CGAL::Object o = CGAL::make_object( pair );
      res.push_back( o );
    }
  }

public: // methods
  ArrangementPainterOstream& operator<<( const X_monotone_curve_2& curve )
  {
    std::cout<<"In operator<< Arr_algebraic_segment_traits_2 X_monotone_curve_2"<<std::endl;

    Bbox_2 bbox = curve.bbox();

    std::cout<<bbox.xmin()<<"\t";
    std::cout<<bbox.xmax()<<"\t";
    std::cout<<bbox.ymin()<<"\t";
    std::cout<<bbox.ymax()<<"\t"<<std::endl;

    if (!std::isinf(bbox.xmin()) 
      && !std::isinf(bbox.xmax())
      && !std::isinf(bbox.ymin())
      && !std::isinf(bbox.ymax()))
    {
      int n = 100;

      std::pair<double, double>* app_pts =
        new std::pair<double, double>[n + 1];
      std::pair<double, double>* end_pts =
        this->polyline_approximation(curve, bbox.xmin(), bbox.xmax(), n, app_pts);
      std::pair<double, double>* p_curr = app_pts;
      std::pair<double, double>* p_next = p_curr + 1;
      int count = 0;
      do
      {
        QPointF p1( p_curr->first, p_curr->second );
        QPointF p2( p_next->first, p_next->second );

        this->qp->drawLine( p1, p2 );
        p_curr++;
        p_next++;
        ++count;
      } while ( p_next != end_pts );

      std::cout<<"count\t"<<count<<std::endl;
      delete [] app_pts;

      std::cout<<"Leaving Arr_algebraic_segment_traits_2 X_monotone_curve_2"<<std::endl;
      return *this;
    }

    std::cout<<"Bbox_2 has infinite bounding box\n";
    std::vector< Intersection_point_2 > seperation_points = this->visibleParts(curve);
    int seperation_points_size = seperation_points.size( );

    for ( unsigned int i = 0; i < seperation_points_size; ++i )
    {
      if (i >= seperation_points_size || (i+1) >= seperation_points_size)
      {
        break;
      }

      Intersection_point_2 p1 = seperation_points[i];
      Intersection_point_2 p2 = seperation_points[i+1];
      int n = 100;

      std::pair<double, double>* app_pts =
        new std::pair<double, double>[n + 1];
      std::pair<double, double>* end_pts =
        this->polyline_approximation(curve, CGAL::to_double(p1.x()), CGAL::to_double(p2.x()), n, app_pts);
      std::pair<double, double>* p_curr = app_pts;
      std::pair<double, double>* p_next = p_curr + 1;
      int count = 0;
      do
      {
        QPointF p1( p_curr->first, p_curr->second );
        QPointF p2( p_next->first, p_next->second );

        this->qp->drawLine( p1, p2 );
        p_curr++;
        p_next++;
        ++count;
      } while ( p_next != end_pts );

      std::cout<<"count\t"<<count<<std::endl;
      delete [] app_pts;
    }
#if 0
    CGAL::Bbox_2 bb = curve.bbox( );

    // QGraphicsView* view = this->scene->views( ).first( );
    double xmin, xmax, ymin, ymax;
    xmin = bb.xmin( );
    xmax = bb.xmax( );
    ymin = bb.ymin( );
    ymax = bb.ymax( );
    // xmin = view->mapFromScene( bb.xmin( ), bb.ymin( ) ).x( );
    // xmax = view->mapFromScene( bb.xmax( ), bb.ymin( ) ).x( );

    std::cout<<"xmin\txmax\tymin\tymax"<<std::endl;
    std::cout<<xmin<<"\t"<<xmax<<"\t"<<ymin<<"\t"<<ymax<<std::endl;
#endif


    std::cout<<"Leaving Arr_algebraic_segment_traits_2 X_monotone_curve_2"<<std::endl;
    return *this;
  }

  std::pair<double, double>* polyline_approximation(const X_monotone_curve_2& curve,
    double source_x, double target_x, const int& num_of_points, 
    std::pair<double, double>* target_memory)
  {
    CGAL::Bbox_2 bb = curve.bbox( );
    Arr_compute_y_at_x_2<Traits> arr_compute_y_at_x_2;

    int cnt = 0;

    double xmin = source_x;
    double xmax = target_x;

    std::cout<<"In polyline_approximation\t";
    std::cout<<xmin<<"\t";
    std::cout<<xmax<<std::endl;

    double interval = (xmax - xmin)/(num_of_points-1);
    
    double x_cur = xmin;
    while(cnt < num_of_points)
    {

      double y_cur = arr_compute_y_at_x_2.approx(curve, x_cur);

      target_memory[cnt++] = std::pair<double, double>(x_cur, y_cur);
      x_cur += interval;
    }

    return target_memory + num_of_points;
  }

  ArrangementPainterOstream& operator<<( const Point_2& p )
  {
    std::cout<<"In operator<< Arr_algebraic_segment_traits_2 Point_2"<<std::endl;

    QPointF qpt = this->convert( p );

    QPen savePen = this->qp->pen( );
    this->qp->setBrush( QBrush( savePen.color( ) ) );
    double radius = savePen.width( ) / 2.0;
    radius /= this->scale;

    this->qp->drawEllipse( qpt, radius, radius );

    this->qp->setBrush( QBrush( ) );
    this->qp->setPen( savePen );

    // }
#if 0

    QPainter *ppnt = &ws.get_painter();
    QPen old_pen = ppnt->pen();
    ppnt->setPen(QPen(Qt::NoPen));

    unsigned sz = CGAL_REND_PT_RADIUS;
    ppnt->drawEllipse(coord.first - sz, ws.height() - coord.second - sz,
                      sz*2, sz*2);
    ppnt->setPen(old_pen);
#endif

    return *this;
  }

  template < typename T >
  ArrangementPainterOstream& operator<<( const T& p )
  {
    (*(static_cast< Superclass* >(this)) << p);
    return *this;
  }

protected:
  QRectF clipRect;
  Traits traits;
  Construct_x_monotone_segment_2 construct_x_monotone_segment_2;

};


} // namespace Qt
} // namespace CGAL

#endif // CGAL_QT_ARRANGEMENT_PAINTER_OSTREAM_H
