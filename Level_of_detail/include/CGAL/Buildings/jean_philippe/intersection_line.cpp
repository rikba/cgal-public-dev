#include "support_plane_objects.h"
#include "support_plane.h"
#include "universe.h"
#include <list>

namespace JPTD {
Intersection_Line::Intersection_Line(const int _id_plane, const CGAL_Line_2 & _line, int intersected)
	: Support_Plane_Object(_id_plane),
	line (_line),
	_a (line.a()),
	_b (line.b()),
	_c (line.c()),
	hint_a (to_double(_a)),
	hint_b (to_double(_b)),
	hint_c (to_double(_c))
{
	// If the id of the plane represented by the line is lower than 5,
	// then it corresponds to a facet of the bounding box.
	is_border = (intersected <= 5);

	// Considers by default that the line is inside the bounding polygon of the associated plane
	is_inside = true;

	planes = std::list<int>(1, intersected);
	segments = std::list<Segment *>();
}



Intersection_Line::~Intersection_Line()
{

}



void Intersection_Line::mark_as_intersected(int intersected)
{
	planes.push_back(intersected);
}



bool Intersection_Line::intersects(const int id_plane) const
{
	for (std::list<int>::const_iterator it_p = planes.begin() ; it_p != planes.end() ; it_p++) {
		if ((*it_p) == id_plane) {
			return true;
		}
	}
	return false;
}



void Intersection_Line::set_inside(const bool _is_inside)
{
	is_inside = _is_inside;
}



bool Intersection_Line::includes(const CGAL_Point_2 & M) const
{
	return (_a * M.x() + _b * M.y() + _c == 0);
}



Sign Intersection_Line::sign(const CGAL_Point_2 & pt) const 
{
	FT t = _a * pt.x() + _b * pt.y() + _c;
	if (t > 0) {
		return PLUS;
	} else if (t < 0) {
		return MINUS;
	} else {
		return ZERO;
	}
}



Sign Intersection_Line::sign(Polygon_Vertex* v, const FT t) const
{
	int K = 0;

	CGAL_Point_2 M = v->pt(t);
	FT dist = _a * M.x() + _b * M.y() + _c;

	while (dist == 0) {
		// While it is not clear whether the vertex v, at time t, 
		// is clearly on one or another side of I, we exponentially increment
		// the argument of function v->pt()

		K += 1;
		M = v->pt(t - int(pow(10, K)));
		dist = _a * M.x() + _b * M.y() + _c;
	}

	if (dist > 0) {
		return PLUS;
	} else {
		return MINUS;
	}
}



bool Intersection_Line::is_parallel(Intersection_Line* I) const
{
	assert(this != I);

	const CGAL_Line_2 & I_line = I->line;
	return (CGAL::determinant(line.to_vector(), I_line.to_vector()) == 0);
}



bool Intersection_Line::exist_segments_including_point_outside_intersections(const CGAL_Point_2 & V_t, const FT & t) const
{
	// We suppose that this line is not a border of the bounding box.
	// Here, we determine if there exists a couple of segments, of positive and negative signs, that include V_t at time t.

	bool po_s_plus = false, po_s_minus = false;

	for (std::list<Segment*>::const_iterator it_s = segments.begin() ; it_s != segments.end() ; it_s++) {
		if (Polygon_Segment* po_s = dynamic_cast<Polygon_Segment*>(*it_s)) {
			
			// If the segment is of positive sign and if there already exists a segment of positive sign that includes V_t, continues
			// Same operation if the segment has negative sign
			if ((po_s->sign == PLUS && po_s_plus) || (po_s->sign == MINUS && po_s_minus)) {
				continue;
			}

			// Performs the operation
			if (po_s->includes_point_on_support_line(V_t, t)) {
				if (po_s->sign == PLUS) {
					po_s_plus = true;
				} else {
					po_s_minus = true;
				}
			}

			// Early exit if possible
			if (po_s_plus && po_s_minus) return true;
		}
	}

	return false;
}



bool Intersection_Line::exist_segments_including_point_at_intersection(const CGAL_Point_2 & V_t, const Constraint & C, const FT & t) const
{
	// This is the same code as above,
	// except that we do not call the same method of the Polygon_Segment :
	// this time we make use of the fact that V_t is at the intersection this line and C.first

	bool po_s_plus = false, po_s_minus = false;

	for (std::list<Segment*>::const_iterator it_s = segments.begin() ; it_s != segments.end() ; it_s++) {
		if (Polygon_Segment* po_s = dynamic_cast<Polygon_Segment*>(*it_s)) {
			
			if ((po_s->sign == PLUS && po_s_plus) || (po_s->sign == MINUS && po_s_minus)) {
				continue;
			}

			if (po_s->includes_point_at_intersection(V_t, C, t)) {
				if (po_s->sign == PLUS) {
					po_s_plus = true;
				} else {
					po_s_minus = true;
				}
			}

			if (po_s_plus && po_s_minus) return true;
		}
	}

	return false;
}



bool Intersection_Line::exist_segments_including_point_at_intersection(const CGAL_Point_2 & V_t, const std::list<Constraint> & C_limits, const FT & t) const
{
	// This is the same code as above
	// It takes into account the case of a multi-line intersection

	bool po_s_plus = false, po_s_minus = false;
		for (std::list<Segment*>::const_iterator it_s = segments.begin() ; it_s != segments.end() ; it_s++) {
		if (Polygon_Segment* po_s = dynamic_cast<Polygon_Segment*>(*it_s)) {
			
			if ((po_s->sign == PLUS && po_s_plus) || (po_s->sign == MINUS && po_s_minus)) {
				continue;
			}

			if (po_s->includes_point_at_intersection(V_t, C_limits, t)) {
				if (po_s->sign == PLUS) {
					po_s_plus = true;
				} else {
					po_s_minus = true;
				}
			}

			if (po_s_plus && po_s_minus) return true;
		}
	}

	return false;
}



bool Intersection_Line::exists_segment_adjacent_to_edge(Polygon_Edge* e) const
{
	if (is_border) return true;

	// e = (v1 v2), where v1 and v2 are double-constrained vertices
	// We get the constraints C_1 and C_2 that correspond to such vertices,
	// in which this object is not involved

	Constraint C_1 = e->v1->get_constraint();
	if (C_1.first == this) C_1 = e->v1->get_second_constraint();

	Constraint C_2 = e->v2->get_constraint();
	if (C_2.first == this) C_2 = e->v2->get_second_constraint();

	// Gets vertices

	const CGAL_Point_2 & V_1 = e->v1->M;
	const CGAL_Point_2 & V_2 = e->v2->M;
	
	// Gets lines that are concurrent with this object, C_1.first and C_2.first
	std::list<Intersection_Line*> CL_1, CL_2, CL;

	Support_Plane* SP = Universe::map_of_planes[this->id_plane];
	SP->get_concurrent_lines(this, C_1.first, CL_1);
	SP->get_concurrent_lines(this, C_2.first, CL_2);

	std::copy(CL_1.begin(), CL_1.end(), std::back_inserter(CL));
	std::copy(CL_2.begin(), CL_2.end(), std::back_inserter(CL));

	// Loops on all segments until find one that contains the edge

	for (std::list<Segment*>::const_iterator it_s = segments.begin() ; it_s != segments.end() ; it_s++) {
		if (Polygon_Segment* po_s = dynamic_cast<Polygon_Segment*>(*it_s)) {
			if (po_s->includes_edge(V_1, V_2, C_1, C_2, CL)) {
				// Early exit
				return true;
			}
		}
	}

	return false;
}


const FT & Intersection_Line::a() const 
{ 
	return _a; 
}


const FT & Intersection_Line::b() const 
{ 
	return _b;
}


const FT & Intersection_Line::c() const 
{ 
	return _c; 
}
}