#include "support_plane_objects.h"
#include "support_plane.h"
#include "universe.h"

namespace JPTD {

Polygon_Edge::Polygon_Edge(const int _id_plane, Polygon_Vertex* _v1, Polygon_Vertex* _v2)
	: Support_Plane_Object(_id_plane),
	v1 (_v1),
	v2 (_v2)
{
	v1->add(this);
	v2->add(this);

	is_constrained = false;
	constraint = Constraint();
}



Polygon_Edge::Polygon_Edge(const int _id_plane, Polygon_Vertex* _v1, Polygon_Vertex* _v2, const Constraint & C)
	: Polygon_Edge(_id_plane, _v1, _v2)
{
	is_constrained = true;
	constraint = C;
}



Polygon_Edge::Polygon_Edge(const int _id_plane, Polygon_Vertex* _v1, Polygon_Vertex* _v2, Intersection_Line* I, Sign epsilon)
	: Polygon_Edge(_id_plane, _v1, _v2)
{
	is_constrained = true;
	constraint = Constraint(I, epsilon);
}



Polygon_Edge::~Polygon_Edge()
{
	v1->remove(this);
	v2->remove(this);
}



bool Polygon_Edge::is_constrained_by(Intersection_Line* I) const
{
	return (is_constrained && constraint.first == I);
}



Polygon_Vertex* Polygon_Edge::intersection(Intersection_Line* I, Sign s, const FT & t, const int K, Event_Flags flags) const
{
	CGAL_Point_2 M;
	CGAL_Vector_2 dM;
	intersection_pt_dir(I, t, M, dM);

	return intersection(I, s, t, M, dM, K, flags);
}



Polygon_Vertex* Polygon_Edge::intersection(Intersection_Line* I, Sign s, const FT & t, const CGAL_Point_2 & M, const CGAL_Vector_2 & dM, const int K, Event_Flags flags) const
{
	Polygon_Vertex* v = nullptr;
	Constraint C = Constraint(I, s);
	if (!is_constrained) {
		v = new Polygon_Vertex(I->id_plane, t, M, dM, C, nullptr, K, flags);
	} else {
		v = new Polygon_Vertex(I->id_plane, t, M, C, constraint);
	}

	return v;
}



void Polygon_Edge::intersection_pt_dir(Intersection_Line* I, const FT & t, CGAL_Point_2 & M, CGAL_Vector_2 & dM) const
{
	// We assume that edge propagates within the frame of a homeomorphic transformation,
	// that's why the speed of a point that represents the intersection of an edge e and a line I is constant.
	// We suppose that, at the time when this function is called, the edge is colliding I.
	
	CGAL_Point_2 v1_t = v1->pt(t), v2_t = v2->pt(t);

	const CGAL_Line_2 & L = I->line;
	CGAL_Segment_2 S_t (v1_t, v2_t);

	// Step 1.
	// Computes the intersection of L and S_t
	
	// For better efficiency, we do not compute a intersection of a line L and a segment S_t,
	// but the intersection of two lines : L and S_t.supporting_line()

	CGAL_Point_2 M_t;
	bool M_t_exists = false;

	CGAL_Line_2 S_t_line = S_t.supporting_line();

	CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Line_2)>::type object_1 = CGAL::intersection(L, S_t_line);
	if (object_1) {
		if (const CGAL_Point_2* ptr = boost::get<CGAL_Point_2>(&*object_1)) {
			M_t = (*ptr);
			M_t_exists = true;
		}
	}

	if (!M_t_exists) {
		throw std::logic_error("Error : unexpected behavior of CGAL::intersection.");
	}

	// Step 2.
	// Computes the support line of the edge at time t + 1, and its intersection with L.
	// This lets us compute the speed dM of the seeked intersection point, and finally the initial point M.

	CGAL_Point_2 v1_u = v1_t + v1->dM;
	CGAL_Point_2 v2_u = v2_t + v2->dM;
	CGAL_Line_2 S_u (v1_u, v2_u);

	CGAL_Point_2 M_u;
	bool M_u_exists = false;

	CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Line_2)>::type object_2 = CGAL::intersection(L, S_u);
	if (object_2) {
		if (const CGAL_Point_2* ptr = boost::get<CGAL_Point_2>(&*object_2)) {
			M_u = (*ptr);
			M_u_exists = true;
		}
	}

	if (!M_u_exists) {
		throw std::logic_error("Error : unexpected behavior of CGAL::intersection.");
	}

	dM = M_u - M_t;
	M = M_t - t * dM;
}



void Polygon_Edge::intersection_pt_dir(Intersection_Line* I, Polygon_Vertex* v1, Polygon_Vertex* v2, const FT & t, const CGAL_Point_2 & v1_t, CGAL_Point_2 & M, CGAL_Vector_2 & dM)
{
	const CGAL_Point_2 v2_t = v2->pt(t);
	const CGAL_Line_2 & L = I->line;

	// The first part of the previous function is now useless : M_t = v1_t.

	// We compute the support line of the edge (v1 v2) at time t + 1, and its intersection with L.
	// This lets us compute the speed dM of the seeked intersection point, and finally the initial point M.

	CGAL_Point_2 v1_u = v1_t + v1->dM;
	CGAL_Point_2 v2_u = v2_t + v2->dM;
	CGAL_Line_2 S_u (v1_u, v2_u);

	CGAL_Point_2 M_u;
	bool M_u_exists = false;

	CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Line_2)>::type object_2 = CGAL::intersection(L, S_u);
	if (object_2) {
		if (const CGAL_Point_2* ptr = boost::get<CGAL_Point_2>(&*object_2)) {
			M_u = (*ptr);
			M_u_exists = true;
		}
	}

	if (!M_u_exists) {
		throw std::logic_error("Error : unexpected behavior of CGAL::intersection.");
	}

	dM = M_u - v1_t;
	M = v1_t - t * dM;
}



bool Polygon_Edge::is_adjacent_to_segment() const
{
	// We assume that vertices v1 and v2 are stopped
	assert(!v1->is_active && !v2->is_active);
	Intersection_Line* I = constraint.first;
	assert(I != nullptr);

	// Gets a definition of the edge at t = + \infty, and the coordinates of its midpoint
	const CGAL_Point_2 & A = v1->M, & B = v2->M;
	const CGAL_Point_2 M = CGAL::midpoint(A, B);

	for (std::list<Segment *>::iterator it_s = I->segments.begin() ; it_s != I->segments.end() ; it_s++) {
		if (Planar_Segment* pl_s = dynamic_cast<Planar_Segment*>(*it_s)) {
			// The segment is adjacent to an edge of the Support_Plane
			return true;
		} 
		
		else if (Polygon_Segment* po_s = dynamic_cast<Polygon_Segment*>(*it_s)) {
			// We check if the point M computed before belongs to the segment
			const CGAL_Point_2 S = po_s->origin(), T = po_s->end();
			const CGAL_Segment_2 ST (S, T);
			//if (ST.has_on(M)) return true;

			FT x_s = S.x(), x_t = T.x(), x = M.x();
			if (x_s != x_t) {
				if ((x_s <= x && x <= x_t) || (x_t <= x && x <= x_s)) return true;
			} else {
				FT y_s = S.y(), y_t = T.y(), y = M.y();
				if ((y_s <= y && y <= y_t) || (y_t <= y && y <= y_s)) return true;
			}
		}
	}

	return false;
}

}