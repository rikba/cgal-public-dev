#include "support_plane_objects.h"
#include "support_plane.h"
#include "universe.h"
#include "event_queue.h"

namespace JPTD {

Polygon_Segment::Polygon_Segment(const int _id_plane, const FT & _t_init, const Constraint & C_support)
	: Segment(_id_plane, C_support.first)
{
	/*if (id_object == 69874) {
		std::cout << "break" << std::endl;
	}*/

	// A private constructor, that factorizes the code common to most constructors.

	Universe::map_of_planes[_id_plane]->segments[id_object] = this;

	sign = C_support.second;

	t_init = _t_init;
	t_stop = FLT_MAX;

	Tr = nullptr;
	Tr_previous = std::list<Segment_Translation *>();

	C_init.first = nullptr;
	C_stop.first = nullptr;
	C_crossed = std::list<Intersection_Line *>();

	opposite = nullptr;
}



Polygon_Segment::Polygon_Segment(const int _id_plane, const FT & _t_init, const Constraint & C_support, Polygon_Vertex* & v, const CGAL_Point_2 & _O, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA)
	: Polygon_Segment (_id_plane, _t_init, C_support)
{
	// An object of class Polygon_Segment can be divided into three different parts :
	//
	// a) An initial segment [O A_0] where A_0 = Tr_previous.front()->A. [O A_0] is null if the segment is created at t > 0.
	// b) A list of segments [A_i A_{i + 1}] where each segment corresponds to the path followed by the ray 
	//    between t->int_start and t->int_end starting from t->A, where t is an element of Tr_previous.
	// c) Finally, the half-line [Tr->A (Tr->A + t * Tr->dA)). It doesn't exist if the segment is not active anymore.

	Tr = new Segment_Translation(_t_init, _A, _dA);
	Tr_previous.push_back(new Segment_Translation(_t_init, _O, _A));

	// This construtor is typically called during the initialization process,
	// or when a Polygon_Edge with two unconstrained vertices intersects an Intersection_Line.

	v->set_as_guided_segment(this);
}



Polygon_Segment::Polygon_Segment(const int _id_plane, const FT & _t_init, const Constraint & C_support, Polygon_Vertex* & v, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA)
	: Polygon_Segment (_id_plane, _t_init, C_support)
{
	// This constructor is typically called as an unconstrained Polygon_Vertex intersects an Intersection_Line.

	Tr = new Segment_Translation(_t_init, _A, _dA);

	v->set_as_guided_segment(this);
}



Polygon_Segment::Polygon_Segment(const int _id_plane, const FT & _t_init, const Constraint & _C_init, const Constraint & C_support, Polygon_Vertex* & v, 
	const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA)
	: Polygon_Segment (_id_plane, _t_init, C_support)
{
	// This constructor is typically called as an constrained Polygon_Vertex intersects an Intersection_Line.

	Tr = new Segment_Translation(_t_init, _A, _dA);

	C_init = _C_init;

	v->set_as_guided_segment(this);
}



Polygon_Segment::Polygon_Segment(const int _id_plane, const FT & _t_init, const Constraint & _C_init, const Constraint & C_support, Polygon_Vertex* & v, 
	const CGAL_Point_2 & _O, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA)
	: Polygon_Segment (_id_plane, _t_init, C_support)
{
	// This constructor builds Polygon_Segments which are created at time t.
	// Contrary to before, at the time of its creation, the Polygon_Segment is a segment [OA] and not a single point {A}.
	// It is typically called when a Polygon_Edge with one constrained vertex intersects an Intersection_Line.

	Tr = new Segment_Translation(_t_init, _A, _dA);
	Tr_previous.push_back(new Segment_Translation(_t_init, _O, _A));

	C_init = _C_init;

	v->set_as_guided_segment(this);
}



Polygon_Segment::Polygon_Segment(const int _id_plane, const FT & t, const Constraint & _C_init, const Constraint & C_support, const Constraint & _C_stop, 
	const CGAL_Point_2 & A, const CGAL_Point_2 & B)
	: Polygon_Segment (_id_plane, t, C_support)
{
	// This constructor builds Polygon_Segments that are constant.
	// At the time of their creation, their initial and final points and constraints are already determined.
	// It is typically called when a Polygon_Edge with two constrained vertices intersects an Intersection_Line.

	t_stop = t;

	Tr_previous.push_back(new Segment_Translation(t, A, B));

	C_init = _C_init;
	C_stop = _C_stop;
}



Polygon_Segment::~Polygon_Segment()
{
	if (Tr != nullptr) {
		delete Tr;
		Tr = nullptr;
	}

	for (std::list<Segment_Translation*>::iterator it = Tr_previous.begin() ; it != Tr_previous.end() ; it++) delete (*it);
	Tr_previous.clear();
}



bool Polygon_Segment::stopped() const 
{ 
	return (t_stop < FLT_MAX); 
}



FT Polygon_Segment::t_min() const 
{ 
	return (Tr_previous.empty() ? Tr->t_int_start : Tr_previous.front()->t_int_start); 
}



CGAL_Point_2 Polygon_Segment::origin() const
{
	if (Tr_previous.empty()) {
		return Tr->A;
	} else {
		return Tr_previous.front()->A;
	}
}



CGAL_Point_2 Polygon_Segment::end() const
{
	assert(stopped());
	
	return Tr_previous.back()->B;
}



#pragma warning(push)
#pragma warning(disable:4715)
CGAL_Point_2 Polygon_Segment::pt(FT t) const
{
	// We first check if the requested t corresponds to the current interval
	if (Tr != nullptr && Tr->t_int_start <= t && t <= Tr->t_int_end) {
		return Tr->A + (t - Tr->t_int_start) * Tr->dA;
	} 
	
	// If it not the case, then we loop on all the previous intervals 
	// until finding the one that contains t 
	else {
		assert(Tr_previous.size() > 0);
		assert(t >= Tr_previous.front()->t_int_start);

		// Loops on different intervals
		std::list<Segment_Translation*>::const_iterator it;
		for (it = Tr_previous.begin() ; it != Tr_previous.end() ; it++) {
			Segment_Translation* T = (*it);
			if (T->type == PROGRESSIVE) {
				// If the translation is progressive, 
				// i.e. if a polygon progressively collides with another polygon
				if (T->t_int_start <= t && t <= T->t_int_end) {
					return T->A + (t - T->t_int_start) * T->dA;
				}
			} else {
				// If the translation is instantaneous,
				// i.e. if this segment has been created at t = 0,
				// or if an edge of polygon has intersected orthogonally another polygon at t > 0
				// By convention it seems more logical to take the last point of the interval,
				// the one that is the closest to the trajectory of the segment in the next 
				// temporal interval
				if (T->t_int_end == t) return T->B;
			}
		}

		assert(it != Tr_previous.end());
	}
}
#pragma warning(pop)


void Polygon_Segment::update_translation(const FT & t, const CGAL_Point_2 & A_t, const CGAL_Vector_2 & dA_t)
{
	// Indicates that the current translation is no longer valid after time t,
	// and inserts Tr into the list of previous translations
	Tr->set_end(t);
	Tr_previous.push_back(Tr);

	// Creates a new translation
	Tr = new Segment_Translation(t, A_t, dA_t);
}


void Polygon_Segment::stop(const Constraint C, const FT & t)
{
	// Indicates that the current translation is no longer valid after time t,
	// and inserts Tr into the list of previous translations
	Tr->set_end(t);

	Tr_previous.push_back(Tr);

	Tr = nullptr;

	C_stop = C;
	t_stop = t;
}


#if 0
bool Polygon_Segment::checks_if_belongs(const CGAL_Point_2 & V_t, const FT & t)
{
	std::list<Constraint> C_limits;

	return checks_if_belongs(V_t, C_limits, t);
}


bool Polygon_Segment::checks_if_belongs(const CGAL_Point_2 & V_t, const Constraint C, const FT & t)
{
	std::list<Constraint> C_limits;
	C_limits.push_back(C);

	return checks_if_belongs(V_t, C_limits, t);
}


bool Polygon_Segment::checks_if_belongs(const CGAL_Point_2 & V_t, const std::list<Constraint> & C_limits, const FT & t)
{
	// We want to determine if the point V_t, which represents the location of a vertex v at time t, belongs to the segment.
	// We denote by A the origin of the segment and B the position of its other end at time t.
	// There are two main possibilities.

	// Case 1 : V_t is positioned at one end of the segment.
	// This time we make use the list of constraints C_limits,
	// and test if one of them corresponds to C_init or C_stop.

	for (std::list<Constraint>::const_iterator it_c = C_limits.begin() ; it_c != C_limits.end() ; it_c++) {
		Constraint C = (*it_c);

		// If C refers to the same line of the initial or stop constraint of the segment
		// then the inclusion of V_t in [AB] only depends on which side of the line is pointed by C and C_init or C_stop
		if (C.first == C_init.first) {
			return (C.second == C_init.second);
		} else if (C.first == C_stop.first) {
			return (C.second == C_stop.second);
		}
	}

	// Case 2 : V_t belongs to the interior of the segment. 
	// This assertion is easy to check : we numerically test if A.x <= V_t.x <= B.x or A.y <= V_t.y <= B.y

	CGAL_Point_2 A = origin();
	CGAL_Point_2 B = (stopped() ? end() : pt(t));
	CGAL_Segment_2 AB (A, B);

	return AB.has_on(V_t);
}
#endif


// The next function only concern opposite bidirectional segments,
// which are initialized with a C_init = nullptr.

// In order to improve the later performances of the algorithm,
// we would like to find the closest lines to the origins of bidirectional segments,
// and set C_init anyway to a consistent value.



void Polygon_Segment::set_as_opposite_bidirectional_segments(Polygon_Segment* s_1, Polygon_Segment* s_2)
{
	s_1->opposite = s_2;
	s_2->opposite = s_1;
}



bool Polygon_Segment::exists_opposite_segment_without_initinal_constraint() const
{
	// This function is for preventing useless and expensive computations.

	return (opposite != nullptr && opposite->C_init.first == nullptr);
}



Polygon_Segment* Polygon_Segment::get_opposite() const
{
	return opposite;
}



void Polygon_Segment::set_pseudo_init_constraint(const CGAL_Point_2 & A, std::list<Intersection_Line*> & L)
{
	assert(C_init.first == nullptr);

	// We consider two opposite bidirectional segments :
	// - A is the extremum of the first one (this object),
	// - L the list of Intersection_Lines that cut the second.

	// We can set an initial constraint for the first segment,
	// using the line that is the closest to it.
	// To detect this line, we use the coordinates of A.

	FT dist_min = FLT_MAX;
	FT absolute_dist_min = FLT_MAX;
	Intersection_Line* argmin = nullptr;

	for (std::list<Intersection_Line*>::iterator it_l = L.begin() ; it_l != L.end() ; it_l++) {
		Intersection_Line* J = (*it_l);

		const CGAL_Line_2 & J_line = J->line;
		FT dist = J_line.a() * A.x() + J_line.b() * A.y() + J_line.c();
		FT absolute_dist = CGAL::abs(dist);
		if (absolute_dist < absolute_dist_min) {
			dist_min = dist;
			absolute_dist_min = absolute_dist;
			argmin = J;
		}
	}
	
	C_init = Constraint(argmin, dist_min < 0 ? MINUS : PLUS);
}



void Polygon_Segment::set_pseudo_init_constraint(const Constraint C_pseudo_init)
{
	assert(C_init.first == nullptr);
	C_init = C_pseudo_init;
}



Constraint Polygon_Segment::get_pseudo_init_constraint() const
{
	return C_init;
}



void Polygon_Segment::insert_as_crossed_line(Intersection_Line* I)
{
	C_crossed.push_back(I);
}



void Polygon_Segment::insert_as_crossed_lines(const std::list<Intersection_Line*> & L)
{
	std::copy(L.begin(), L.end(), std::back_inserter(C_crossed));
}



bool Polygon_Segment::includes_point_on_support_line(const CGAL_Point_2 & M, const FT & t) const
{
	// We assume that M is a point of the support line of this segment.
	// Here, we determine if at time t, M is contained by this segment,
	// which is the case if A.x <= M.x <= B.x, or A.y <= M.y <= B.y,
	// where A and B represent the segments' ends.

	CGAL_Point_2 A = origin();
	CGAL_Point_2 B = (stopped() ? end() : pt(t));
	
	FT x_a = A.x(), x_b = B.x(), x = M.x();
	if (x_a != x_b) {
		return ((x_a <= x && x <= x_b) || (x_b <= x && x <= x_a));
	} else {
		FT y_a = A.y(), y_b = B.y(), y = M.y();
		return ((y_a <= y && y <= y_b) || (y_b <= y && y <= y_a));
	}
}



bool Polygon_Segment::includes_point_at_intersection(const CGAL_Point_2 & M, const Constraint & C, const FT & t) const
{
	// Like before, we assume that M is a point of the support line of this segment.
	// Here, we know that M is the intersection of the 'support' line and C.first,
	// and we also know on which side (or halfplane) M is located, via C.second

	Intersection_Line* J = C.first;
	Sign J_eps = C.second;

	// If J is crossed by the segment, then M belongs to it
	std::list<Intersection_Line*>::const_iterator it_l = std::find(C_crossed.cbegin(), C_crossed.cend(), J);
	if (it_l != C_crossed.cend()) return true;

	// Compares J to the lines that delimit the segment.
	// If they are not nullptr, then the segment includes M 
	// iff the initial or stop constraint have the same sign as C

	Intersection_Line *I_init = C_init.first, *I_stop = C_stop.first;
	if (I_init == J) {
		return (C_init.second == J_eps);
	} else if (I_stop == J) {
		return (C_stop.second == J_eps);
	}

	if (I_init != nullptr && I_stop != nullptr) {
		// M is necessarily located outside the point because J is not I_init, I_stop or in C_crossed
		return false;

	} else {
		// Unfortunately we have to perform a numerical computation
		// M belongs to the line, so if A and B represent the segments' ends,
		// we should have A.x <= M.x <= B.x (or A.y <= M.y <= B.y)
		return includes_point_on_support_line(M, t);
	}
}



bool Polygon_Segment::includes_point_at_intersection(const CGAL_Point_2 & M, const std::list<Constraint> & C_limits, const FT & t) const
{
	// Same function as before, 
	// for the case when several lines intersect in M

	for (std::list<Constraint>::const_iterator it_c = C_limits.begin() ; it_c != C_limits.end() ; it_c++) {
		Constraint C = (*it_c);

		Intersection_Line* J = C.first;
		Sign J_eps = C.second;

		// If J is crossed by the segment, then M belongs to it
		std::list<Intersection_Line*>::const_iterator it_l = std::find(C_crossed.cbegin(), C_crossed.cend(), J);
		if (it_l != C_crossed.cend()) return true;

		// Compares J to the lines that delimit the segment.
		// If they are not nullptr, then the segment includes M 
		// iff the initial or stop constraint have the same sign as C

		Intersection_Line *I_init = C_init.first, *I_stop = C_stop.first;
		if (I_init == J) {
			return (C_init.second == J_eps);
		} else if (I_stop == J) {
			return (C_stop.second == J_eps);
		}
	}

	// Numerical computation
	return includes_point_on_support_line(M, t);
}



bool Polygon_Segment::includes_edge(const CGAL_Point_2 & V_1, const CGAL_Point_2 & V_2, const Constraint & C_1, const Constraint & C_2,
	const std::list<Intersection_Line*> & CL) const
{
	// This function is called when all segments no longer propagate,
	// as we group adjacent polygons of the set to define the facets of the partition.
	
	// We assume that e = (v1 v2) is delimited by two constraints C_1 and C_2.

	Intersection_Line *I_1 = C_1.first, *I_2 = C_2.first;
	Sign eps_1 = C_1.second, eps_2 = C_2.second;

	// First of all we determine if C_1 and C_2 belong to the list of lines crossed by the segment.
	// If so, we return true.

	if (std::find(C_crossed.cbegin(), C_crossed.cend(), I_1) != C_crossed.cend()) return true;
	if (std::find(C_crossed.cbegin(), C_crossed.cend(), I_2) != C_crossed.cend()) return true;

	// After that we check if I_1 or I_2 are used to define C_init.
	// If so, we compare the sign of the constraint to C_1 and C_2.

	Intersection_Line* I_init = C_init.first;
	if (I_init == I_1) {
		return (C_init.second == eps_1);
	} else if (I_init == I_2) {
		return (C_init.second == eps_2);
	}

	// Same for C_stop.

	Intersection_Line* I_stop = C_stop.first;
	if (I_stop == I_1) {
		return (C_stop.second == eps_1);
	} else if (I_stop == I_2) {
		return (C_stop.second == eps_2);
	}

	// Finally, we focus on lines which are concurrent with this->support, I_1 and I_2.
	// If such lines don't exist, then we return false because we have already listed all the possible cases.

	if (CL.empty()) {
		return false;
	} else {
		if (std::find(CL.begin(), CL.end(), I_init) != CL.end() && std::find(CL.begin(), CL.end(), I_stop) != CL.end()) {
			// Arithmetic computation
			CGAL_Point_2 A = origin(), B = end();
			CGAL_Point_2 M = CGAL::midpoint(V_1, V_2);

			FT x_a = A.x(), x_b = B.x(), x = M.x();
			if (x_a != x_b) {
				return ((x_a <= x && x <= x_b) || (x_b <= x && x <= x_a));
			} else {
				FT y_a = A.y(), y_b = B.y(), y = M.y();
				return ((y_a <= y && y <= y_b) || (y_b <= y && y <= y_a));
			}

		} else {
			return false;
		}
	}
}


std::list<Intersection_Line*>::const_iterator Polygon_Segment::crossed_lines_begin() const
{
	return C_crossed.cbegin();
}


std::list<Intersection_Line*>::const_iterator Polygon_Segment::crossed_lines_end() const
{
	return C_crossed.cend();
}



void Polygon_Segment::check() const 
{
	// Assuming the segment is stopped, performs some checks

	CGAL_Point_2 A = origin();
	CGAL_Point_2 B = end();

	// Test 1 : C_init and C_stop are not null
	Intersection_Line* I_init = C_init.first;
	Intersection_Line* I_stop = C_stop.first;
	Intersection_Line* I_supp = support;

	bool test_1 = (I_init != nullptr && I_stop != nullptr);
	if (!test_1) who_is();
	assert(test_1);

	CGAL_Point_2 A_ref, B_ref;
	bool A_ref_exists = false, B_ref_exists = false;
	CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Line_2)>::type object_1 = CGAL::intersection(I_supp->line, I_init->line);
	if (object_1) {
		if (const CGAL_Point_2* ptr = boost::get<CGAL_Point_2>(&*object_1)) {
			A_ref = *ptr;
			A_ref_exists = true;
		}
	}
	CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Line_2)>::type object_2 = CGAL::intersection(I_supp->line, I_stop->line);
	if (object_2) {
		if (const CGAL_Point_2* ptr = boost::get<CGAL_Point_2>(&*object_2)) {
			B_ref = *ptr;
			B_ref_exists = true;
		}
	}

	assert(A_ref_exists && B_ref_exists);

	// Test 2 : B is the intersection of I_supp and I_stop
	bool test_2 = (B == B_ref);
	if (!test_2) {
		who_is();
		std::cout << "B = " << B << std::endl;
		std::cout << "B_ref = " << B_ref << std::endl;
	}
	assert(test_2);

	bool test_3 = true, test_4 = true, test_5 = true;

	if (opposite == nullptr) {
		// Test 3A : If there exists no segment going in the opposite direction, A is the intersection of I_supp and I_init
		test_3 = (A == A_ref);
		if (!test_3) {
			who_is();
			std::cout << "A = " << A << std::endl;
			std::cout << "A_ref = " << A_ref << std::endl;
		}
		assert(test_3);
	}

	else {
		// Test 3B : If there exists a segment going in the opposite direction, then A is included between
		// the intersection of I_supp and I_init, and all lines of C_crossed (or I_stop)
		for (std::list<Intersection_Line*>::const_iterator it_l = C_crossed.begin() ; it_l != C_crossed.end() ; it_l++) {
			CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Line_2)>::type object_3 = CGAL::intersection(I_supp->line, (*it_l)->line);
			if (object_3) {
				if (const CGAL_Point_2* ptr = boost::get<CGAL_Point_2>(&*object_3)) {
					CGAL_Point_2 K = *ptr;
					CGAL_Segment_2 AK_ref(A_ref, K);
					test_4 = (AK_ref.has_on(A));
					if (!test_4) {
						who_is();
						std::cout << "AK = " << AK_ref << std::endl;
						std::cout << "A = " << A << std::endl;
					}
					assert(test_4);
				}
			}
		}

		CGAL_Segment_2 AB_ref(A_ref, B_ref);
		test_5 = (AB_ref.has_on(A));
		if (!test_5) {
			who_is();
			std::cout << "AB = " << AB_ref << std::endl;
			std::cout << "A = " << A << std::endl;
		}
		assert(test_5);
	}
}



void Polygon_Segment::who_is() const
{
	std::cout << "*** Segment at fault : " << std::endl;
	std::cout << "id : " << id_object << " [" << id_plane << "] on line " << support->id_object << " which has " << support->planes.size() << " planes : " << std::endl;
	for (int i : support->planes) {
		std::cout << "# " << i << std::endl;
	}
	std::cout << "t_init = " << t_init << std::endl;
	std::cout << "t_stop = " << t_stop << std::endl;
	std::cout << "opposite = " << opposite << std::endl;
	if (C_init.first == nullptr) {
		std::cout << "C_init = nullptr" << std::endl;
	} else {
		std::cout << "C_init = " << C_init.first->id_object << ", " << (C_init.second == 1 ? "PLUS" : "MINUS") << std::endl;
	}
	if (C_stop.first == nullptr) {
		std::cout << "C_stop = nullptr" << std::endl;
	} else {
		std::cout << "C_stop = " << C_stop.first->id_object << ", " << (C_stop.second == 1 ? "PLUS" : "MINUS") << std::endl;
	}
	std::cout << "|C_crossed| = " << C_crossed.size() << std::endl;
	std::cout << "|Tr| = " << Tr_previous.size() << std::endl;
}
}