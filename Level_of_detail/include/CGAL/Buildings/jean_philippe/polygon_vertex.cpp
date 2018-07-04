#include "support_plane_objects.h"
#include "support_plane.h"
#include "universe.h"
#include "parameters.h"
#include "event_queue.h"
#include "stats.h"

namespace JPTD {
using CGAL::to_double;



Polygon_Vertex::Polygon_Vertex(const int _id_plane, const FT & _t_init, const CGAL_Point_2 & _M, const CGAL_Vector_2 & _dM, Event_Flags flags, bool is_moving)
	: Support_Plane_Object(_id_plane)
{
	Universe::map_of_planes[_id_plane]->vertices[id_object] = this;

	M = _M;
	dM = _dM;
	t_init = _t_init;
	t_stop = _t_init;

	k_th_interval = 0;
	queued_events = 0;

	if (is_moving) {
		// Sets the lower and upper bounds for the time tau,
		// which is the time required by this vertex to move by a distance D
		const FT dm_sq = dM.squared_length();
		tau_inf = FT(sqrt(CGAL::to_double(Universe::params->D_inf_2 / dm_sq)));
		tau_sup = FT(sqrt(CGAL::to_double(Universe::params->D_sup_2 / dm_sq)));

		// Adds this vertex to the table of moving objects.
		Universe::map_of_objects[id_object] = this;
		++Universe::moving_objects;
	} else {
		tau_inf = tau_sup = 0;
	}

	e1 = e2 = nullptr;
	crossed_lines = 0;

	is_active = is_moving;
	
	constraints = std::list<Constraint>();

	guided_segments = std::list<Polygon_Segment*>();
	polygon = nullptr;

	paired_vertex = Companion();

	if (flags & SCHEDULE) schedule_events();
}



Polygon_Vertex::Polygon_Vertex(const int _id_plane, const FT & _t_init, const CGAL_Point_2 & _M, const CGAL_Vector_2 & _dM, const Constraint & C, Intersection_Line* I_discarded, Event_Flags flags)
	: Polygon_Vertex(_id_plane, _t_init, _M, _dM, NO_SCHEDULE, true)
{
	assert(C.second != ZERO);
	constraints.push_back(C);

	if (flags & SCHEDULE) schedule_events(I_discarded);
}



Polygon_Vertex::Polygon_Vertex(const int _id_plane, const FT & _t_init, const CGAL_Point_2 & _M, const CGAL_Vector_2 & _dM, const Constraint & C, const std::list<Intersection_Line*> & I_discarded, Event_Flags flags)
	: Polygon_Vertex(_id_plane, _t_init, _M, _dM, NO_SCHEDULE, true)
{
	assert(C.second != ZERO);
	constraints.push_back(C);

	if (flags & SCHEDULE) schedule_events(I_discarded);
}



Polygon_Vertex::Polygon_Vertex(const int _id_plane, const FT & _t_init, const CGAL_Point_2 & _M, const Constraint & C_1, const Constraint & C_2)
	: Polygon_Vertex(_id_plane, _t_init, _M, CGAL_Vector_2(0, 0), NO_SCHEDULE, false)
{
	assert(C_1.second != ZERO);
	assert(C_2.second != ZERO);

	constraints.push_back(C_1);
	constraints.push_back(C_2);
}



Polygon_Vertex::Polygon_Vertex(Polygon_Vertex* v_ts, Event_Flags flags)
	: Polygon_Vertex(v_ts->id_plane, v_ts->t_init, v_ts->M, v_ts->dM, NO_SCHEDULE, v_ts->is_active)
{
	// Constructs a vertex v_os which is the same vertex as v_ts,
	// except that it is on the other side of line I = C_ts.first

	Constraint C_ts = v_ts->get_constraint();
	Constraint C_os = Constraint(C_ts.first, (C_ts.second == PLUS ? MINUS : PLUS));
	constraints.push_back(C_os);

	set_paired_vertices(v_ts, this);

	// The value of the flags will not differ from the moment when v_ts has been itself created
	// So, if intersections are computed, we first copy t_border and I_border from v_ts,
	// and we duplicate the events involving v_ts.

	if (flags & SCHEDULE) schedule_events(v_ts);

	// Guided segments and owning polygon aren't known so far
}



Polygon_Vertex::~Polygon_Vertex()
{
	if (is_active) {
		++KP_Stats::life_expectancy_terms;
		KP_Stats::life_expectancy_lines += crossed_lines;
		//CGAL_Vector_2 D = (t_stop - t_init) * dM;
		//KP_Stats::life_expectancy_distance += sqrt(to_double(D.squared_length()));

		++KP_Stats::schedule_events_vertices;
		KP_Stats::schedule_events_lines += lines_to_events.size();

		--Universe::moving_objects;
	}

	reset_paired_vertex();

	deschedule_events();

	std::map<int, Polygon_Vertex*>::iterator it_v = Universe::map_of_planes[id_plane]->vertices.find(id_object);
	Universe::map_of_planes[id_plane]->vertices.erase(it_v);

	if (is_active) {
		it_v = Universe::map_of_objects.find(id_object);
		Universe::map_of_objects.erase(it_v);
	}
}



CGAL_Point_2 Polygon_Vertex::pt(const FT & t) const
{
	return M + t * dM;
}



void Polygon_Vertex::add(Polygon_Edge* e)
{
	if (e1 == nullptr) {
		e1 = e;
	} else {
		assert(e2 == nullptr);
		e2 = e;
	}
}



void Polygon_Vertex::remove(Polygon_Edge* e)
{
	if (e1 == e) {
		e1 = nullptr;
	} else {
		assert(e2 != nullptr);
		e2 = nullptr;
	}
}



bool Polygon_Vertex::unconstrained() const
{
	return constraints.empty();
}



bool Polygon_Vertex::is_constrained_by(Intersection_Line* I) const
{
	for (std::list<Constraint>::const_iterator it_l = constraints.begin() ; it_l != constraints.end() ; it_l++) {
		if (it_l->first == I) return true;
	}
	return false;
}



Constraint Polygon_Vertex::get_constraint() const
{
	assert(!constraints.empty());
	return constraints.front();
}



Constraint Polygon_Vertex::get_second_constraint() const
{
	assert(!is_active);
	return constraints.back();
}



Sign Polygon_Vertex::sign_of_constraint(Intersection_Line* I) const
{
	assert(!constraints.empty());
	
	std::list<Constraint>::const_iterator it_l = constraints.begin();
	do {
		if (it_l->first == I) break;
	} while (++it_l != constraints.end());
	assert(it_l != constraints.end());

	return it_l->second;
}



Polygon_Vertex* Polygon_Vertex::get_constrained_neighbor(Intersection_Line* I_0, const std::list<Intersection_Line*> & I) const
{
	// This function is called as we try to factorize events of the queue.
	// Typically, this vertex intersects a set of lines I and we try to determine
	// if it meets a neighbor, constrained by one of the lines of I, at the same time.

	// I_0 represents the possible constraint of v. 
	// If e1 or e2, the edges of v, lie on the line I_0 then we already know the vertices at their ends
	// are constrained by I_0, so it's not necessary to perform the test with elements of I.

	if (!e1->is_constrained_by(I_0)) {
		Polygon_Vertex* v1 = (e1->v1 == this ? e1->v2 : e1->v1);
		if (v1->is_constrained_neighbor(I)) return v1;
	}

	if (!e2->is_constrained_by(I_0)) {
		Polygon_Vertex* v2 = (e2->v1 == this ? e2->v2 : e2->v1);
		if (v2->is_constrained_neighbor(I)) return v2;
	}

	return nullptr;
}



bool Polygon_Vertex::is_constrained_neighbor(const std::list<Intersection_Line*> & I) const
{
	for (std::list<Intersection_Line*>::const_iterator it_l = I.begin() ; it_l != I.end() ; it_l++) {
		if (is_constrained_by(*it_l)) {
			return true;
		}
	}

	return false;
}



void Polygon_Vertex::stop(const FT & _t_stop)
{
	t_stop = _t_stop;
}



Polygon_Vertex* Polygon_Vertex::get_neighbor_intersecting_identical_line(Intersection_Line* I_0, const std::list<Intersection_Line*> & I, const FT t) const
{
	// This function is called as we try to factorize events of the queue.
	// Typically this vertex intersects a set of lines I, and here we try to determine
	// if one of the neighbors v1 and v2 (ends of e1 and e2) intersect one of the lines in I.

	// Assumption : v1 and v2 are NOT constrained by any of the lines of I,
	// and therefore, edges e1 and e2 are not of length zero.

	if (!e1->is_constrained_by(I_0)) {
		Polygon_Vertex* v1 = (e1->v1 == this ? e1->v2 : e1->v1);
		CGAL_Point_2 M = v1->pt(t);

		for (std::list<Intersection_Line*>::const_iterator it_l = I.begin() ; it_l != I.end() ; it_l++) {
			if ((*it_l)->includes(M)) return v1;
		}
	}

	if (!e2->is_constrained_by(I_0)) {
		Polygon_Vertex* v2 = (e2->v1 == this ? e2->v2 : e2->v1);
		CGAL_Point_2 M = v2->pt(t);
		
		for (std::list<Intersection_Line*>::const_iterator it_l = I.begin() ; it_l != I.end() ; it_l++) {
			if ((*it_l)->includes(M)) return v2;
		}
	}

	return nullptr;
}



void Polygon_Vertex::transfer_segments(Polygon_Vertex* v_dest, const FT & t) 
{
	// Gets the 3D position (a) and the speed (b - a) of the vertex v_dest, which has just been created
	const CGAL_Point_2 V_t = v_dest->pt(t);
	const CGAL_Point_2 V_u = V_t + v_dest->dM;

	const CGAL_Point_3 a = Universe::map_of_planes[v_dest->id_plane]->backproject(V_t);
	const CGAL_Point_3 b = Universe::map_of_planes[v_dest->id_plane]->backproject(V_u);

	for (std::list<Polygon_Segment*>::iterator it_s = guided_segments.begin() ; it_s != guided_segments.end() ; it_s++) {
		Polygon_Segment* s = (*it_s);

		Support_Plane* SP = Universe::map_of_planes[s->id_plane];

		// Updates the speed of the segment
		const CGAL_Point_2 A = SP->project(a);
		const CGAL_Point_2 B = SP->project(b);
		s->update_translation(t, A, B - A);

		// From now on, s is now going to move along with v_dest
		v_dest->set_as_guided_segment(s);
	}

	guided_segments.clear();
}



void Polygon_Vertex::stop_segments(Intersection_Line* L_ik, const FT & t)
{
	// Segments carried by this vertex no longer propagate after time t, when v(t) intersects a line I.
	// This vertex v belongs to plane i, and I represents the plane k : it is denoted as L_ik.
	// A segment lives on plane j, on line L_ji. We compute the its constraint C_jk.
	
	const int i = id_plane;
	const int k = L_ik->planes.front();

	for (std::list<Polygon_Segment*>::iterator it_s = guided_segments.begin() ; it_s != guided_segments.end() ; it_s++) {
		Polygon_Segment* s = (*it_s);

		const int j = s->id_plane;
		Support_Plane* P_j = Universe::map_of_planes[j];
		
		Intersection_Line* L_jk = P_j->get_line(k);
		Sign eps_jk = L_jk->sign(s->origin());
		assert(eps_jk != ZERO);

		Constraint C_jk (L_jk, eps_jk);
		s->stop(C_jk, t);
		
		// If the segment is bidirectional, 
		// the stop constraint for s may be used for delimiting the opposite segment.

		if (s->exists_opposite_segment_without_initinal_constraint()) {
			Polygon_Segment* s_opposite = s->get_opposite();
			s_opposite->set_pseudo_init_constraint(C_jk);
		}
	}

	guided_segments.clear();
}



void Polygon_Vertex::extend_segments(Intersection_Line* L_ik) const
{
	// This is the opposite function to before.
	// A constraint vertex v intersects a line I but it can propagate beyond that line,
	// so all segments carried by this segment also keep propagating.
	// This vertex v belongs to plane i, and I represents the plane k : it is denoted as L_ik.

	const int i = id_plane;

	for (std::list<Polygon_Segment*>::const_iterator it_s = guided_segments.begin() ; it_s != guided_segments.end() ; it_s++) {
		Polygon_Segment* s = (*it_s);

		const int j = s->id_plane;
		Support_Plane* P_j = Universe::map_of_planes[j];
		
		// We always insert the line that corresponds to the plane k in the list of lines intersected by s

		std::set<Intersection_Line*> _CL;
		for (std::list<int>::iterator it_p = L_ik->planes.begin() ; it_p != L_ik->planes.end() ; it_p++) {
			int k = L_ik->planes.front();
			Intersection_Line* L_jk = P_j->get_line(k);
			_CL.insert(L_jk);
		}
		
		std::list<Intersection_Line*> CL(_CL.begin(), _CL.end());
		for (std::list<Intersection_Line*>::iterator it_l = CL.begin() ; it_l != CL.end() ; it_l++) {
			s->insert_as_crossed_line(*it_l);
		}
		Intersection_Line* L_jk = CL.front();

		// If the segment is bidirectional,
		// L_jk may be used to delimit the opposite segment and set a pseudo-init constraint.
		if (s->exists_opposite_segment_without_initinal_constraint()) {
			Polygon_Segment* s_opposite = s->get_opposite();
			Sign eps_jk = L_jk->sign(s_opposite->origin());
			Constraint C_jk (L_jk, eps_jk);

			s_opposite->set_pseudo_init_constraint(C_jk);
		}
	}
}



void Polygon_Vertex::indicate_line_initially_crossed_by_segments(Intersection_Line* L_ik) const
{
	// This function is called during the initialization process, and is very similar to before, 
	// except that we don't try to set an initial constraint to the unidirectional segments
	// because such a constraint already exists.

	const int i = id_plane;

	for (std::list<Polygon_Segment*>::const_iterator it_s = guided_segments.begin() ; it_s != guided_segments.end() ; it_s++) {
		Polygon_Segment* s = (*it_s);

		const int j = s->id_plane;
		Support_Plane* P_j = Universe::map_of_planes[j];
		
		// We always insert the line that corresponds to the plane k in the list of lines intersected by s

		std::set<Intersection_Line*> _CL;
		for (std::list<int>::iterator it_p = L_ik->planes.begin() ; it_p != L_ik->planes.end() ; it_p++) {
			int k = L_ik->planes.front();
			Intersection_Line* L_jk = P_j->get_line(k);
			_CL.insert(L_jk);
		}
		
		std::list<Intersection_Line*> CL(_CL.begin(), _CL.end());
		for (std::list<Intersection_Line*>::iterator it_l = CL.begin() ; it_l != CL.end() ; it_l++) {
			s->insert_as_crossed_line(*it_l);
		}
		Intersection_Line* L_jk = CL.front();
	}
}



void Polygon_Vertex::copy_crossed_lines(Polygon_Vertex* v_os) const
{
	// This function is called during the initialization process.
	// We suppose that we know the lines intersected by segments 
	// carried by a vertex that is on one side of a line.

	// Now, we want to copy these informations to segments carried by the same vertex,
	// but that is located on the other side of the same line.

	std::list<Polygon_Segment*>::const_iterator it_s1 = guided_segments.begin();
	std::list<Polygon_Segment*>::iterator it_s2 = v_os->guided_segments.begin();

	while (it_s1 != guided_segments.end() && it_s2 != v_os->guided_segments.end()) {
		Polygon_Segment* s1 = (*it_s1);
		Polygon_Segment* s2 = (*it_s2);
		
		for (std::list<Intersection_Line*>::const_iterator it_l = s1->crossed_lines_begin() ; it_l != s1->crossed_lines_end() ; ++it_l) {
			s2->insert_as_crossed_line(*it_l);
		}

		++it_s1;
		++it_s2;
	}
}



void Polygon_Vertex::set_as_guided_segment(Polygon_Segment* s)
{
	guided_segments.push_back(s);
}



void Polygon_Vertex::set_as_member_of_polygon(Polygon* _polygon)
{
	polygon = _polygon;
}



Polygon* Polygon_Vertex::get_polygon()
{
	return polygon;
}



void Polygon_Vertex::schedule_events()
{
	std::list<Intersection_Line*> discarded;
	schedule_events(discarded);
}



void Polygon_Vertex::schedule_events(Intersection_Line* I_redirect)
{
	std::list<Intersection_Line*> discarded(1, I_redirect);
	schedule_events(discarded);
}



void Polygon_Vertex::schedule_events(const std::list<Intersection_Line*> & discarded)
{
	// Each vertex has a map that associates lines to events.
	// A null event indicates that the events has already been processed, or cannot happen.
	// For this reason, we insert associate null events to the lines that must be discarded
	// from the algorithm.

	if (!discarded.empty()) {	
		for (std::list<Intersection_Line*>::const_iterator it_l = discarded.begin() ; it_l != discarded.end() ; it_l++) {
			Intersection_Line* I = (*it_l);
			if (I != nullptr) lines_to_events[I->id_object] = nullptr;
		}
	}

	reschedule_events();
}



void Polygon_Vertex::reschedule_events()
{
	// We assume that this vertex is active and there are no events associated to this vertex in the queue.
	// If it is not the case, returns.

	if (queued_events != 0) return;

	// We assume there are still lines that it can intersect.
	// In this function, we are going to compute a set of events occuring in I_k = [tk_min, tk_max],
	// where tk_min = t_init + k * tau, tk_max = t_init + (k + 1) * tau.

	// The variable tau is the time when the vertex moves by a distance of D (a parameter).
	// As it is impossible to compute it exactly, we use two bounds tau_inf and tau_sup
	// to set t_min and t_max.

	// Given our initial assumption, we expect to compute at least one event that we insert in E_VL.
	// We iteratively loop on intervals I_k until E_VL gets no longer empty.

	Support_Plane* SP = Universe::map_of_planes[id_plane];

	std::list<Event_Vertex_Line*> E_VL;
	while (E_VL.empty()) {

		// We define the current interval I_k.

		FT t_min = t_init + k_th_interval * tau_inf;
		FT t_max = t_init + (k_th_interval + 1) * tau_sup;

		// We compute v->pt(t_min) and v->pt(t_max).
		// We define J, a list of lines that might have been crossed in the meanwhile.

		std::list<std::tuple<Intersection_Line*, bool, FT> > J;
		SP->search_lines_in_neighborhood(this, t_min, t_max, J);

		for (std::list<std::tuple<Intersection_Line*, bool, FT> >::iterator it_l = J.begin(); it_l != J.end(); it_l++) {
			const std::tuple<Intersection_Line*, bool, FT> & J_tuple = (*it_l);

			Intersection_Line* I = std::get<0>(J_tuple);
			bool t_already_known = std::get<1>(J_tuple);
			Event_Vertex_Line* e_vl = nullptr;

			// For each line I, we query the look-up table lines_to_events.
			// The non-existence of an entry implies that the line is not discarded
			// or that no intersection with this line has been computed before.
			// When we compute an event, we increment by anticipation the number
			// of queued events involving this vertex.

			if (lines_to_events.find(I->id_object) == lines_to_events.end()) {
				if (t_already_known) {
					e_vl = make_event(I, std::get<2>(J_tuple));
				} else {
					e_vl = make_event(I);
				}
				lines_to_events[I->id_object] = e_vl;

				if (e_vl != nullptr) {
					// e_vl can be null is the intersection takes place behind the vertex
					E_VL.push_back(e_vl);
					++queued_events;
				}
			}
		}

		// Iterates on k.
		++k_th_interval;
	}

	// Inserts events in the map.
	Universe::event_queue->push(E_VL);
}



void Polygon_Vertex::decrement_queued_events(const int n)
{
	queued_events -= n;
}



Event_Vertex_Line* Polygon_Vertex::get_event_for_line(const int I_object) const
{
	// This function is called to access an event that has been computed before

	std::map<int, Event_Vertex_Line*>::const_iterator it_e = lines_to_events.find(I_object);
	assert(it_e != lines_to_events.end());

	return it_e->second;
}



void Polygon_Vertex::delete_event(Event_Vertex_Line* e_vl)
{
	// This function is called to delete a specified event, typically when it is processed,
	// and set the local reference to such event to nullptr.

	std::map<int, Event_Vertex_Line*>::iterator it_e = lines_to_events.find(e_vl->intersected);
	assert(it_e != lines_to_events.end() && it_e->second == e_vl);
	
	delete e_vl;
	it_e->second = nullptr;
}



void Polygon_Vertex::delete_events(const std::list<Event_Vertex_Line*> & E_VL)
{
	for (std::list<Event_Vertex_Line*>::const_iterator it_e = E_VL.begin() ; it_e != E_VL.end() ; it_e++) {
		delete_event(*it_e);
	}
}



void Polygon_Vertex::copy_events(Polygon_Vertex* v_ts)
{
	// Copies local events that were not existing before

	if (lines_to_events.size() != v_ts->lines_to_events.size()) {
		std::map<int, Event_Vertex_Line*>::iterator it_e_ts;

		for (it_e_ts = v_ts->lines_to_events.begin() ; it_e_ts != v_ts->lines_to_events.end() ; it_e_ts++) {
			int I_id_object = it_e_ts->first;
			Event_Vertex_Line* e_ts = it_e_ts->second;

			// Missing event : makes a local copy if e_ts is not nullptr
			if (lines_to_events.find(I_id_object) == lines_to_events.end()) {
				if (e_ts == nullptr) {
					lines_to_events[I_id_object] = nullptr;
				} else {
					lines_to_events[I_id_object] = new Event_Vertex_Line(id_object, I_id_object, e_ts->t_intersectant);
				}
			}
		}
	}
}



void Polygon_Vertex::schedule_events(Polygon_Vertex* v_ts)
{
	// We assume that some events have been computed for the master vertex v_ts.
	// We just read the look-up table of events and make a deep copy of all vertices.

	for (std::map<int, Event_Vertex_Line*>::iterator it_e = v_ts->lines_to_events.begin() ; it_e != v_ts->lines_to_events.end() ; ++it_e) {
		int index = it_e->first;
		Event_Vertex_Line* e_ts = it_e->second;

		if (lines_to_events.find(index) == lines_to_events.end()) {
			if (e_ts == nullptr) {
				lines_to_events[index] = nullptr;
			} else {
				Event_Vertex_Line* e_os = new Event_Vertex_Line(id_object, e_ts->intersected, e_ts->t_intersectant);
				lines_to_events[index] = e_os;
				Universe::event_queue->push(e_os, e_ts->queue_iterator);
			}
		}
	}

	queued_events = v_ts->queued_events;
}



void Polygon_Vertex::deschedule_events()
{
	// Removes events from the queue before destroying them.

	for (std::map<int, Event_Vertex_Line*>::iterator it_e = lines_to_events.begin() ; it_e != lines_to_events.end() ; it_e++) {
		Event_Vertex_Line* e_vl = it_e->second;
		if (e_vl != nullptr) {
			Universe::event_queue->erase(e_vl);
			delete it_e->second;
		}
	}

	queued_events = 0;
}



std::pair<FT, bool> Polygon_Vertex::get_intersection_time(Intersection_Line* I) const
{
	// We assume that v is active.
	// We are going to determine the time t when the vertex v(t) = M + t * dM intersects I.

	// When the direction of the line is collinear to the direction of the vertex, 
	// there is no intersection and we return FLT_MAX.

	const FT & a = I->a(), b = I->b(), c = I->c();

	FT t_den = a * dM.x() + b * dM.y();
	if (t_den == 0) {
		// The direction of the vector is collinear to the line
		return std::make_pair(0, false);
	}

	FT t_num = -(a * M.x() + b * M.y() + c);
	
	FT t = t_num / t_den;
	return std::make_pair(t, true);
}



Event_Vertex_Line* Polygon_Vertex::make_event(Intersection_Line* I) const
{
	std::pair<FT, bool> R = get_intersection_time(I);

	if (!R.second || R.first < t_init) {
		return nullptr;
	} else {
		return new Event_Vertex_Line(id_object, I->id_object, R.first);
	}
}



Event_Vertex_Line* Polygon_Vertex::make_event(Intersection_Line* I, const FT & t) const
{
	if (t < t_init) {
		return nullptr;
	} else {
		return new Event_Vertex_Line(id_object, I->id_object, t);
	}
}



bool Polygon_Vertex::represents_same_intersection(Polygon_Vertex* v1, Polygon_Vertex* v2) const
{
	return (represents_same_intersection(v1) || represents_same_intersection(v2));
}



bool Polygon_Vertex::represents_same_intersection(Polygon_Vertex* v) const
{
	// We assume that all vertices are stopped. Therefore they have two constraints.
	const Constraint C_1 = get_constraint(), C_2 = get_second_constraint();
	const Constraint C_v1 = v->get_constraint(), C_v2 = v->get_second_constraint();

	// We get the indices of the lines that correspond to such constraints.
	const Intersection_Line *I_1 = C_1.first, *I_2 = C_2.first;
	const Intersection_Line *I_v1 = C_v1.first, *I_v2 = C_v2.first;

	// If (C_1, C_2) and (C_v1, C_v2) represent the same lines, then this vertex and v represent the same intersection
	if ((I_1 == I_v1 && I_2 == I_v2) || (I_1 == I_v2 && I_2 == I_v1)) return true;

	// If we couldn't simultaneously match (C_1, C_2) and (C_v1, C_v2), 
	// then we get triplets of concurrent lines in which we find (C_1 and C_2)
	// If C_v1 or C_v2 is part of such a triplet, then we return true.
	std::list<Intersection_Line *> I_L;
	Universe::map_of_planes[v->id_plane]->get_concurrent_lines(I_1, I_2, I_L);

	for (Intersection_Line* L : I_L) {
		if (I_v1 == L || I_v2 == L) return true;
	}

	// The search failed
	return false;
}



void Polygon_Vertex::set_paired_vertices(Polygon_Vertex* v_ts, Polygon_Vertex* v_os)
{
	// We create a link between two vertices v_ts and v_os.
	// v_os is considered as a copy of v_ts, located on the other side of a line I.

	// In this function, we provide the algorithm the information according to what
	// events of v_os don't have to be computed : they are copied from v_ts as well
	// since v_ts and v_os propagate at the same speed.

	// The duplication will take place at the insertion of events (v_ts, I) in
	// the queue and matrices of events.

	v_ts->set_paired_vertex(Companion(v_os, true));  // v_ts is the master
	v_os->set_paired_vertex(Companion(v_ts, false)); // v_os is the slave
}



bool Polygon_Vertex::are_paired_vertices(Polygon_Vertex* v_ts, Polygon_Vertex* v_os)
{
	if (!v_ts->has_paired_vertex() || !v_os->has_paired_vertex()) return false;

	return (v_ts->paired_vertex.first == v_os && v_os->paired_vertex.first == v_ts);
}



void Polygon_Vertex::set_paired_vertex(const Companion & _paired_vertex)
{
	paired_vertex = _paired_vertex;
}



void Polygon_Vertex::reset_paired_vertex()
{
	Polygon_Vertex* v_os = paired_vertex.first;
	if (v_os != nullptr) {
		v_os->set_paired_vertex(Companion());
		paired_vertex = Companion();
	}
}



bool Polygon_Vertex::has_paired_vertex() const
{
	return (paired_vertex.first != nullptr);
}



bool Polygon_Vertex::is_independent() const
{
	return (paired_vertex.first == nullptr || paired_vertex.second);
}



Polygon_Vertex* Polygon_Vertex::get_master_vertex() const
{
	return paired_vertex.first;
}
}