#pragma once
#include "defs.h"
#include "event.h"
#include "support_plane_object.h"
#include "segment_translation.h"
#include <list>
#include <map>
#include <vector>

namespace JPTD {

class Intersection_Line;
class Polygon_Vertex;
class Polygon_Edge;
class Polygon;
class Polygon_Tree;
class Polygon_Group;
class Segment;
class Planar_Segment;
class Polygon_Segment;

typedef std::pair<Intersection_Line*, Sign> Constraint;
typedef std::pair<Polygon_Vertex*, bool> Companion;


class Intersection_Line : public Support_Plane_Object {
public:
	Intersection_Line(const int _id_plane, const CGAL_Line_2 & _line, int intersected);

	~Intersection_Line();

	void mark_as_intersected(int intersected);

	bool intersects(const int id_plane) const;

	void set_inside(const bool _is_inside);

	Sign sign(const CGAL_Point_2 & pt) const;

	Sign sign(Polygon_Vertex* v, const FT t) const;


	bool includes(const CGAL_Point_2 & M) const;

	bool is_parallel(Intersection_Line* I) const;



	bool exist_segments_including_point_outside_intersections(const CGAL_Point_2 & V_t, const FT & t) const;

	bool exist_segments_including_point_at_intersection(const CGAL_Point_2 & V_t, const Constraint & C, const FT & t) const;

	bool exist_segments_including_point_at_intersection(const CGAL_Point_2 & V_t, const std::list<Constraint> & C_limits, const FT & t) const;

	bool exists_segment_adjacent_to_edge(Polygon_Edge* e) const;


	const FT & a() const;

	const FT & b() const;

	const FT & c() const;

public:
	bool is_border;
	bool is_inside;

	const CGAL_Line_2 line;
	const FT _a;
	const FT _b;
	const FT _c;

	const double hint_a;
	const double hint_b;
	const double hint_c;

	std::list<int> planes;
	std::list<Segment*> segments;
};



class Polygon_Vertex : public Support_Plane_Object {
public:
	Polygon_Vertex(const int _id_plane, const FT & _t_init, const CGAL_Point_2 & _M, const CGAL_Vector_2 & _dM, const int _K, Event_Flags flags, bool is_moving);

	Polygon_Vertex(const int _id_plane, const FT & _t_init, const CGAL_Point_2 & _M, const CGAL_Vector_2 & _dM, const Constraint & C, Intersection_Line* I_discarded, const int _K, Event_Flags flags);

	Polygon_Vertex(const int _id_plane, const FT & _t_init, const CGAL_Point_2 & _M, const CGAL_Vector_2 & _dM, const Constraint & C, const std::list<Intersection_Line*> & I_discarded, const int _K, Event_Flags flags);

	Polygon_Vertex(const int _id_plane, const FT & _t_init, const CGAL_Point_2 & _M, const Constraint & C_1, const Constraint & C_2);

	Polygon_Vertex(Polygon_Vertex* v_ts, Event_Flags flags);

	~Polygon_Vertex();


	void add(Polygon_Edge* e);

	void remove(Polygon_Edge* e);

	CGAL_Point_2 pt(const FT & t) const;


	bool unconstrained() const;

	bool is_constrained_by(Intersection_Line* I) const;

	Constraint get_constraint() const;

	Constraint get_second_constraint() const;

	Sign sign_of_constraint(Intersection_Line* I) const;

	Polygon_Vertex* get_constrained_neighbor(Intersection_Line* I_0, const std::list<Intersection_Line*> & I) const;

protected:
	bool is_constrained_neighbor(const std::list<Intersection_Line*> & I) const;

public:
	Polygon_Vertex* get_neighbor_intersecting_identical_line(Intersection_Line* I_0, const std::list<Intersection_Line*> & I, const FT t) const;

public:
	void set_as_guided_segment(Polygon_Segment* s);

	void set_as_member_of_polygon(Polygon* _polygon);

	Polygon* get_polygon();

	void stop(const FT & t_stop);

	void transfer_segments(Polygon_Vertex* v_dest, const FT & t);

	void stop_segments(Intersection_Line* I, const FT & t);

	void extend_segments(Intersection_Line* I) const;


	void indicate_line_initially_crossed_by_segments(Intersection_Line* I) const;

	void copy_crossed_lines(Polygon_Vertex* v_os) const;


	void schedule_events(Polygon_Vertex* v_ts);

	void schedule_events();

	void schedule_events(Intersection_Line* I_redirect);

	void schedule_events(const std::list<Intersection_Line*> & discarded);

	void reschedule_events();

	void decrement_queued_events(const int n);


	Event_Vertex_Line* get_event_for_line(const int I_object) const;

	void delete_event(Event_Vertex_Line* e_vl);

	void delete_events(const std::list<Event_Vertex_Line*> & E_VL);


	std::pair<FT, bool> get_intersection_time(Intersection_Line* I) const;

	bool represents_same_intersection(Polygon_Vertex* v1, Polygon_Vertex* v2) const;

	static void set_paired_vertices(Polygon_Vertex* v_ts, Polygon_Vertex* v_os);

	static bool are_paired_vertices(Polygon_Vertex* v_ts, Polygon_Vertex* v_os);

	void set_paired_vertex(const Companion & _paired_vertex);

	void reset_paired_vertex();

	bool has_paired_vertex() const;

	bool is_independent() const;

	Polygon_Vertex* get_master_vertex() const;

protected:
	Event_Vertex_Line* make_event(Intersection_Line* I) const;

	Event_Vertex_Line* make_event(Intersection_Line* I, const FT & t) const;

	void copy_events(Polygon_Vertex* v_ts);

	void deschedule_events();

	bool represents_same_intersection(Polygon_Vertex* v) const;

public:
	CGAL_Point_2 M;
	CGAL_Vector_2 dM;
	FT t_init;
	FT t_stop;
	int K;

protected:
	std::map<int, Event_Vertex_Line*> lines_to_events;
	int queued_events;

protected:
	int k_th_interval;
	FT tau_inf;
	FT tau_sup;

public:
	Polygon_Edge* e1;
	Polygon_Edge* e2;

	bool is_active;
	int crossed_lines;

protected:
	std::list<Constraint> constraints;
	std::list<Polygon_Segment*> guided_segments;
	Polygon* polygon;

	Companion paired_vertex;
};


class Polygon_Edge : public Support_Plane_Object {
public:
	Polygon_Edge(const int _id_plane, Polygon_Vertex* _v1, Polygon_Vertex* _v2);

	Polygon_Edge(const int _id_plane, Polygon_Vertex* _v1, Polygon_Vertex* _v2, const Constraint & C);

	Polygon_Edge(const int _id_plane, Polygon_Vertex* _v1, Polygon_Vertex* _v2, Intersection_Line* I, Sign epsilon);

	~Polygon_Edge();

	bool is_constrained_by(Intersection_Line* I) const;

	// Different versions of intersection_pt_dir :
	// - we don't know if v1 or v2 already intersects I.
	void intersection_pt_dir(Intersection_Line* I, const FT & t, CGAL_Point_2 & M, CGAL_Vector_2 & dM) const;

	// - we know that v1 intersects I in V1_t.
	static void intersection_pt_dir(Intersection_Line* I, Polygon_Vertex* v1, Polygon_Vertex* v2, const FT & t, const CGAL_Point_2 & V1_t, CGAL_Point_2 & M, CGAL_Vector_2 & dM);

	Polygon_Vertex* intersection(Intersection_Line* I, Sign s, const FT & t, const int K, Event_Flags flags) const;

	Polygon_Vertex* intersection(Intersection_Line* I, Sign s, const FT & t, const CGAL_Point_2 & M, const CGAL_Vector_2 & dM, const int K, Event_Flags flags) const;

	bool is_adjacent_to_segment() const;

public:
	Polygon_Vertex* v1;
	Polygon_Vertex* v2;

	bool is_constrained;
	Constraint constraint;
};


class Segment : public Support_Plane_Object {
protected:
	Segment(const int _id_plane, Intersection_Line* _support);

public:
	virtual ~Segment();

	static bool closed_segment_includes(const CGAL_Point_2 & M, const CGAL_Point_2 & A, const CGAL_Point_2 & B);

	static bool half_closed_segment_includes(const CGAL_Point_2 & M, const CGAL_Point_2 & A, const CGAL_Point_2 & B);

public:
	Intersection_Line* support;	
};


class Planar_Segment : public Segment {
public:
	Planar_Segment(const int _id_plane, Intersection_Line* _support, CGAL_Point_2 & _A, CGAL_Point_2 & _B);

	~Planar_Segment();

	bool checks_if_belongs(const CGAL_Point_2 & V) const;

public:
	CGAL_Point_2 A;
	CGAL_Point_2 B;
};


class Polygon_Segment : public Segment {
protected:
	Polygon_Segment(const int _id_plane, const FT & _t_init, const Constraint & C_support);

public:
	Polygon_Segment(const int _id_plane, const FT & _t_init, const Constraint & C_support, Polygon_Vertex* & v, const CGAL_Point_2 & _O, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA);

	Polygon_Segment(const int _id_plane, const FT & _t_init, const Constraint & C_support, Polygon_Vertex* & v, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA);

	Polygon_Segment(const int _id_plane, const FT & _t_init, const Constraint & _C_init, const Constraint & C_support, Polygon_Vertex* & v, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA);

	Polygon_Segment(const int _id_plane, const FT & _t_init, const Constraint & _C_init, const Constraint & C_support, Polygon_Vertex* & v, const CGAL_Point_2 & _O, const CGAL_Point_2 & _A, const CGAL_Vector_2 & _dA);

	Polygon_Segment(const int _id_plane, const FT & t, const Constraint & _C_init, const Constraint & _C_support, const Constraint & _C_stop, const CGAL_Point_2 & _A, const CGAL_Point_2 & _B);

	~Polygon_Segment();


	bool stopped() const;

	FT t_min() const;

	CGAL_Point_2 origin() const;

	CGAL_Point_2 end() const;

	CGAL_Point_2 pt(FT t) const;



	void update_translation(const FT & t, const CGAL_Point_2 & A_t, const CGAL_Vector_2 & dA_t);
		
	void stop(const Constraint C, const FT & t);

	
	//bool checks_if_belongs(const CGAL_Point_2 & V_t, const FT & t);

	//bool checks_if_belongs(const CGAL_Point_2 & V_t, const Constraint C, const FT & t);

	//bool checks_if_belongs(const CGAL_Point_2 & V_t, const std::list<Constraint> & C_limits, const FT & t);



	static void set_as_opposite_bidirectional_segments(Polygon_Segment* s_1, Polygon_Segment* s_2);

	bool exists_opposite_segment_without_initinal_constraint() const;

	Polygon_Segment* get_opposite() const;

	void set_pseudo_init_constraint(const CGAL_Point_2 & A, std::list<Intersection_Line*> & L);

	void set_pseudo_init_constraint(const Constraint C_pseudo_init);

	Constraint get_pseudo_init_constraint() const;

	void insert_as_crossed_line(Intersection_Line* I);

	void insert_as_crossed_lines(const std::list<Intersection_Line*> & L);

	bool includes_point_on_support_line(const CGAL_Point_2 & M, const FT & t) const;

	bool includes_point_at_intersection(const CGAL_Point_2 & M, const Constraint & C, const FT & t) const;

	bool includes_point_at_intersection(const CGAL_Point_2 & M, const std::list<Constraint> & C_limits, const FT & t) const;

	bool includes_edge(const CGAL_Point_2 & V_1, const CGAL_Point_2 & V_2, const Constraint & C_1, const Constraint & C_2, const std::list<Intersection_Line*> & CL) const;

	std::list<Intersection_Line*>::const_iterator crossed_lines_begin() const;

	std::list<Intersection_Line*>::const_iterator crossed_lines_end() const;

	void check() const;

	void who_is() const;

public:
	Sign sign;

	FT t_init;
	FT t_stop;

	Segment_Translation* Tr;
	std::list<Segment_Translation*> Tr_previous;

protected:
	Constraint C_init;
	Constraint C_stop;
	
	std::list<Intersection_Line*> C_crossed;
	Polygon_Segment* opposite;
};

}