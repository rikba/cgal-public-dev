#pragma once
#include "support_plane_objects.h"
#include "partition_objects.h"
#include "partition_vertex_octree.h"

namespace JPTD {

#define TEST_DIRECTION_S 0

class Partition
{
public:
	Partition();

	~Partition();

	void build();
	
protected:
	typedef std::pair<CGAL_Direction_2, Partition_HalfEdge> Direction_H;

#if TEST_DIRECTION_S
	typedef std::pair<CGAL_Direction_2, Partition_Side> Direction_S;
#else
	typedef std::pair<double, Partition_Side> Direction_S;
#endif

	void build_facets_for_support_planes();

	Partition_Vertex* get_partition_vertex(Polygon_Vertex* v, std::map<Polygon_Vertex*, Partition_Vertex*> & conversions,
		std::map<int, std::pair<Partition_Edge*, Partition_Edge*> > & adjacent_edges);

	Partition_Edge* get_partition_edge(Partition_Vertex* v1, Partition_Vertex* v2, 
		std::map<int, std::pair<Partition_Edge*, Partition_Edge*> > & adjacent_edges);

	void get_sequence_of_partition_edges(std::map<int, std::pair<Partition_Edge*, Partition_Edge*> > & adjacent_edges, std::list<Partition_Edge*> & sorted_edges);


	void build_missing_vertices_for_bounding_box();

	void build_missing_edges_for_bounding_box();

	void build_facets_for_bounding_box();

	void build_halfedges(const int F_i, Partition_Edge* e, std::vector<std::vector<Direction_H> > & directions, bool** & queued);

	void add_halfedge(std::vector<std::vector<Direction_H> > & D, int id_v, const CGAL_Vector_2 & u_theta, Partition_Edge* e, bool v1_v2);

	//void add_halfedge(std::vector<std::vector<Direction_H> > & D, int id_v, double theta, Partition_Edge* e, bool v1_v2);
	

	void loop_and_build_facets(const int F_i, std::vector<std::vector<Direction_H> > & directions, bool** & queued);


	void remove_bivalent_vertices();


	void build_polyhedrons();

	void index_facets_sides(std::vector<std::vector<Direction_S> > & D);

	void get_local_frame(Partition_Edge* e, Partition_Facet* f, Partition_Vertex* & O, CGAL_Vector_3 & v_i, CGAL_Vector_3 & v_j, CGAL_Vector_3 & v_k);

#if TEST_DIRECTION_S
	void add_side(std::vector<Direction_S> & D, const CGAL_Direction_2 & u_theta, Partition_Facet* F, bool positive_side_comes_first);

	CGAL_Direction_2 get_angle_of_projected_vertex(const CGAL_Point_3 & A, Partition_Vertex* O, const CGAL_Vector_3 & v_i, const CGAL_Vector_3 & v_j, const CGAL_Vector_3 & v_k);
#else
	void add_side(std::vector<Direction_S> & D, double theta, Partition_Facet* f, bool positive_side_comes_first);

	double get_angle_of_projected_vertex(const CGAL_Point_3 & A, Partition_Vertex* O, const CGAL_Vector_3 & v_i, const CGAL_Vector_3 & v_j, const CGAL_Vector_3 & v_k);
#endif

	void init_table_of_queued_sides(bool** & queued);

	void delete_table_of_queued_sides(bool** & queued);

	void loop_and_build_polyhedrons(std::vector<std::vector<Direction_S> > & directions, bool** & queued);

	void debug();

public:
	KINETIC_PARTITION_API void ply_facets(const std::string & filename) const;

	KINETIC_PARTITION_API void ply_individual_polyhedron(const std::string & filename, const int P) const;


	KINETIC_PARTITION_API void get_all_vertices(std::list<Partition_Vertex*> & V) const;

	KINETIC_PARTITION_API void get_all_vertices_sorted_by_identifier(std::vector<Partition_Vertex*> & V) const;


	KINETIC_PARTITION_API std::list<Partition_Facet*>::const_iterator planar_facets_begin(const int id) const;

	KINETIC_PARTITION_API std::list<Partition_Facet*>::const_iterator planar_facets_end(const int id) const;

	KINETIC_PARTITION_API std::list<Partition_Facet*>::iterator planar_facets_begin(const int id);

	KINETIC_PARTITION_API std::list<Partition_Facet*>::iterator planar_facets_end(const int id);

	
	KINETIC_PARTITION_API std::list<Partition_Polyhedron*>::iterator polyhedrons_begin();

	KINETIC_PARTITION_API std::list<Partition_Polyhedron*>::const_iterator polyhedrons_begin() const;

	KINETIC_PARTITION_API std::list<Partition_Polyhedron*>::iterator polyhedrons_end();

	KINETIC_PARTITION_API std::list<Partition_Polyhedron*>::const_iterator polyhedrons_end() const;

	KINETIC_PARTITION_API size_t polyhedrons_size() const;

	Partition_Vertex_Octree* octree;
	std::vector<std::list<Partition_Facet*> > facets;

private:
	std::vector<double> dims;

	std::list<Partition_Edge*> edges;
	std::list<Partition_Polyhedron*> polyhedrons;
};

}