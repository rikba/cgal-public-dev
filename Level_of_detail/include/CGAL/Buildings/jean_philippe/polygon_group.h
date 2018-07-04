#pragma once
#include "polygon.h"

namespace JPTD {
class Polygon_Group
{
public:
	Polygon_Group(Polygon* P);

	~Polygon_Group();

	//void append(Polygon* P);

	void append(Polygon* P, std::list<Polygon_Edge*> & E_remove, std::list<Polygon_Edge*> & E_not_insert);

	//bool is_adjacent_to(Polygon* P);

	bool is_adjacent_to(Polygon* P, std::list<Polygon_Edge*> & E_remove, std::list<Polygon_Edge*> & E_not_insert);

	//void merge(Polygon_Group* G);

	void merge(Polygon_Group* G, std::list<Polygon_Edge*> E_remove, std::list<Polygon_Edge*> & E_not_insert);

	void get_convex_hull(std::vector<CGAL_Point_2> & H);

	void get_convex_hull(std::list<std::tuple<int, int, CGAL_Point_2> > & H);

	void get_extended_convex_hull(std::set<Polygon_Vertex*> & V);

	void get_extended_convex_hull(std::list<std::tuple<int, int, CGAL_Point_2> > & H);

	inline std::list<Polygon_Edge*>::iterator borders_begin() { return borders.begin(); }
	inline std::list<Polygon_Edge*>::iterator borders_end() { return borders.end(); }

private:
	void remove_edges_from_borders(std::list<Polygon_Edge*> & E);

	void add_selected_edges_to_borders(Polygon* P, std::list<Polygon_Edge*> & E);

	void add_selected_edges_to_borders(Polygon_Group* G, std::list<Polygon_Edge*> & E);

private:
	std::list<Polygon*> polygons;
	std::list<std::tuple<int, int, int, CGAL_Point_2> > intersections;
	std::list<Polygon_Edge*> borders;
};
}