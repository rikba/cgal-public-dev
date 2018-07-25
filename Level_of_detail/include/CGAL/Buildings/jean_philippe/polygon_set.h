#pragma once
#include <set>
#include "polygon_cell.h"

namespace JPTD {

class Polygon_Set
{
public:
	Polygon_Set(const std::map<int, Intersection_Line*> & L);

	~Polygon_Set();

	void insert(const Signature & S, Polygon* P);

	bool exists(const Signature & S, const int seed);

	Polygon* get_adjacent_polygon(Polygon* P_ts, Intersection_Line* I);

	std::vector<bool> get_adjacent_polygons_signature(Polygon* P, Intersection_Line* I);

	std::vector<bool> get_adjacent_polygons_signature(Polygon* P, Intersection_Line* I_1, Intersection_Line* I_2);

	void get_signature_of_adjacent_cell(std::vector<bool> & S, Intersection_Line* I);

	void get_polygons(std::list<Polygon*> & P);

	inline std::map<Signature, Polygon_Cell*, Vector_Bool_Comparator>::iterator cells_begin() { return cells.begin(); }
	inline std::map<Signature, Polygon_Cell*, Vector_Bool_Comparator>::iterator cells_end() { return cells.end(); }

	void get_polygon_description(std::list<std::list<CGAL_Point_3> > & polygons, std::list<CGAL_Color> & colors, const double t);

	std::map<Intersection_Line*, int> dictionary;
protected:
	std::map<Signature, Polygon_Cell*, Vector_Bool_Comparator> cells;
};

}