#include "polygon_cell.h"

namespace JPTD {

Polygon_Cell::Polygon_Cell(const Signature & S, Polygon* P)
{
	signature = S;
	push(P);
}


Polygon_Cell::~Polygon_Cell()
{
	for (std::list<Polygon*>::iterator it_p = polygons.begin() ; it_p != polygons.end() ; it_p++) {
		delete (*it_p);
	}

	polygons.clear();
}


Signature Polygon_Cell::make_signature(const CGAL_Point_2 & O, const std::map<int, Intersection_Line*> & L)
{
	// Initializes a vector
	Signature S(L.size(), false);

	int entry = -1;
	for (std::map<int, Intersection_Line*>::const_iterator it_l = L.begin() ; it_l != L.end() ; it_l++) {
		const CGAL_Line_2 & l = it_l->second->line;
		S[++entry] = (l.a() * O.x() + l.b() * O.y() + l.c() > 0);
	}

	return S;
}


Polygon* Polygon_Cell::get_one_polygon() const
{
	return polygons.front();
}


void Polygon_Cell::push(Polygon* P)
{	
	polygons.push_back(P); 
	P->set_cell(this); 
}


size_t Polygon_Cell::size() const
{ 
	return polygons.size(); 
}


const std::list<Polygon*>::iterator Polygon_Cell::polygons_begin()
{ 
	return polygons.begin(); 
}


const std::list<Polygon*>::iterator Polygon_Cell::polygons_end()
{
	return polygons.end();
}

}