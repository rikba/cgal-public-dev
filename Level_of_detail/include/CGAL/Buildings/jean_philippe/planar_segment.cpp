#include "support_plane_objects.h"
#include "support_plane.h"
#include "universe.h"

namespace JPTD {
Planar_Segment::Planar_Segment(const int _id_plane, Intersection_Line* _support, CGAL_Point_2 & _A, CGAL_Point_2 & _B)
	: Segment (_id_plane, _support)
{
	Universe::map_of_planes[_id_plane]->borders[id_object] = this;

	A = _A;
	B = _B;
}


Planar_Segment::~Planar_Segment()
{
	//std::map<int, Planar_Segment*>::iterator it_s = Universe::map_of_planes[id_plane]->borders.find(id_object);
	//Universe::map_of_planes[id_plane]->borders.erase(it_s);
}


bool Planar_Segment::checks_if_belongs(const CGAL_Point_2 & V) const
{
	// Returns true if V is on [AB], false otherwise

	CGAL_Segment_2 AB (A, B);
	return AB.has_on(V);
}
}