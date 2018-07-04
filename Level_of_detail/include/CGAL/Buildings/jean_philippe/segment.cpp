#include "support_plane_objects.h"

namespace JPTD {
Segment::Segment(const int _id_plane, Intersection_Line* _support)
	: Support_Plane_Object(_id_plane)
{
	support = _support;
	support->segments.push_back(this);
}


Segment::~Segment()
{

}
}