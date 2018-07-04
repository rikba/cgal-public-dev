#include "support_plane_object.h"
#include "vars.h"
#include <iostream>

namespace JPTD {

Support_Plane_Object::Support_Plane_Object(const int _id_plane)
	: id_object (++Counters::id_objects),
	id_plane (_id_plane)
{

}


Support_Plane_Object::~Support_Plane_Object()
{

}

}