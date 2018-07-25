#pragma once

#include "defs.h"
#include <map>
#include <list>
#include <vector>
#include <random>

#include "defs_cgal.h"
#include "octree_base.h"

namespace JPTD {

class Kinetic_Parameters;

class Polygon_Vertex;

class Support_Plane;

class Event_Queue;

class Partition;


namespace Universe
{
	extern std::default_random_engine generator;

	extern Kinetic_Parameters* params;

	extern std::map<int, Polygon_Vertex*> map_of_objects;

	extern std::vector<Support_Plane*> map_of_planes;

	extern int moving_objects;

	extern Event_Queue* event_queue;

	extern Octree_Base<CGAL_Point_3, CGAL_Inexact_Point_3>* point_cloud_octree;
};

}