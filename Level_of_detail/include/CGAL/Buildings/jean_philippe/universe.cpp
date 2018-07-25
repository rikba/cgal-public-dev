#include "universe.h"
#include "parameters.h"
#include <chrono>

namespace JPTD {

#pragma warning(suppress:4244)
std::default_random_engine Universe::generator (std::chrono::system_clock::now().time_since_epoch().count());

Kinetic_Parameters* Universe::params = new Kinetic_Parameters();

std::map<int, Polygon_Vertex*> Universe::map_of_objects;

std::vector<Support_Plane*> Universe::map_of_planes;

int Universe::moving_objects = 0;

Event_Queue* Universe::event_queue = nullptr;

Octree_Base<CGAL_Point_3, CGAL_Inexact_Point_3>* Universe::point_cloud_octree = nullptr;

}