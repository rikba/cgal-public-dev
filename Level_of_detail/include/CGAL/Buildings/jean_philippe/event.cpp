#include "event.h"
#include <iterator>
#include <iostream>

namespace JPTD {

Event_Vertex_Line::Event_Vertex_Line(const int _intersectant, const int _intersected, const FT & _t_intersectant)
	: intersectant (_intersectant),
	intersected (_intersected),
	t_intersectant (_t_intersectant),
	is_queued (false)
{
	// The iterator is left uninitialized.
	// Indeed, the event is not assigned yet to an entry of the queue.
}


Event_Vertex_Line::~Event_Vertex_Line()
{
}

}