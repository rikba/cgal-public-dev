#pragma once
#include <map>
#include <vector>
#include "event.h"

namespace JPTD {

class Event_Queue
{
public:
	Event_Queue();

	~Event_Queue();

public:
	bool has_events() const;

	Event_Vertex_Line* pop();

	void push(Event_Vertex_Line* e_vl);

	void push(const std::list<Event_Vertex_Line*> & E_VL);

	void push(Event_Vertex_Line* e_vl, const std::map<FT, std::list<Event_Vertex_Line*> >::iterator it_push);

	void erase(Event_Vertex_Line* e_vl);

	void erase(const std::list<Event_Vertex_Line*> & E_VL);

	void get_simultaneous_events_for_this_vertex(const int intersectant, const FT & t_intersectant, std::list<Event_Vertex_Line*> & E);

protected:
	void print() const;

protected:
	std::map<FT, std::list<Event_Vertex_Line*> > queue;
};

}
