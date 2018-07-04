#include "partition_objects.h"
#include "vars.h"

namespace JPTD {
Partition_Edge::Partition_Edge(Partition_Vertex* _v1, Partition_Vertex* _v2)
	: id (++Counters::id_partition_edge)
{
	v1 = _v1;
	v2 = _v2;
	v1->push(this);
	v2->push(this);

	std::list<int> P1, P2;
	v1->get_definition(P1);
	v2->get_definition(P2);
	
	// We determine the ids of the planes to which this edge belongs,
	// by removing from one list of indices the elements that don't appear in the other list
	std::list<int>::iterator it_1 = P1.begin();
	while (it_1 != P1.end()) {
		std::list<int>::iterator it_2 = P2.begin();
		while (it_2 != P2.end()) {
			if ((*it_1) == (*it_2)) break;
			it_2++;
		}
		if (it_2 == P2.end()) {
			it_1 = P1.erase(it_1);
		} else {
			it_1++;
		}
	}

	// Assigns ids
	for (std::list<int>::iterator it = P1.begin() ; it != P1.end() ; it++) {
		local_ids[*it] = ++Counters::par_e_local_ids[*it];
	}
}


Partition_Edge::Partition_Edge(Partition_Vertex* _v1, Partition_Vertex* _v2, int F_i, int F_j)
	: id (++Counters::id_partition_edge)
{
	v1 = _v1;
	v2 = _v2;
	v1->push(this);
	v2->push(this);

	local_ids[F_i] = ++Counters::par_e_local_ids[F_i];
	local_ids[F_j] = ++Counters::par_e_local_ids[F_j];
}


Partition_Edge::~Partition_Edge()
{
	v1->pop(this);
	v2->pop(this);
}


void Partition_Edge::push(Partition_Facet* f) 
{ 
	facets.push_back(f); 
}



void Partition_Edge::pop(Partition_Facet* f)
{
	std::list<Partition_Facet*>::iterator it_f;
	for (it_f = facets.begin() ; it_f != facets.end() ; it_f++) {
		if ((*it_f) == f) break;
	}
	assert(it_f != facets.end());
	it_f = facets.erase(it_f);
}


int Partition_Edge::get_index_of_another_plane(const int excluded_id) const
{
	for (std::map<int, int>::const_iterator it = local_ids.cbegin() ; it != local_ids.cend() ; it++) {
		if (it->first != excluded_id) return it->first;
	}
		
	return -1;
}


int Partition_Edge::get_local_id(const int id) const
{
	std::map<int, int>::const_iterator it = local_ids.find(id);
	return (it != local_ids.end() ? it->second : -1);
}


int Partition_Edge::is_edge_of_another_bounding_facet(int F_i) const
{
	for (std::map<int, int>::const_iterator it = local_ids.cbegin(); it != local_ids.cend(); it++) {
		int F_j = it->first;
		if (F_j < 6) {
			if (F_j != F_i) return F_j;
		} else {
			return -1;
		}
	}
	return -1;
}


bool Partition_Edge::reaches(Partition_Vertex* v) const
{
	return (v1 == v || v2 == v);
}


Partition_Vertex* Partition_Edge::source(bool v1_v2) const
{ 
	return v1_v2 ? v1 : v2 ; 
}


Partition_Vertex* Partition_Edge::target(bool v1_v2) const
{ 
	return v1_v2 ? v2 : v1 ; 
}


Partition_Vertex* Partition_Edge::second_vertex(Partition_Vertex* v) const
{ 
	return (v1 == v ? v2 : v1); 
}


std::map<int, int>::const_iterator Partition_Edge::local_ids_begin() const 
{ 
	return local_ids.cbegin(); 
}


std::map<int, int>::const_iterator Partition_Edge::local_ids_end() const 
{ 
	return local_ids.cend(); 
}


std::map<int, int>::iterator Partition_Edge::local_ids_begin() 
{ 
	return local_ids.begin(); 
}


std::map<int, int>::iterator Partition_Edge::local_ids_end() 
{
	return local_ids.end(); 
}


std::list<Partition_Facet*>::const_iterator Partition_Edge::facets_begin() const
{
	return facets.cbegin();
}


std::list<Partition_Facet*>::const_iterator Partition_Edge::facets_end() const
{
	return facets.cend();
}


std::list<Partition_Facet*>::iterator Partition_Edge::facets_begin()
{
	return facets.begin();
}


std::list<Partition_Facet*>::iterator Partition_Edge::facets_end()
{
	return facets.end();
}
}