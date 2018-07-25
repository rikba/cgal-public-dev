#include "partition_objects.h"
#include "vars.h"

namespace JPTD {

Partition_Polyhedron::Partition_Polyhedron(const std::set<Partition_Side> & P)
	: id (++Counters::id_partition_polyhedron)
{
	for (std::set<Partition_Side>::iterator it_f = P.begin() ; it_f != P.end() ; it_f++) {
		facets.push_back(*it_f);
		it_f->first->push(this);
	}
}


Partition_Polyhedron::~Partition_Polyhedron()
{
	for (std::list<Partition_Side>::iterator it_f = facets.begin() ; it_f != facets.end() ; it_f++) {
		it_f->first->pop(this);
	}
}


std::list<Partition_Side>::const_iterator Partition_Polyhedron::facets_begin() const
{
	return facets.cbegin();
}


std::list<Partition_Side>::const_iterator Partition_Polyhedron::facets_end() const
{
	return facets.cend();
}


std::list<Partition_Side>::iterator Partition_Polyhedron::facets_begin()
{
	return facets.begin();
}


std::list<Partition_Side>::iterator Partition_Polyhedron::facets_end()
{
	return facets.end();
}

}