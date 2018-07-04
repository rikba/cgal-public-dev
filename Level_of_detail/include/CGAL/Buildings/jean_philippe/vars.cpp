#include "vars.h"

namespace JPTD {
int Counters::id_planes = -1;
int Counters::id_objects = -1;

int Counters::id_partition_vertex = -1;
int Counters::id_partition_edge = -1;
int Counters::id_partition_facet = -1;
int Counters::id_partition_polyhedron = -1;

std::vector<int> Counters::par_v_local_ids = std::vector<int>();
std::vector<int> Counters::par_e_local_ids = std::vector<int>();
}