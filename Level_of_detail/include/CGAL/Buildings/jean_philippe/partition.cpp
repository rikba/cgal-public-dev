#include "partition.h"
#include "universe.h"
#include "support_plane.h"
#include "polygon_group.h"
#include "vars.h"
#include "ply_out.h"
#include <queue>
#include <CGAL/Direction_2.h>

namespace JPTD {

using CGAL::to_double;

Partition::Partition()
{
	// The six first planes listed in the map have the following equations :
	// (0) x = x_min, (1) x = x_max, (2) y = y_min, (3) y = y_max, (4) z = z_min, (5) z = z_max
	// We loop on these planes to get the definition of these six variables.
	dims = std::vector<double>(6, 0);
	for (int i = 0 ; i < 6 ; i++) {
		dims[i] = to_double(- Universe::map_of_planes[i]->plane.d());
	}

	// Once dimensions of the octree are known, we initialize it
	octree = new Partition_Vertex_Octree(dims);
	
	// We initialize a vector of lists of facets
	// It is going to collect facets and store them by plane index
	size_t n = Universe::map_of_planes.size();
	facets = std::vector<std::list<Partition_Facet*> >(n, std::list<Partition_Facet*>());

	// We also initialize the map that generates object ids for each plane
	Counters::par_v_local_ids = std::vector<int>(n, -1);
	Counters::par_e_local_ids = std::vector<int>(n, -1);
}


Partition::~Partition()
{
	// Deletes polyhedra
	for (std::list<Partition_Polyhedron*>::iterator it_p = polyhedrons.begin(); it_p != polyhedrons.end(); it_p++) delete (*it_p);
	polyhedrons.clear();

	// Deletes planar facets
	for (size_t i = 0; i < facets.size(); i++) {
		for (std::list<Partition_Facet*>::iterator it_f = facets[i].begin() ; it_f != facets[i].end() ; it_f++) {
			delete (*it_f);
		}
	}
	facets.clear();

	// Deletes edges
	for (std::list<Partition_Edge*>::iterator it_e = edges.begin(); it_e != edges.end(); it_e++) delete (*it_e);
	edges.clear();

	// Deletes vertices
	delete octree;
	octree = nullptr;

	dims.clear();
}


void Partition::build()
{
	// Part 1 : support planes of polygons
	// Groups of adjacent polygons inside trees directly define planar facets
	build_facets_for_support_planes();

	// Part 2 : facets of the bounding box
	// We first need to create the missing elements of the bounding box (corner vertices and edges)
	// After that we loop on the edges to get facets
	build_missing_vertices_for_bounding_box();
	build_missing_edges_for_bounding_box();
	build_facets_for_bounding_box();

	// Part 3 : cleans the partition
	debug();
	remove_bivalent_vertices();
	

	// Part 4 : get final polyhedrons
	build_polyhedrons();
}



void Partition::debug()
{
	std::vector<Partition_Vertex*> V;
	octree->get_all_vertices_sorted_by_identifier(V);

	FILE* file;
	
	file = fopen("/Users/danisimo/Documents/pipeline/logs/tmp/bad_input/txt/vertices.txt", "w");
	if (file != NULL) {
		for (int i = 0 ; i < V.size() ; i++) {
			Partition_Vertex* v = V[i];
			const CGAL_Point_3 & M = v->M;
			fprintf(file, "V[%i] : M = (%lf, %lf, %lf), IDS : ", i, to_double(M.x()), to_double(M.y()), to_double(M.z()));
			for (auto it = v->local_ids_begin(); it != v->local_ids_end(); it++) {
				fprintf(file, "(%i, %i) ", it->first, it->second);
			}
			fprintf(file, "\n");
		}
		fclose(file);
	}

	file = fopen("/Users/danisimo/Documents/pipeline/logs/tmp/bad_input/txt/edges.txt", "w");
	if (file != NULL) {
		for (std::list<Partition_Edge*>::iterator it_e = edges.begin() ; it_e != edges.end() ; it_e++) {
			Partition_Edge* e = (*it_e);
			fprintf(file, "E[%i] : v1 = %i, v2 = %i, IDS :", e->id, e->source(true)->id, e->target(true)->id);
			for (auto it = e->local_ids_begin() ; it != e->local_ids_end() ; it++) {
				fprintf(file, "(%i, %i) ", it->first, it->second);
			}
			fprintf(file, "\n");
		}
		fclose(file);
	}

	file = fopen("/Users/danisimo/Documents/pipeline/logs/tmp/bad_input/txt/facets.txt", "w");
	if (file != NULL) {
		for (int i = 0; i < facets.size(); i++) {
			const std::list<Partition_Facet*> & L = facets[i];

			fprintf(file, "PLANE %i\n", i);
			for (std::list<Partition_Facet*>::const_iterator it_f = L.begin() ; it_f != L.end() ; it_f++) {
				Partition_Facet* f = (*it_f);
				fprintf(file, "ID : %i, P : %i, edges : ", f->id, f->p);
				for (std::list<Partition_Edge*>::const_iterator it_e = f->edges_begin() ; it_e != f->edges_end() ; it_e++) {
					fprintf(file, "%i [%i %i] ", (*it_e)->id, (*it_e)->source(true)->id, (*it_e)->target(true)->id);
				}
				fprintf(file, "\n");
			}
		}
		fclose(file);
	}
}


void Partition::build_facets_for_support_planes()
{
	// Polygons are initially stored in trees, but we merge polygons when
	// they are adjacent and not separated by any Polygon_Segment.

	clock_t t_begin = clock();

	for (int p = 6 ; p < Universe::map_of_planes.size() ; p++) {
		Support_Plane* SP = Universe::map_of_planes[p];
		SP->group_final_polygons();
	}
	clock_t t_end = clock();
	//std::cout << "Group polygons : " << double(t_end - t_begin) / CLOCKS_PER_SEC << std::endl;

	// We loop on the different Support_Planes associated to input polygons.
	for (int p = 6 ; p < Universe::map_of_planes.size() ; p++) {
		Support_Plane* SP = Universe::map_of_planes[p];
		std::list<Partition_Facet*> faces_p;

		// We loop on all the obtained groups of convex polygons.
		// We get points that are on the convex hull of the group of polygons.
		for (std::list<Polygon_Group*>::iterator it_g = SP->groups.begin() ; it_g != SP->groups.end() ; it_g++) {
			Polygon_Group* G = (*it_g);
			
			// We loop on the set of Polygon_Edges identified as delimiters of the convex polygon group.
			// We take the two Polygon_Vertices at its ends, interpret them as intersections of 3 planes or more.
			// We insert these points in the octree if they do not exist yet, and build an edge.
			std::map<Polygon_Vertex*, Partition_Vertex*> conversions;
			std::map<int, std::pair<Partition_Edge*, Partition_Edge*> > adjacent_edges;

			for (std::list<Polygon_Edge*>::iterator it_e = G->borders_begin() ; it_e != G->borders_end() ; it_e++) {
				Polygon_Edge* pol_e = (*it_e);
				Polygon_Vertex *pol_v1 = pol_e->v1, *pol_v2 = pol_e->v2;

				// Matches a Partition_Vertex to a Polygon_Vertex, creates it if necessary
				Partition_Vertex* par_v1 = get_partition_vertex(pol_v1, conversions, adjacent_edges);
				Partition_Vertex* par_v2 = get_partition_vertex(pol_v2, conversions, adjacent_edges);

				// Matches a Partition_Edge to a Polygon_Edge
				get_partition_edge(par_v1, par_v2, adjacent_edges);
			}

			// Constructs one facet for G
			// Later, when we will assemble polyhedrons, this facet will be duplicated 
			// depending on the orientation of its normal
			std::list<Partition_Edge*> E;
			get_sequence_of_partition_edges(adjacent_edges, E);
			faces_p.push_back(new Partition_Facet(p, E));
		}

		std::copy(faces_p.begin(), faces_p.end(), std::back_inserter(facets[p]));
	}
}



Partition_Vertex* Partition::get_partition_vertex(Polygon_Vertex* pol_v, std::map<Polygon_Vertex*, Partition_Vertex*> & conversions,
	std::map<int, std::pair<Partition_Edge*, Partition_Edge*> > & adjacent_edges)
{
	// Before performing a more complex operation, we search if there already exists a partition vertex for pol_v
	// This make sense, since, as we process a Polygon_Group, we expect a Polygon_Vertex to appear twice in most cases
	// (Other cases correspond to intersections of 3 lines or more)
	
	std::map<Polygon_Vertex*, Partition_Vertex*>::iterator it_match = conversions.find(pol_v);
	if (it_match != conversions.end()) return it_match->second;

	// Searches for the Partition_Vertex that represents the Polygon_Vertex v,
	// using its definition as intersection of 3 planes or more.
	Support_Plane* SP = Universe::map_of_planes[pol_v->id_plane];
	
	// Gets indices of all planes intersecting in the point
	std::list<int> P;
	const Intersection_Line* I_1 = pol_v->get_constraint().first;
	const Intersection_Line* I_2 = pol_v->get_second_constraint().first;
	SP->get_indices_of_intersecting_planes(I_1, I_2, P);
	
	// Gets the 3D point
	// When we partition is built, vertices no longer move and therefore, its adjacent edges are constraint
	// To get a more precise 3D point, we compute it as the intersection of a plane and a line

	const CGAL_Point_3 M = SP->backproject(pol_v->M);

	// Checks if a vertex with the same definition already exists. 
	// If so, we return it
	// Two reasons why this case may happen : creation by another support plane, or 3 lines of the same plane intersect
	Partition_Vertex* par_v = octree->match(M, P);
	if (par_v == nullptr) {
		par_v = new Partition_Vertex(M, P);
		octree->add(par_v);
	}

	// As it is the first time that pol_v is processed in the current Polygon_Group, we keep trace of this result in the table of known conversions.
	// In the general case when only 2 lines intersect in a given point, the use of this table will spare the algorithm another search.
	conversions[pol_v] = par_v;
	
	// If par_v hadn't been identified as member of the Polygon_Group before, 
	// then we insert an empty entry for it in the table of adjacent edges
	if (adjacent_edges.find(par_v->id) == adjacent_edges.end()) adjacent_edges[par_v->id] = std::make_pair(nullptr, nullptr);

	return par_v;
}



Partition_Edge* Partition::get_partition_edge(Partition_Vertex* v1, Partition_Vertex* v2, 
	std::map<int, std::pair<Partition_Edge*, Partition_Edge*> > & adjacent_edges)
{
	Partition_Edge* e = nullptr;

	// Loops over the edges of v1, to get the edge that links it to v2
	for (std::list<Partition_Edge*>::iterator it_e = v1->edges_begin() ; it_e != v1->edges_end() ; it_e++) {
		if ((*it_e)->reaches(v2)) {
			e = (*it_e);
			break;
		}
	}

	// If we couldn't find an edge e = (v1 v2), then we create it
	if (e == nullptr) {
		e = new Partition_Edge(v1, v2);
		edges.push_back(e);
	}

	// Sets e as adjacent edge to v1 and v2
	if (adjacent_edges[v1->id].first == nullptr) {
		adjacent_edges[v1->id].first = e;
	} else {
		adjacent_edges[v1->id].second = e;
	}

	if (adjacent_edges[v2->id].first == nullptr) {
		adjacent_edges[v2->id].first = e;
	} else {
		adjacent_edges[v2->id].second = e;
	}

	return e;
}



void Partition::get_sequence_of_partition_edges(std::map<int, std::pair<Partition_Edge*, Partition_Edge*> > & adjacent_edges, 
	std::list<Partition_Edge*> & sorted_edges)
{
	for (auto it_e = adjacent_edges.begin() ; it_e != adjacent_edges.end() ; it_e++) {
		assert(it_e->second.first != nullptr && it_e->second.second != nullptr);
	}

	int v_init = adjacent_edges.begin()->first, v_curr = v_init;
	Partition_Edge* e_prev = nullptr;
	sorted_edges.clear();

	do {
		// Gets the first of the two edges associated to v_curr
		// Compares it to the last processed edge, to loop on a new vertex
		Partition_Edge* e = adjacent_edges[v_curr].first;
		if (e == e_prev) {
			e = adjacent_edges[v_curr].second;
		}
		sorted_edges.push_back(e);

		// Iterates
		Partition_Vertex *e_v1 = e->source(true), *e_v2 = e->target(true);
		v_curr = (e_v1->id == v_curr ? e_v2->id : e_v1->id);
		e_prev = e;

	} while (v_curr != v_init);
}



void Partition::build_missing_vertices_for_bounding_box()
{
	int ind_corners[8][3] = { {0, 2, 4}, {0, 3, 4}, {1, 3, 4}, {1, 2, 4}, {0, 2, 5}, {0, 3, 5}, {1, 3, 5}, {1, 2, 5} };
	for (int c = 0; c < 8; c++) {
		int i = ind_corners[c][0], j = ind_corners[c][1], k = ind_corners[c][2];

		// Constructs a 3D point at the intersection of 3 planes
		CGAL_Point_3 M(dims[i], dims[j], dims[k]);

		// Fills in a list representing the indices of the 3 intersecting planes
		std::list<int> P;
		P.push_back(i); P.push_back(j); P.push_back(k);

		// Formulates a query : if it results a non null result, then the vertex already exists
		Partition_Vertex* v = nullptr;
		if (v = octree->partial_match(M, P)) {
			continue;
		} else {
			v = new Partition_Vertex(M, P);
			octree->add(v);
		}
	}
}



void Partition::build_missing_edges_for_bounding_box()
{
	int ind_edges[12][2] = { {2, 4}, {3, 4}, {2, 5}, {3, 5}, {0, 4}, {1, 4}, {0, 5}, {1, 5}, {0, 2}, {1, 2}, {0, 3}, {1, 3} };
	for (int l = 0 ; l < 12 ; l++) {

		// Formulates box edges as intersections of two planes
		int F_i = ind_edges[l][0], F_j = ind_edges[l][1];

		// Gets vertices at the intersection of F_i and F_j
		std::list<Partition_Vertex*> l_V;
		octree->get_vertices_for_bounding_box_edge(F_i, F_j, l_V);

		// Sorts vertices
		std::vector<Partition_Vertex*> V;
		std::copy(l_V.begin(), l_V.end(), std::back_inserter(V));

		bool (*comparative_function) (Partition_Vertex*, Partition_Vertex*);
		if (l <= 3) {
			comparative_function = sort_by_x_coordinate;
		} else if (l <= 7) {
			comparative_function = sort_by_y_coordinate;
		} else {
			comparative_function = sort_by_z_coordinate;
		}
		std::sort(V.begin(), V.end(), comparative_function);

		// Builds vertices between each couple (V[i] V[i + 1]), upon condition that this edge doesn't exist yet
		for (int i = 0; i < V.size() - 1; i++) {
			Partition_Vertex *v1 = V[i], *v2 = V[i + 1];

			std::list<Partition_Edge*>::iterator it_e = v1->edges_begin();
			while (it_e != v1->edges_end()) {
				if ((*it_e)->reaches(v2)) break;
				it_e++;
			}

			if (it_e == v1->edges_end()) {
				edges.push_back(new Partition_Edge(v1, v2));
			}
		}
	}
}



void Partition::build_facets_for_bounding_box()
{
	// Let us loop on the different facets of the bounding box.
	// For each facet, we want to obtain the list of Partition_Vertices and Partition_Edges that belong to it.
	// This way we initialize a vector of directions and halfedges for each vertex.
	// Finally, we build a list of planar Partition_Facets, later assembled as Partition_Polyhedrons.

	for (int p = 0 ; p < 6 ; p++) {

		// Gets the list of vertices that belongs to facet p
		std::list<Partition_Vertex*> V;
		octree->get_vertices_for_bounding_box_facet(p, V);

		// Inits containers
		const int n = 1 + Counters::par_v_local_ids[p], m = 1 + Counters::par_e_local_ids[p];

		std::vector<std::vector<Direction_H> > D (n, std::vector<Direction_H>());
		bool **queued = new bool*[m];
		for (int i = 0; i < m; i++) {
			queued[i] = new bool[2];
			queued[i][0] = queued[i][1] = false;
		}

		// Loops on each vertex to get edges that belong to the same facet
		for (std::list<Partition_Vertex*>::iterator it_v = V.begin() ; it_v != V.end() ; it_v++) {
			Partition_Vertex* v = (*it_v);

			// For each edge e, gets its id in the plane p
			// If it differs from -1 then we take it into account
			for (std::list<Partition_Edge*>::iterator it_e = v->edges_begin() ; it_e != v->edges_end() ; it_e++) {
				Partition_Edge* e = (*it_e);
				int ind_e = e->get_local_id(p);
				if (ind_e == -1) continue;

				// If e->v1->id < e->v2->id, then two halfedges are created,
				// except when the edge is adjacent to two facets of the bounding box
				if (v->id < e->second_vertex(v)->id) {
					build_halfedges(p, e, D, queued);
				}
			}
		}

		// Loops, constructs facets
		loop_and_build_facets(p, D, queued);

		// Deletes containers
		D.clear();
		for (int i = 0; i < m; i++) delete[] queued[i];
		delete[] queued;
	}
}



void Partition::build_halfedges(const int F_i, Partition_Edge* e, std::vector<std::vector<Direction_H> > & directions, bool** & queued)
{
	Support_Plane* SP = Universe::map_of_planes[F_i];

	// We assume that e = (v_1 v_2)
	Partition_Vertex *v1 = e->source(true), *v2 = e->target(true);
	const int id_v1 = v1->get_local_id(F_i), id_v2 = v2->get_local_id(F_i);
	const int id_e = e->get_local_id(F_i);

	// We get the line equation of e
	int P = e->get_index_of_another_plane(F_i); assert(P != -1);
	Intersection_Line* I = SP->get_line(P);
	
	// Compute the angle made by I and determines if it is the angle of (v1 v2) or (v2 v1)
	const CGAL_Line_2 I_l = I->line;

	// n(I) = (a, b) => u(I) = (-b, a) => mes(u) = atan2(a, -b)
	//const double theta = atan2(to_double(-I_l.a()), to_double(I_l.b())); 
	//const CGAL_Vector_2 I_dir (I_l.b(), -I_l.a());
	const CGAL_Vector_2 v_12 = SP->project(v2->M) - SP->project(v1->M);
	const CGAL_Vector_2 v_21 = -v_12;

	/*double alpha_1, alpha_2;
	if (I_dir * v_12 > 0) {
		// (v1 v2) is in the same direction as (cos(theta), sin(theta))
		// So, v1 receives theta and v2 receives theta + PI or theta - PI
		alpha_1 = theta;
		alpha_2 = (theta > 0 ? theta - PI : theta + PI);
	} else {
		alpha_2 = theta;
		alpha_1 = (theta > 0 ? theta - PI : theta + PI);
	}*/

	// Constructs halfedge that are inserted at the correct entry in the table of directions
	add_halfedge(directions, id_v1, v_12, e, true);
	add_halfedge(directions, id_v2, v_21, e, false);

	const CGAL_Vector_3 v12 = v2->M - v1->M;
	const double dx = to_double(v12.x()), dy = to_double(v12.y()), dz = to_double(v12.z());
	int F_j = e->is_edge_of_another_bounding_facet(F_i);
	if (F_j != -1) {
		if (F_i <= 1) {
			// If the current bounding box facet is a plane X = M
			if (F_j == 5) {
				// Top facet : we only keep the halfedge s.t. dy > 0 
				if (dy > 0) queued[id_e][0] = true; else queued[id_e][1] = true;
			} else if (F_j == 4) {
				// Bottom facet : we only keep the halfedge s.t. dy < 0
				if (dy < 0) queued[id_e][0] = true; else queued[id_e][1] = true;
			} else if (F_j == 2) {
				// Left facet : we only keep the halfedge s.t. dz > 0
				if (dz > 0) queued[id_e][0] = true; else queued[id_e][1] = true;
			} else if (F_j == 3) {
				// Right facet : we only keep the halfedge s.t. dz < 0
				if (dz < 0) queued[id_e][0] = true; else queued[id_e][1] = true;
			}
		} 
		
		else if (F_i <= 3) {
			// If the current bounding box facet is a plane Y = M
			if (F_j == 5) {
				// Top facet : we only keep the halfedge s.t. dx < 0 
				if (dx < 0) queued[id_e][0] = true; else queued[id_e][1] = true;
			} else if (F_j == 4) {
				// Bottom facet : we only keep the halfedge s.t. dx > 0 
				if (dx > 0) queued[id_e][0] = true; else queued[id_e][1] = true;
			} else if (F_j == 1) {
				// Left facet : we only keep the halfedge s.t. dz > 0 
				if (dz > 0) queued[id_e][0] = true; else queued[id_e][1] = true;
			} else if (F_j == 0) {
				// Right facet : we only keep the halfedge s.t. dz < 0 
				if (dz < 0) queued[id_e][0] = true; else queued[id_e][1] = true;
			}
		} 
		
		else {
			// If the current bounding box facet is a plane Z = M
			if (F_j == 0) {
				// Top facet : we only keep the halfedge s.t. dy > 0 
				if (dy > 0) queued[id_e][0] = true; else queued[id_e][1] = true;
			} else if (F_j == 1) {
				// Bottom facet : we only keep the halfedge s.t. dy < 0 
				if (dy < 0) queued[id_e][0] = true; else queued[id_e][1] = true;
			} else if (F_j == 2) {
				// Left facet : we only keep the halfedge s.t. dx < 0
				if (dx < 0) queued[id_e][0] = true; else queued[id_e][1] = true;
			} else if (F_j == 3) {
				// Right facet : we only keep the halfedge s.t. dx > 0
				if (dx > 0) queued[id_e][0] = true; else queued[id_e][1] = true;
			}
		}
	}
}



/*void Partition::add_halfedge(std::vector<std::vector<Direction_H> > & D, int id_v, double theta, Partition_Edge* e, bool v1_v2)
{
	std::vector<Direction_H>::iterator it;
	for (it = D[id_v].begin() ; it != D[id_v].end() ; it++) {
		if (it->first > theta) break;
	}
	D[id_v].insert(it, std::make_pair(theta, std::make_pair(e, v1_v2)));
}*/



void Partition::add_halfedge(std::vector<std::vector<Direction_H> > & D, int id_v, const CGAL_Vector_2 & u_theta, Partition_Edge* e, bool v1_v2)
{
	CGAL_Direction_2 d_theta(u_theta);

	std::vector<Direction_H>::iterator it;
	for (it = D[id_v].begin() ; it != D[id_v].end() ; it++) {
		if (it->first > d_theta) break;
	}

	D[id_v].insert(it, std::make_pair(d_theta, std::make_pair(e, v1_v2)));
}



void Partition::loop_and_build_facets(const int F_i, std::vector<std::vector<Direction_H> > & directions, bool** & queued)
{
	std::queue<Partition_HalfEdge> queue;
	
	// In this algorithm, we are going to loop on halfedges to obtain planar facets using a FIFO queue.
	// Halfedges (v1 v2) have previously been added to directions[v1] with the angle made by (v1, v2) in the frame F_i.

	// When an edge e = (v1 v2) is inside the facet F_i, both halfedges are considered and can be queued.
	// However, when the edge is shared by another bounding box facet, only one halfedge can be used :
	// indeed, the use of the opposite halfedge would make us loop outside the facet.

	// For each bounding box facet F_i, the indices of top/bottom/left/right facets are known.
	// Like in the 2D version of the kinetic partition algorithm, edges on the top/bottom/left/right facets define a cycle. 
	// It is supposed to be a clockwise cycle, but due to the 3D orientation of the local frames of the bounding box facets,
	// the cycle may actually be counterclockwise.

	// When we process the halfedge (v1 v2), we normally iterate on the next halfedge (by increasing orientation)
	// starting from the vertex v2. But if the outer cycle is counterclockwise, then it is the previous one we should use.
	
	// To determine if the outer cycle is clockwise or counterclockwise, we get the "top left" corner C of each facet,
	// and determine if the halfedge (C ...) on the top border can be queued. If so the sequence is clockwise.


	// Step 1.
	// Initialization procedure

	// Gets the top left corner
	
	int neighbor_facets[3][4] = { {5, 4, 2, 3}, {5, 4, 1, 0}, {0, 1, 2, 3} };
	int F_t = neighbor_facets[F_i >> 1][0];
	int F_l = neighbor_facets[F_i >> 1][2];

	std::vector<int> P;
	P.push_back(F_i); P.push_back(F_t); P.push_back(F_l);
	std::sort(P.begin(), P.end());

	std::list<int> P_l;
	std::copy(P.begin(), P.end(), std::back_inserter(P_l));

	const CGAL_Point_3 M (dims[P[0]], dims[P[1]], dims[P[2]]);
	Partition_Vertex* C = octree->partial_match(M, P_l);

	// Gets the edge on the top border

	bool outer_cycle_is_clockwise;
	Partition_Edge* e = nullptr;
	int F_i_ind = -1, F_t_ind = -1;

	for (std::list<Partition_Edge*>::iterator it_e = C->edges_begin() ; it_e != C->edges_end() ; it_e++) {
		e = (*it_e);
		F_i_ind = e->get_local_id(F_i), F_t_ind = e->get_local_id(F_t);
		if (F_i_ind != -1 && F_t_ind != -1) {
			break;
		}
	}

	// Tests if the halfedge (C ...) is already marked as queued, in order to prevent the algorithm from using it

	if (e->source(true) == C) {
		if (!queued[F_i_ind][1]) {
			outer_cycle_is_clockwise = true;
			queue.push(std::make_pair(e, true));
		} else {
			outer_cycle_is_clockwise = false;
			queue.push(std::make_pair(e, false));
		}
	} else {
		if (!queued[F_i_ind][1]) {
			outer_cycle_is_clockwise = false;
			queue.push(std::make_pair(e, true));
		} else {
			outer_cycle_is_clockwise = true;
			queue.push(std::make_pair(e, false));
		}
	}

	// Step 2.
	// Main loop

	// As long as there remains halfedges to explore
	while (queue.size() > 0) {
	
		// We set a reference to the first element of the queue,
		// which represents the first iterated edge for obtaining a new facet
		Partition_HalfEdge h_e = queue.front();
		Partition_HalfEdge h_f = h_e;
		std::list<Partition_Edge*> F;

		// If this element hasn't been already processed simultaneously with other elements
		if (!queued[h_f.first->get_local_id(F_i)][h_f.second]) {
			queued[h_f.first->get_local_id(F_i)][h_f.second] = true;

			do {
				// Finds the next edge in clockwise order
				// Adds it to the facet that we are currently defining
				Partition_Vertex* v = h_f.first->target(h_f.second);
				int id_v = v->get_local_id(F_i);
				int path = 0;
				int n = int(directions[id_v].size());
				for (int i = 0 ; i < n ; i++) {
					if (directions[id_v][i].second.first == h_f.first) {
						path = (outer_cycle_is_clockwise ? (i + 1) % n : (i + n - 1) % n);
						break;
					}
				}

				h_f = directions[id_v][path].second;
				assert(v == h_f.first->source(h_f.second));
				F.push_back(h_f.first);

				// The current half-edge shouldn't be queued
				queued[h_f.first->get_local_id(F_i)][h_f.second] = true;

				// Insert the opposite half-edge into the queue, if it is not already done
				if (!queued[h_f.first->get_local_id(F_i)][!h_f.second]) {
					queue.push(std::make_pair(h_f.first, !h_f.second));
				}
			} while (h_f != h_e);

			// Inserts one facet
			facets[F_i].push_back(new Partition_Facet(F_i, F));
		}

		// Removes the first element of the queue, since we are done with it
		queue.pop();
	}
}



void Partition::remove_bivalent_vertices()
{
	// Gets all vertices of the partition
	std::list<Partition_Vertex*> V;
	octree->get_all_vertices_sorted_by_identifier(V);

	for (std::list<Partition_Vertex*>::iterator it_v = V.begin() ; it_v != V.end() ; it_v++) {
		assert((*it_v)->connectivity() >= 2);
	}

	// Maps every vertex and every edge to an iterator in the list that contains all of them
	std::map<int, std::list<Partition_Vertex*>::iterator> V_map;
	std::map<int, std::list<Partition_Edge*>::iterator> E_map;

	for (typename std::list<Partition_Vertex*>::iterator it_v = V.begin(); it_v != V.end(); ++it_v) {
		V_map[(*it_v)->id] = it_v;
		std::cout << (*it_v)->id << std::endl;
	}
	for (typename std::list<Partition_Edge*>::iterator it_e = edges.begin(); it_e != edges.end(); ++it_e) E_map[(*it_e)->id] = it_e;

	// Loops on every vertex
	std::list<Partition_Vertex*>::iterator it_v_main = V.begin();
	while (it_v_main != V.end()) {
		Partition_Vertex* v_main = (*it_v_main);
		if (v_main->connectivity() > 2) {
			++it_v_main;
		} 
		
		else {
			//std::cout << *v_main << std::endl;

			std::list<Partition_Edge*> edges_to_delete;
			std::list<Partition_Vertex*> vertices_to_delete;

			// If we find a vertex that is bound to only two edges, then we should remove it from the partition
			// That's why we loop on the vertices that delimit these edges. However, these vertices may be bivalent as well.
			// So, finally, we iterate on v1 and v2 until finding vertices that are not trivalent.
			// Until then, we save vertices and edges to delete

			Partition_Edge *e1 = v_main->edges_front(), *e2 = v_main->edges_back();
			Partition_Vertex *v1 = e1->second_vertex(v_main), *v2 = e2->second_vertex(v_main);

			// We define a sequence of vertices (v1 .. v_main .. v2)
			
			edges_to_delete.push_back(e1);
			while (v1->connectivity() == 2) {
				vertices_to_delete.push_back(v1);
				if (v1->edges_front() == e1) {
					e1 = v1->edges_back();
				} else {
					e1 = v1->edges_front();
				}
				edges_to_delete.push_back(e1);
				v1 = e1->second_vertex(v1);
			}

			edges_to_delete.push_back(e2);
			while (v2->connectivity() == 2) {
				vertices_to_delete.push_back(v2);
				if (v2->edges_front() == e2) {
					e2 = v2->edges_back();
				} else {
					e2 = v2->edges_front();
				}
				edges_to_delete.push_back(e2);
				v2 = e2->second_vertex(v2);
			}

			// Starts deleting references to the discarded vertices and edges
			// First of all we remove the edges to delete from the facets in which there are included

			std::list<Partition_Facet*> F;
			Partition_Edge* e = edges_to_delete.front();
			std::copy(e->facets_begin(), e->facets_end(), std::back_inserter(F));

			for (std::list<Partition_Facet*>::iterator it_f = F.begin() ; it_f != F.end() ; it_f++) {
				for (std::list<Partition_Edge*>::iterator it_e = edges_to_delete.begin() ; it_e != edges_to_delete.end() ; it_e++) {
					(*it_f)->pop(*it_e);
				}
			}

			for (std::list<Partition_Edge*>::iterator it_e = edges_to_delete.begin() ; it_e != edges_to_delete.end() ; it_e++) {
				Partition_Edge* e = (*it_e);
				edges.erase(E_map[e->id]);
				E_map.erase(e->id);
				delete e;
			}

			for (std::list<Partition_Vertex*>::iterator it_v = vertices_to_delete.begin() ; it_v != vertices_to_delete.end() ; it_v++) {
				Partition_Vertex* v = (*it_v);

				// std::cout << v->id << std::endl;

				// if (std::find(V_map.at(v->id) == V_map.end()) continue;

				V.erase(V_map.at(v->id));
				V_map.erase(v->id);
				octree->remove(v, true);
			}

			// We can now create a big edge (v1 v2), that we add to all the facets where edges have been removed

			Partition_Edge* e_12 = new Partition_Edge(v1, v2);
			edges.push_back(e_12);
			E_map[e_12->id] = --edges.end();

			for (std::list<Partition_Facet*>::iterator it_f = F.begin() ; it_f != F.end() ; it_f++) (*it_f)->push(e_12);

			// Finally iterates

			octree->remove(v_main, true);
			V_map.erase(v_main->id);
			it_v_main = V.erase(it_v_main);
		}
	}

	// Resets indices

	Counters::id_partition_vertex = -1;
	Counters::id_partition_edge = -1;
	for (std::list<Partition_Vertex*>::iterator it_v = V.begin(); it_v != V.end(); it_v++) (*it_v)->id = ++Counters::id_partition_vertex;
	for (std::list<Partition_Edge*>::iterator it_e = edges.begin(); it_e != edges.end(); it_e++) (*it_e)->id = ++Counters::id_partition_edge;
}



void Partition::build_polyhedrons()
{
	// Part 1.
	// For each edge, initializes circular vectors representing the orientations of all facets sides

	std::vector<std::vector<Direction_S> > directions;
	index_facets_sides(directions);

	/*for (int i = 0 ; i < Universe::map_of_planes.size() ; i++) {
		const CGAL_Plane & P_i = Universe::map_of_planes[i]->plane;
		std::cout << "P[" << i << "] : a = " << to_double(P_i.a()) << ", b = " << to_double(P_i.b()) << ", c = " << to_double(P_i.c()) << ", d = " << to_double(P_i.d()) << std::endl;
	}

	for (int i = 0 ; i < directions.size() ; i++) {
		std::vector<Direction_S> & directions_i = directions[i];
		std::cout << std::endl << "Edge " << i << std::endl;
		for (int j = 0 ; j < directions_i.size() ; j++) {
			Direction_S & d_ij = directions_i[j];
			std::cout << "theta = " << d_ij.first << ", side = (" << d_ij.second.first->id << ", " << (d_ij.second.second ? "TRUE)" : "FALSE)") << std::endl;
		}
	}*/

	// Part 2.
	// Initializes a table of sides to queue
	bool** queued;
	init_table_of_queued_sides(queued);

	// Part 3.
	// Main loop
	loop_and_build_polyhedrons(directions, queued);
	delete_table_of_queued_sides(queued);
}



void Partition::init_table_of_queued_sides(bool** & queued)
{
	// Initializes a matrix with as many rows as existing facets
	// Both sides of Partition_Facets can be queued, except those which are on the facets on the bounding box,
	// where we consider only the sides whose orthogonal vectors point towards the center of the box.
	
	const int n = 1 + Counters::id_partition_facet;

	queued = new bool*[n];
	for (int i = 0; i < n; i++) {
		queued[i] = new bool[2];
		queued[i][0] = queued[i][1] = false;
	}

	for (int p = 0 ; p < 6 ; p++) {
		for (std::list<Partition_Facet*>::iterator it_f = facets[p].begin() ; it_f != facets[p].end() ; it_f++) {
			Partition_Facet* f = (*it_f);
			
			// A Partition_Side (f, b) with b = true means that its normal is the normal of the support plane
			// Given the initial equations for the bounding box facets, we keep the positive (resp. negative) sides when p is odd (resp. even)

			if (p % 2 == 0) {
				// Excludes the negative side
				queued[f->id][0] = true;
			} else {
				// Excludes the positive side
				queued[f->id][1] = true;
			}
		}
	}
}



void Partition::delete_table_of_queued_sides(bool** & queued)
{
	const int n = 1 + Counters::id_partition_facet;

	for (int i = 0 ; i < n ; i++) {
		assert(queued[i][0] && queued[i][1]);
		delete[] queued[i];
	}
	
	delete[] queued;
	queued = nullptr;
}



void Partition::index_facets_sides(std::vector<std::vector<Direction_S> > & D)
{
	// Initializes D
	D = std::vector<std::vector<Direction_S> >(edges.size(), std::vector<Direction_S>());

	// Loops on all edges and apply the indexation procedure
	for (std::list<Partition_Edge*>::iterator it_e = edges.begin() ; it_e != edges.end() ; it_e++) {
		Partition_Edge* e = (*it_e);
		std::vector<Direction_S> & D_edge = D[e->id];

		// Gets a reference planar facet
		std::list<Partition_Facet*>::iterator it_f = e->facets_begin();
		Partition_Facet* F_ref = (*it_f);
		
		// We define a local frame such that O is one of the two vertices of e,
		// v_i is included in P_ref, v_j is orthogonal to P_ref, and v_k = v_i ^ v_j follows the direction of e.
		Partition_Vertex* O = nullptr;
		CGAL_Vector_3 v_i, v_j, v_k;
		get_local_frame(e, F_ref, O, v_i, v_j, v_k);

		// We are now ready to assign angles to the different planar facets containing e.
		// The idea is to compute the angle of the half-line obtained by intersection of any facet F and P_ref = (O->M, v_k).
		// To this end, we find a vertex v of F that doesn't define e, project it on E and measure the angle (v_i, vec(O, v_proj).
		// Finally, we insert a couple of values in D.

#if TEST_DIRECTION_S
		CGAL_Direction_2 i(1, 0);
		D_edge.push_back(std::make_pair(i, std::make_pair(F_ref, false)));
		D_edge.push_back(std::make_pair(i, std::make_pair(F_ref, true)));

		while (++it_f != e->facets_end()) {
			Partition_Facet* F = (*it_f);

			if (F->p == F_ref->p) {
				// Particular case : when the facet is on the same plane as F_ref, 
				// we already know that we should insert it at PI.
				add_side(D_edge, -i, F, true);
			
			} else {
				Support_Plane* SP = Universe::map_of_planes[F->p];

				// Finds a vertex to project on E
				Partition_Vertex* v = F->get_projectible_vertex(e);
				assert(v != nullptr);

				// Gets A, the coordinates of v, and obtains B by translating A along the normal vector of SP.
				const CGAL_Point_3 & A = v->M;
				const CGAL_Point_3 B = v->M + CGAL_Vector_3(SP->plane.a(), SP->plane.b(), SP->plane.c());

				// Projects A to get an indexaction value for F in the circulator
				CGAL_Direction_2 u_theta_a = get_angle_of_projected_vertex(A, O, v_i, v_j, v_k);

				// Now, in which order should we insert the sides (F, true) and (F, false) ?
				// The answer is given by another projection (v + SP(F)->orthonormal_vector)
				CGAL_Direction_2 u_theta_b = get_angle_of_projected_vertex(B, O, v_i, v_j, v_k);

				bool front_side_is_positive = u_theta_b.counterclockwise_in_between(u_theta_a, -u_theta_a);

				// Finally indexes F
				add_side(D_edge, u_theta_a, F, front_side_is_positive);
			}
		}
#else
		D_edge.push_back(std::make_pair(0, std::make_pair(F_ref, true)));

		while (++it_f != e->facets_end()) {
			Partition_Facet* F = (*it_f);

			if (F->p == F_ref->p) {
				// Particular case : when the facet is on the same plane as F_ref, 
				// we already know that we should insert it at PI.
				add_side(D_edge, PI, F, true);

			} else {
				Support_Plane* SP = Universe::map_of_planes[F->p];

				// Finds a vertex to project on E
				Partition_Vertex* v = F->get_projectible_vertex(e);
				assert(v != nullptr);

				// Gets A, the coordinates of v, and obtains B by translating A along the normal vector of SP.
				const CGAL_Point_3 & A = v->M;
				const CGAL_Point_3 B = v->M + CGAL_Vector_3(SP->plane.a(), SP->plane.b(), SP->plane.c());

				// Projects A to get an indexaction value for F in the circulator
				double theta = get_angle_of_projected_vertex(A, O, v_i, v_j, v_k);

				// Now, in which order should we insert the sides (F, true) and (F, false) ?
				// The answer is given by another projection (v + SP(F)->orthonormal_vector)
				double theta_b = get_angle_of_projected_vertex(B, O, v_i, v_j, v_k);

				// We set a boolean for telling the algorithm if, we will later loop on the circulator,
				// the positive side of F is the one that should appear first
				bool front_side_is_positive = false;
				if (theta < PI) {
					front_side_is_positive = (theta_b < theta || theta_b > theta + PI);
				} else {
					front_side_is_positive = (theta - PI < theta_b && theta_b < theta);
				}

				// Finally indexes F
				add_side(D_edge, theta, F, front_side_is_positive);
			}
		}

		D_edge.push_back(std::make_pair(2 * PI, std::make_pair(F_ref, false)));
#endif
	}
}


#if TEST_DIRECTION_S
void Partition::add_side(std::vector<Direction_S> & D, const CGAL_Direction_2 & u_theta, Partition_Facet* F, bool positive_side_comes_first)
{
	std::vector<Direction_S>::iterator it;
	for (it = D.begin() ; it != D.end() ; it++) {
		if (it->first > u_theta) break;
	}

	if (positive_side_comes_first) {
		// If the positive side of F must appear first,
		// then the negative side of F is the first element to insert
		it = D.insert(it, std::make_pair(u_theta, std::make_pair(F, false)));
		it = D.insert(it, std::make_pair(u_theta, std::make_pair(F, true)));
	} 
	
	else {
		it = D.insert(it, std::make_pair(u_theta, std::make_pair(F, true)));
		it = D.insert(it, std::make_pair(u_theta, std::make_pair(F, false)));
	}
}
#else
void Partition::add_side(std::vector<Direction_S> & D, double theta, Partition_Facet* F, bool positive_side_comes_first)
{
	std::vector<Direction_S>::iterator it;
	for (it = D.begin() ; it != D.end() ; it++) {
		if (it->first > theta) break;
	}

	if (positive_side_comes_first) {
		// If the positive side of F must appear first,
		// then the negative side of F is the first element to insert
		it = D.insert(it, std::make_pair(theta, std::make_pair(F, false)));
		it = D.insert(it, std::make_pair(theta, std::make_pair(F, true)));
	} 
	
	else {
		it = D.insert(it, std::make_pair(theta, std::make_pair(F, true)));
		it = D.insert(it, std::make_pair(theta, std::make_pair(F, false)));
	}
}
#endif


void Partition::get_local_frame(Partition_Edge* e, Partition_Facet* f, Partition_Vertex* & O, CGAL_Vector_3 & v_i, CGAL_Vector_3 & v_j, CGAL_Vector_3 & v_k)
{
	Support_Plane* SP = Universe::map_of_planes[f->p];
	const CGAL_Plane & H_SP = SP->plane;
	const CGAL_Vector_3 N = H_SP.orthogonal_vector();

	// We consider that f is the reference side for computing the orientations of all sides adjacent to e.
	// In the circulator that we will later obtain, the side (f, true) makes an angle of 0 with the local frame, 
	// and the side (f, false) makes an angle of 2 * PI.

	Partition_Vertex *v_1 = e->source(true), *v_2 = e->target(true);
	CGAL_Line_3 e_line (v_1->M, v_2->M);
	CGAL_Vector_3 e_line_dir = e_line.to_vector();
	CGAL_Vector_3 e_ortho_dir = CGAL::cross_product(N, e_line_dir);

	// Searches for a vertex v_3 not aligned with v_1 and v_2
	Partition_Vertex* v_3 = f->get_projectible_vertex(e);
	
	CGAL_Vector_3 v_12 = v_2->M - v_1->M, v_13 = v_3->M - v_1->M;
	CGAL_Vector_3 v_n = CGAL::cross_product(v_12, v_13);

	// Tests if (v_12, v_13, v_12 ^ v_13) is correctly oriented, as we assume to be on the positive side of f
	// A correct orientation means that v_12 ^ v_13 = K * N, with K > 0
	// (i.e. v_2, v_1, v_3 make a right turn on the side of f that has a positive orientation)

	FT K = (CGAL::abs(N.x()) > 0 ? v_n.x() / N.x() : (CGAL::abs(N.y()) > 0 ? v_n.y() / N.y() : v_n.z() / N.z()));

	if (K > 0) {
		// If so we center the frame in v_1
		O = v_1;
		
		CGAL_Line_3 e_ortho_line(v_1->M, e_ortho_dir);
		CGAL_Point_3 v_3_proj = e_ortho_line.projection(v_3->M);
		v_k = CGAL::cross_product(v_3_proj - v_1->M, N);
	} else {
		// Otherwise we center the frame in v_2
		O = v_2;

		CGAL_Line_3 e_ortho_line(v_2->M, e_ortho_dir);
		CGAL_Point_3 v_3_proj = e_ortho_line.projection(v_3->M);
		v_k = CGAL::cross_product(v_3_proj - v_2->M, N);
	}

	v_j = N;
	v_i = CGAL::cross_product(v_j, v_k);

	// In the end we get a local and direct frame (but NOT normalized) centered in one of the two edges of e,
	// s.t. v_j is the normal vector to the positive side of f, v_k is collinear to e, and v_i is included in f
}


#if TEST_DIRECTION_S
CGAL_Direction_2 Partition::get_angle_of_projected_vertex(const CGAL_Point_3 & M, Partition_Vertex* v_O, const CGAL_Vector_3 & v_i, const CGAL_Vector_3 & v_j, const CGAL_Vector_3 & v_k)
{
	// First we find t s.t. v->M + t * v_k belongs to (O, v_k)
	const CGAL_Point_3 & O = v_O->M;
	const FT a = v_k.x(), b = v_k.y(), c = v_k.z(), d = -a * O.x() - b * O.y() - c * O.z();
	const FT t = -(a * M.x() + b * M.y() + c * M.z() + d) / (a * a + b * b + c * c);
	const CGAL_Point_3 M_t = M + t * v_k;

	// There exists a couple (l, m) s.t. V_t = O + l * v_i + m * v_j
	// We obtain a system with 3 equations and two unknown variables
	FT l, m;
	const FT det_xy = (v_i.x() * v_j.y() - v_i.y() * v_j.x());
	const FT det_xz = (v_i.x() * v_j.z() - v_i.z() * v_j.x());
	const FT det_yz = (v_i.y() * v_j.z() - v_i.z() * v_j.y());
	assert((CGAL::abs(det_xy) > 0) || (CGAL::abs(det_xz) > 0) || (CGAL::abs(det_yz) > 0));

	if (CGAL::abs(det_xy) > 1e-6) {
		const FT dx = (M_t.x() - O.x()), dy = (M_t.y() - O.y());
		l = ((dx * v_j.y() - dy * v_j.x()) / det_xy);
		m = ((v_i.x() * dy - v_i.y() * dx) / det_xy);
	} else if (CGAL::abs(det_xz) > 1e-6) {
		const FT dx = (M_t.x() - O.x()), dz = (M_t.z() - O.z());
		l = ((dx * v_j.z() - dz * v_j.x()) / det_xz);
		m = ((v_i.x() * dz - v_i.z() * dx) / det_xz);
	} else {
		const FT dy = (M_t.y() - O.y()), dz = (M_t.z() - O.z());
		l = ((dy * v_j.z() - dz * v_j.y()) / det_yz);
		m = ((v_i.y() * dz - v_i.z() * dy) / det_yz);
	}

	// (l, m) are coordinates of V_t in the 2D frame (O, v_i, v_j))
	// We finally compute an angle in [0, 2 * PI] representing mes(v_i, vec(O, V_t))
	return CGAL_Direction_2(l, m);
}
#else
double Partition::get_angle_of_projected_vertex(const CGAL_Point_3 & M, Partition_Vertex* v_O, const CGAL_Vector_3 & v_i, const CGAL_Vector_3 & v_j, const CGAL_Vector_3 & v_k)
{
	// First we find t s.t. v->M + t * v_k belongs to (O, v_k)
	const CGAL_Point_3 & O = v_O->M;
	const FT a = v_k.x(), b = v_k.y(), c = v_k.z(), d = -a * O.x() - b * O.y() - c * O.z();
	const FT t = -(a * M.x() + b * M.y() + c * M.z() + d) / (a * a + b * b + c * c);
	const CGAL_Point_3 M_t = M + t * v_k;

	// There exists a couple (l, m) s.t. V_t = O + l * v_i + m * v_j
	// We obtain a system with 3 equations and two unknown variables
	double l, m;
	const FT det_xy = (v_i.x() * v_j.y() - v_i.y() * v_j.x());
	const FT det_xz = (v_i.x() * v_j.z() - v_i.z() * v_j.x());
	const FT det_yz = (v_i.y() * v_j.z() - v_i.z() * v_j.y());
	assert((CGAL::abs(det_xy) > 0) || (CGAL::abs(det_xz) > 0) || (CGAL::abs(det_yz) > 0));

	if (CGAL::abs(det_xy) > 1e-6) {
		const FT dx = (M_t.x() - O.x()), dy = (M_t.y() - O.y());
		l = to_double((dx * v_j.y() - dy * v_j.x()) / det_xy);
		m = to_double((v_i.x() * dy - v_i.y() * dx) / det_xy);
	} else if (CGAL::abs(det_xz) > 1e-6) {
		const FT dx = (M_t.x() - O.x()), dz = (M_t.z() - O.z());
		l = to_double((dx * v_j.z() - dz * v_j.x()) / det_xz);
		m = to_double((v_i.x() * dz - v_i.z() * dx) / det_xz);
	} else {
		const FT dy = (M_t.y() - O.y()), dz = (M_t.z() - O.z());
		l = to_double((dy * v_j.z() - dz * v_j.y()) / det_yz);
		m = to_double((v_i.y() * dz - v_i.z() * dy) / det_yz);
	}

	// (l, m) are coordinates of V_t in the 2D frame (O, v_i, v_j))
	// We finally compute an angle in [0, 2 * PI] representing mes(v_i, vec(O, V_t))
	double theta = atan2(m, l);
	if (theta < 0) theta += 2 * PI;
	return theta;
}
#endif


void Partition::loop_and_build_polyhedrons(std::vector<std::vector<Direction_S> > & directions, bool** & queued)
{
	std::queue<Partition_Side> outer_queue;

	// Initialization of the outer queue :
	// Takes the first Partition_Facet found on the first Support_Plane (defined by x = x_min)

	Partition_Facet* f_init = *facets[0].begin();
	outer_queue.push(std::make_pair(f_init, true));

	// As long as the remains Partition_Sides to explore

	while (!outer_queue.empty()) {	
		Partition_Side & S_0 = outer_queue.front();
	
		if (!queued[S_0.first->id][S_0.second]) {

			// We have obtained the first side of a new polyhedron
			// We initialize an inner queue with it : it represents the list of facets to loop on
			std::set<Partition_Side> P;
			std::queue<Partition_Side> inner_queue;
			inner_queue.push(S_0);

			while (!inner_queue.empty()) {
				Partition_Side & S = inner_queue.front();
				
				if (!queued[S.first->id][S.second]) {
					queued[S.first->id][S.second] = true;
					
					Partition_Facet* f = S.first;
					P.insert(S);

					// Inserts the opposite side into the outer queue
					if (!queued[S.first->id][!S.second]) {
						outer_queue.push(std::make_pair(S.first, !S.second));
					}

					// Loops on the adjacent facets of F
					// Inserts them into the inner queue if they have not already been reached before
					for (std::list<Partition_Edge*>::iterator it_e = f->edges_begin() ; it_e != f->edges_end() ; it_e++) {
						const int id_e = (*it_e)->id, n = int(directions[id_e].size());
						
						int j = 0;
						for (int i = 0 ; i < n ; i++) {
							Partition_Side & S_i = directions[id_e][i].second;
							if (S_i.first == S.first && S_i.second == S.second) {
								j = i;
								break;
							}
						}

						// Accesses to the two sides just before and after i in the circulator
						// One is the opposite side of S, one is adjacent to S inside the polyhedron P
						int j_prev = (j + n - 1) % n, j_next = (j + 1) % n;
						Partition_Side & S_prev = directions[id_e][j_prev].second, & S_next = directions[id_e][j_next].second;

						// Gets the adjacent side and pushes it, if not already done
						Partition_Side & S_adj = (S_prev.first == S.first ? S_next : S_prev);
						if (!queued[S_adj.first->id][S_adj.second]) {
							inner_queue.push(S_adj);
						}
					}
				}

				inner_queue.pop();
			}

			// Creates a new polyhedron with the 
			polyhedrons.push_back(new Partition_Polyhedron(P));
		}

		outer_queue.pop();
	}
}



void Partition::ply_facets(const std::string & filename) const
{
	std::list<CGAL_Point_3> P;
	std::list<std::list<int> > F;
	std::list<CGAL_Color> C;

	// Gets all 3D points

	std::vector<Partition_Vertex*> V;
	octree->get_all_vertices_sorted_by_identifier(V);
	assert(V.size() > 0);
	for (size_t i = 0; i < V.size(); i++) P.push_back(V[i]->M);

	// Loops on the different support planes
	for (size_t i = 0 ; i < facets.size() ; i++) {

		// Given a plane, loops on the facets
		for (std::list<Partition_Facet*>::const_iterator it_f = facets[i].begin() ; it_f !=  facets[i].end() ; it_f++) {
			Partition_Facet* f = (*it_f);
			if (f->p < 6) continue;

			// Given a facet, gets a sequence of vertices
			std::list<Partition_Vertex*> f_ver;
			f->get_circular_sequence_of_vertices(f_ver, true);

			// Converts vertices to indices
			std::list<int> f_ind;
			for (std::list<Partition_Vertex*>::iterator it_v = f_ver.begin(); it_v != f_ver.end(); it_v++) f_ind.push_back((*it_v)->id);

			F.push_back(f_ind);

			// Prints a color for this facet
			C.push_back(Universe::map_of_planes[i]->color);
		}
	}

	Ply_Out::print_plain_colorful_facets(filename, P, F, C, 25);
}



void Partition::ply_individual_polyhedron(const std::string & filename, const int K) const
{
	// Containers for PLY file
	std::list<CGAL_Point_3> P;
	std::list<std::list<int> > F;

	// Accesses the polyhedron K
	std::list<Partition_Polyhedron*>::const_iterator it_p = polyhedrons.begin();
	for (int i = 0; i < K; i++) ++it_p;
	Partition_Polyhedron* PH = (*it_p);

	// Lists vertices for all sides of the polyhedron
	// Reverts the orientation for a nicer display
	std::list<std::list<Partition_Vertex*> > sequences_per_side;
	std::set<Partition_Vertex*> vertices_used;

	for (std::list<Partition_Side>::iterator it_s = PH->facets_begin() ; it_s != PH->facets_end() ; it_s++) {
		Partition_Facet* f = it_s->first;

		std::list<Partition_Vertex*> f_ver;
		f->get_circular_sequence_of_vertices(f_ver, !it_s->second);

		for (std::list<Partition_Vertex*>::iterator it_v = f_ver.begin() ; it_v != f_ver.end() ; it_v++) {
			vertices_used.insert(*it_v);
		}

		sequences_per_side.push_back(f_ver);
	}

	// This polyhedron doesn't use all the vertices of the octree, and we don't need to print them all,
	// but in the PLY file, vertices start from 0. That's why we set a conversion map.
	//std::vector<Partition_Vertex*> V;
	//octree->get_all_vertices_sorted_by_identifier(V);

	std::map<Partition_Vertex*, int> conversions;
	for (std::set<Partition_Vertex*>::iterator it_v = vertices_used.begin() ; it_v != vertices_used.end() ; it_v++) {
		Partition_Vertex* v = (*it_v);
		P.push_back(v->M);
		conversions[v] = int(P.size()) - 1;
	}

	for (std::list<std::list<Partition_Vertex*> >::iterator it_lv = sequences_per_side.begin() ; it_lv != sequences_per_side.end() ; it_lv++) {
		std::list<Partition_Vertex*> & seq = (*it_lv);
		std::list<int> seq_indices;
		for (std::list<Partition_Vertex*>::iterator it_v = seq.begin() ; it_v != seq.end() ; it_v++) {
			seq_indices.push_back(conversions[*it_v]);
		}
		F.push_back(seq_indices);
	}
	
	Ply_Out::print_colorless_facets(filename, P, F);
}



void Partition::get_all_vertices(std::list<Partition_Vertex*> & V) const
{
	octree->get_all_vertices(V);
}


void Partition::get_all_vertices_sorted_by_identifier(std::vector<Partition_Vertex*> & V) const
{
	octree->get_all_vertices_sorted_by_identifier(V);
}


std::list<Partition_Facet*>::const_iterator Partition::planar_facets_begin(const int id) const
{
	return facets[id].cbegin();
}


std::list<Partition_Facet*>::const_iterator Partition::planar_facets_end(const int id) const
{
	return facets[id].cend();
}


std::list<Partition_Facet*>::iterator Partition::planar_facets_begin(const int id)
{
	return facets[id].begin();
}


std::list<Partition_Facet*>::iterator Partition::planar_facets_end(const int id)
{
	return facets[id].end();
}


size_t Partition::polyhedrons_size() const
{
	return polyhedrons.size();
}


std::list<Partition_Polyhedron*>::const_iterator Partition::polyhedrons_begin() const 
{ 
	return polyhedrons.cbegin(); 
}


std::list<Partition_Polyhedron*>::iterator Partition::polyhedrons_begin()
{
	return polyhedrons.begin();
}


std::list<Partition_Polyhedron*>::const_iterator Partition::polyhedrons_end() const 
{ 
	return polyhedrons.cend(); 
}


std::list<Partition_Polyhedron*>::iterator Partition::polyhedrons_end()
{
	return polyhedrons.end();
}

}