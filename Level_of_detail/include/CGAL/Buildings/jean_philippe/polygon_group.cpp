#include "polygon_group.h"
#include "support_plane_objects.h"
#include <CGAL/convex_hull_2.h>

namespace JPTD {

Polygon_Group::Polygon_Group(Polygon* P)
{
	// We assume that this function is called when all vertices have stopped propagating
	// Therefore vertices have exactly two constraints

	polygons.push_back(P);
	
	for (std::list<Polygon_Vertex*>::iterator it_v = P->vertices.begin() ; it_v != P->vertices.end() ; it_v++) {
		Polygon_Vertex* v = (*it_v);
		
		const Constraint C_1 = v->get_constraint();
		const Constraint C_2 = v->get_second_constraint();
		const CGAL_Point_2 pt = v->M;

		int l_min, l_max;
		if (C_1.first->id_object < C_2.first->id_object) {
			l_min = C_1.first->id_object;
			l_max = C_2.first->id_object;
		} else {
			l_min = C_2.first->id_object;
			l_max = C_1.first->id_object;
		}

		intersections.push_back(std::make_tuple(l_min, l_max, 1, pt));
	}

	// We initialize a map that tells us if a certain edge is part of the convex hull of the group
	// For the moment, all edges are part of the convex hull

	for (std::list<Polygon_Edge*>::iterator it_e = P->edges.begin(); it_e != P->edges.end() ; it_e++) {
		borders.push_back(*it_e);
	}
}


Polygon_Group::~Polygon_Group()
{
}

/*
void Polygon_Group::append(Polygon* P)
{
	// P is added to the list of polygons belonging to this group
	polygons.push_back(P);

	// We loop on all vertices of P, which are defined as the crossing of two Intersection_Lines (l_1, l_2)
	// We search for such an entry in the list of known intersections, if we find such an element the number of occurences is incremented
	// Otherwise we add a new entry to this list
	for (std::map<int, Polygon_Vertex*>::iterator it_v = P->vertices.begin() ; it_v != P->vertices.end() ; it_v++) {
		Polygon_Vertex* v = it_v->second;
		
		const Constraint C_1 = v->get_constraint();
		const Constraint C_2 = v->get_second_constraint();
		const CGAL_Point_2 pt = v->M;

		int l_min, l_max;
		if (C_1.first->id_object < C_2.first->id_object) {
			l_min = C_1.first->id_object;
			l_max = C_2.first->id_object;
		} else {
			l_min = C_2.first->id_object;
			l_max = C_1.first->id_object;
		}

		bool construct_tuple = true;
		for (std::tuple<int, int, int, CGAL_Point_2> & T : intersections) {
			if (std::get<0>(T) == l_min && std::get<1>(T) == l_max) {
				std::get<2>(T) += 1;
				construct_tuple = false;
				break;
			}
		}

		if (construct_tuple) {
			intersections.push_back(std::make_tuple(l_min, l_max, 1, pt));
		}
	}
}
*/


void Polygon_Group::remove_edges_from_borders(std::list<Polygon_Edge*> & E)
{
	// We loop on the list of edges currently listed as borders
	// If one of them is actually in E, then this edge is removed from the list of borders
	// We also remove it from E, since it is no longer useful to compare other edges to that value

	std::list<Polygon_Edge*>::iterator it_e1 = borders.begin(), it_e2;

	while (it_e1 != borders.end() && !E.empty()) {
		bool exclude = false;
		for (it_e2 = E.begin() ; it_e2 != E.end() ; it_e2++) {
			if ((*it_e1) == (*it_e2)) {
				exclude = true;
				E.erase(it_e2);
				break;
			}
		}
		if (exclude) {
			it_e1 = borders.erase(it_e1);
		} else {
			it_e1++;
		}
	}
}


void Polygon_Group::add_selected_edges_to_borders(Polygon* P, std::list<Polygon_Edge*> & E)
{
	// We insert all edges of P to the list of borders, except those which are listed in E
	for (std::list<Polygon_Edge*>::iterator it_e1 = P->edges.begin(); it_e1 != P->edges.end() ; it_e1++) {
		bool include = true;

		std::list<Polygon_Edge*>::iterator it_e2 = E.begin();
		while (it_e2 != E.end()) {
			if ((*it_e1) == (*it_e2)) {
				include = false;
				E.erase(it_e2);
				break;
			}
			it_e2++;
		}

		if (include) borders.push_back(*it_e1);
	}
}


void Polygon_Group::add_selected_edges_to_borders(Polygon_Group* G, std::list<Polygon_Edge*> & E)
{
	// We insert all elements of G to the list of borders, except those which are listed in E
	for (std::list<Polygon_Edge*>::iterator it_e1 = G->borders.begin(); it_e1 != G->borders.end() ; it_e1++) {
		bool include = true;

		std::list<Polygon_Edge*>::iterator it_e2 = E.begin();
		while (it_e2 != E.end()) {
			if ((*it_e1) == (*it_e2)) {
				include = false;
				E.erase(it_e2);
				break;
			}
			it_e2++;
		}

		if (include) borders.push_back(*it_e1);
	}
}


void Polygon_Group::append(Polygon* P, std::list<Polygon_Edge*> & edges_to_remove, std::list<Polygon_Edge*> & edges_not_to_insert)
{
	// P is added to the list of polygons belonging to this group
	polygons.push_back(P);

	// Updates the borders of the current region
	remove_edges_from_borders(edges_to_remove);
	add_selected_edges_to_borders(P, edges_not_to_insert);
}


/*
bool Polygon_Group::is_adjacent_to(Polygon* P)
{
	// Convenience function for testing the adjacency of P to any of the polygons of this group
	
	for (std::list<Polygon*>::iterator it_p = polygons.begin() ; it_p != polygons.end() ; it_p++) {
		if ((*it_p)->is_adjacent_to(P)) {
			return true;
		}
	}

	return false;
}
*/


bool Polygon_Group::is_adjacent_to(Polygon* P, std::list<Polygon_Edge*> & edges_to_remove, std::list<Polygon_Edge*> & E_not_insert)
{
	// Convenience function for testing the adjacency of P to any of the polygons of this group
	edges_to_remove.clear();
	E_not_insert.clear();
	
	for (std::list<Polygon*>::iterator it_p = polygons.begin() ; it_p != polygons.end() ; it_p++) {
		Polygon_Edge *e, *f;

		// Every time we find an edge of P that is adjacent to any polygon of this group, we update E and P_E
		if ((*it_p)->is_adjacent_to(P, e, f)) {
			edges_to_remove.push_back(e);
			E_not_insert.push_back(f);
		}
	}

	return (!edges_to_remove.empty());
}



void Polygon_Group::merge(Polygon_Group* G, std::list<Polygon_Edge*> edges_to_remove, std::list<Polygon_Edge*> & edges_not_to_insert)
{
	// Adds all the polygons listed in G to this list of polygons
	std::move(G->polygons.begin(), G->polygons.end(), std::back_inserter(polygons));

	// Remove a subset of edges from the list of borders,
	// and adds those of G, except those listed in the aforementioned container
	remove_edges_from_borders(edges_to_remove);
	add_selected_edges_to_borders(G, edges_not_to_insert);

	// Deletes G
	delete G;
}


void Polygon_Group::get_convex_hull(std::vector<CGAL_Point_2> & H)
{
	// Gets the list of the intersection points between lines 
	// that are not shared by two polygons or more inside the group

	std::vector<CGAL_Point_2> V;
	V.reserve(intersections.size());

	for (std::tuple<int, int, int, CGAL_Point_2> & T : intersections) {
		if (std::get<2>(T) == 1) {
			V.push_back(std::get<3>(T));
		}
	}

	H.clear();
	CGAL::convex_hull_2(V.begin(), V.end(), std::back_inserter(H));
}


void Polygon_Group::get_convex_hull(std::list<std::tuple<int, int, CGAL_Point_2> > & H)
{
	// Gets the list of intersection points between two Intersection_Lines,
	// that are not shared by two polygons or more inside the group.

	// Warning : if 3 lines (L_1, L_2 and L_3) are concurrent then the same intersection point
	// may be duplicated, as we insert (L_1, L_2), (L_1, L_3), etc.

	H.clear();
	for (std::tuple<int, int, int, CGAL_Point_2> & T : intersections) {
		if (std::get<2>(T) == 1) {
			H.push_back(std::make_tuple(std::get<0>(T), std::get<1>(T), std::get<3>(T)));
		}
	}
}


void Polygon_Group::get_extended_convex_hull(std::set<Polygon_Vertex *> & V)
{
	// We loop on the list of edges marked as borders and get their vertices
	// We do not care if the vertices don't really represent the convex hull,
	// and are in fact located between vertices of the exact convex hull

	V.clear();
	for (std::list<Polygon_Edge*>::iterator it_e = borders.begin() ; it_e != borders.end() ; it_e++) {
		V.insert((*it_e)->v1);
		V.insert((*it_e)->v2);
	}
}


void Polygon_Group::get_extended_convex_hull(std::list<std::tuple<int, int, CGAL_Point_2> > & H)
{
	// Gets a set of vertices located on the convex hull
	// The returned set may actually contain more entries than the exact convex hull

	std::set<Polygon_Vertex*> V;
	get_extended_convex_hull(V);

	// Now, loops on the list of vertices we have obtained
	// and inserts points as intersections of two Intersection_Lines

	H.clear();
	for (std::set<Polygon_Vertex*>::iterator it_v = V.begin() ; it_v != V.end() ; it_v++) {
		int i_1 = (*it_v)->get_constraint().first->id_object;
		int i_2 = (*it_v)->get_second_constraint().first->id_object;
		CGAL_Point_2 M = (*it_v)->M;

		H.push_back(std::make_tuple(i_1, i_2, M));
	}
}

}