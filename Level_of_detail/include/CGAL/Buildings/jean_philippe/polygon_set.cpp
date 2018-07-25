#include "polygon_set.h"
#include "support_plane.h"
#include "universe.h"

namespace JPTD {

Polygon_Set::Polygon_Set(const std::map<int, Intersection_Line*> & L)
{
	// Builds a map of correspondances between lines and bits of a signature
	int entry = -1;
	for (std::map<int, Intersection_Line*>::const_iterator it_l = L.begin() ; it_l != L.end() ; it_l++) {
		dictionary[it_l->second] = ++entry;
	}
}


Polygon_Set::~Polygon_Set()
{
	dictionary.clear();

	for (std::map<Signature, Polygon_Cell*, Vector_Bool_Comparator>::iterator it_s = cells.begin() ; it_s != cells.end() ; it_s++) {
		delete it_s->second;
	}
}


void Polygon_Set::insert(const Signature & S, Polygon* P)
{
	// Searches for an entry with the same signature
	// If it does exist, we only add P to the cell
	// (This branch will only be used for initialization of coplanar polygons)
	// Otherwise we create a new cell containing the polygon

	std::map<Signature, Polygon_Cell*, Vector_Bool_Comparator>::iterator it_s = cells.find(S);

	if (it_s != cells.end()) {
		it_s->second->push(P);
	} else {
		cells[S] = new Polygon_Cell(S, P);
	}
}


bool Polygon_Set::exists(const Signature & S, const int seed)
{
	std::map<Signature, Polygon_Cell*, Vector_Bool_Comparator>::iterator it_s = cells.find(S);
	if (it_s == cells.end()) return false;

	Polygon_Cell* C = it_s->second;
	for (std::list<Polygon*>::iterator it_p = C->polygons_begin() ; it_p != C->polygons_end() ; it_p++) {
		Polygon* P = (*it_p);
		if (P->seed == seed) return true;
	}

	return false;

	// return (cells.find(S) != cells.end());
}


Polygon* Polygon_Set::get_adjacent_polygon(Polygon* P_ts, Intersection_Line* I)
{
	// Constructs the signature that corresponds to the polygon adjacent to the one containing v,
	// but located on the other side of the intersection line I, reached by v	
	std::vector<bool> S_os = get_adjacent_polygons_signature(P_ts, I);

	// Searches for such a polygon, returns it
	std::map<Signature, Polygon_Cell*, Vector_Bool_Comparator>::iterator it_s = cells.find(S_os);
	if (it_s == cells.end()) {
		return nullptr;
	} else {
		for (std::list<Polygon*>::iterator it_p = it_s->second->polygons_begin() ; it_p != it_s->second->polygons_end() ; it_p++) {
			Polygon* P = (*it_p);
			if (P->seed == P_ts->seed) {
				return P;
			}
		}
		return nullptr;
		//return it_s->second->get_unique_polygon();
	}
}


std::vector<bool> Polygon_Set::get_adjacent_polygons_signature(Polygon* P, Intersection_Line* I)
{
	// Gets original signature of P
	std::vector<bool> & S_ts = P->get_cell()->get_signature();

	// Sets the bit that corresponds to line I
	std::vector<bool> S_os;
	std::copy(S_ts.begin(), S_ts.end(), std::back_inserter(S_os));

	int b = dictionary[I];
	S_os[b] = !S_os[b];

	// Returns result
	return S_os;
}


std::vector<bool> Polygon_Set::get_adjacent_polygons_signature(Polygon* P, Intersection_Line* I_1, Intersection_Line* I_2)
{
	// Same as before, except that we set two bits
	std::vector<bool> & S_ts = P->get_cell()->get_signature();

	std::vector<bool> S_os;
	std::copy(S_ts.begin(), S_ts.end(), std::back_inserter(S_os));

	int b_1 = dictionary[I_1];
	int b_2 = dictionary[I_2];
	S_os[b_1] = !S_os[b_1];
	S_os[b_2] = !S_os[b_2];

	return S_os;
}


void Polygon_Set::get_signature_of_adjacent_cell(std::vector<bool> & S, Intersection_Line* I)
{
	int b = dictionary[I];
	S[b] = !S[b];
}



void Polygon_Set::get_polygons(std::list<Polygon*> & P)
{
	P.clear();

	// There may be more than one polygon per cell
	// However, assuming that this function is called once all points have stopped propagating,
	// we only take one of those polygons into account (they all have the same definition in the end)

	for (std::map<Signature, Polygon_Cell*, Vector_Bool_Comparator>::iterator it_s = cells.begin() ; it_s != cells.end() ; it_s++) {
		Polygon_Cell* C = it_s->second;
		P.push_back(C->get_one_polygon());
	}
}


void Polygon_Set::get_polygon_description(std::list<std::list<CGAL_Point_3> > & polygons, std::list<CGAL_Color> & colors, const double t)
{
	// Gets the description of all polygons at time t.
	// Their assigned colors are those of their support planes.

	for (std::map<Signature, Polygon_Cell*, Vector_Bool_Comparator>::iterator it_c = cells.begin(); it_c != cells.end(); it_c++) {
		Polygon_Cell* C = it_c->second;
		for (std::list<Polygon*>::iterator it_p = C->polygons_begin(); it_p != C->polygons_end(); it_p++) {

			Polygon* polygon = (*it_p);
			std::list<Polygon_Vertex*> & vertices = polygon->vertices;
			Polygon_Vertex *v_init = *(vertices.begin()), *v_curr = v_init, *v_next = nullptr;
			Polygon_Edge* e_prev = nullptr;

			std::vector<CGAL_Point_3> V;
			size_t n = vertices.size();
			V.reserve(n);

			// Gets the plane equation
			Support_Plane* SP = Universe::map_of_planes[v_init->id_plane];
			const CGAL_Plane & H = SP->plane;
			const CGAL_Vector_3 & h = H.orthogonal_vector();

			// Gets a list of backprojected vertices
			while (true) {
				CGAL_Point_2 M = v_curr->pt(t);
				V.push_back(SP->backproject(M));
				Polygon_Edge* e = (v_curr->e1 == e_prev ? v_curr->e2 : v_curr->e1);
				Polygon_Vertex* v_next = (e->v1 == v_curr ? e->v2 : e->v1);
				if (v_next == v_init) {
					break;
				} else {
					e_prev = e;
					v_curr = v_next;
				}
			}

			// Builds facets
			std::list<CGAL_Point_3> P(1, V[0]);
			bool forward = (CGAL::cross_product(V[1] - V[0], V[2] - V[0]) * h > 0);
			if (forward) {
				for (size_t i = 1; i < n; i++) P.push_back(V[i]);
			} else {
				for (size_t i = n - 1; i > 0; i--) P.push_back(V[i]);
			}
			polygons.push_back(P);
			colors.push_back(SP->color);
		}
	}
}

}