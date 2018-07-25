#include "support_plane.h"
#include "polygon_group.h"
#include "universe.h"
#include "event_queue.h"
#include "parameters.h"
#include "vars.h"
#include "stats.h"
#include <iostream>
#include <CGAL/convex_hull_2.h>

namespace JPTD {

using CGAL::to_double;


Support_Plane::Support_Plane(const CGAL_Plane & _plane)
	: id (++Counters::id_planes),
	// rtree_boxes_width (0.5),
	// rtree_boxes_per_line (Universe::params->boxes),
	plane (_plane)
{
	std::uniform_int_distribution<int> color_picker (100, 255);
	uchar r = color_picker(Universe::generator);
	uchar g = color_picker(Universe::generator);
	uchar b = color_picker(Universe::generator);
	color = CGAL_Color(r, g, b);

	Universe::map_of_planes.push_back(this);

	polygon_set = nullptr;

	x_min = FLT_MAX, x_max = -FLT_MAX;
	y_min = x_min, y_max = x_max;
}


Support_Plane::~Support_Plane()
{
	// rtree_boxes_to_lines.clear();
	// rtree_lines.clear();
	concurrences.clear();

	// Destroys directions
	for (std::vector<Polygon_Directions*>::iterator it_d = polygon_directions.begin(); it_d != polygon_directions.end(); it_d++) {
		delete (*it_d);
	}

	// Destroys groups of polygons
	for (std::list<Polygon_Group*>::iterator it_g = groups.begin(); it_g != groups.end(); it_g++) {
		delete (*it_g);
	}
	groups.clear();

	// Destroys polygons, vertices and edges
	delete polygon_set;
	vertices.clear();

	// Destroys segments
	std::map<int, Polygon_Segment*>::iterator it_po_s;
	std::map<int, Planar_Segment*>::iterator it_pl_s;
	
	for (it_po_s = segments.begin() ; it_po_s != segments.end() ; it_po_s++) {
		delete it_po_s->second;
	}

	for (it_pl_s = borders.begin() ; it_pl_s != borders.end() ; it_pl_s++) {
		delete it_pl_s->second;
	}

	// Destroys lines
	std::map<int, Intersection_Line*>::iterator it_l;
	for (it_l = lines.begin() ; it_l != lines.end() ; it_l++) {
		delete it_l->second;
	}
}



CGAL_Point_2 Support_Plane::project(const CGAL_Point_3 & P) const
{
	// Assuming that P is a point on this supporting plane,
	// returns its 2D coordinates in the local 2D frame.

	return plane.to_2d(P);
}


CGAL_Point_3 Support_Plane::backproject(const CGAL_Point_2 & P) const
{
	// Assuming that P is a point on this supporting plane,
	// returns its 3D coordinates in the global 3D frame.

	return plane.to_3d(P);
}



void Support_Plane::init_intersection_lines(const std::map<int, CGAL_Line_3> & lines_3d)
{
	for (std::map<int, CGAL_Line_3>::const_iterator it_l = lines_3d.begin() ; it_l != lines_3d.end() ; it_l++) {

		// Gets a definition of the 3D line
		const CGAL_Point_3 A = it_l->second.point();
		const CGAL_Point_3 B = A + it_l->second.to_vector();

		// Gets a definition of the projected 2D line
		const CGAL_Point_2 A_proj = plane.to_2d(A);
		const CGAL_Point_2 B_proj = plane.to_2d(B);
		CGAL_Line_2 L (A_proj, B_proj);

		// Checks if this line is equal to another line
		// i.e. if their directions are collinear, and a point of one line belongs to the other
		std::map<int, Intersection_Line*>::iterator it_m;
		for (it_m = lines.begin() ; it_m != lines.end() ; it_m++) {
			const CGAL_Line_2 & M = it_m->second->line;
			if (CGAL::determinant(L.to_vector(), M.to_vector()) == 0 && (M.a() * A_proj.x() + M.b() * A_proj.y() + M.c() == 0)) {
				break;
			}
		}

		if (it_m != lines.end()) {
			// This line is therefore shared by another plane
			it_m->second->mark_as_intersected(it_l->first);
			std::cerr << "Warning : P[" << id << "] : the intersection of three planes is a line" << std::endl;

		} else {
			// We are going to insert a new line, but before that we consider K_proj = A_proj + v_ab.
			// If K = backproject(K_proj) is above L and below the plane represented by L, or vice-versa,
			// then we revert the coefficients of L.

			const CGAL_Point_2 K_proj = A_proj + CGAL_Vector_2(L.a(), L.b());
			const CGAL_Point_3 K = plane.to_3d(K_proj);
			const CGAL_Plane & P = Universe::map_of_planes[it_l->first]->plane;

			const FT k_l = L.a() * K_proj.x() + L.b() * K_proj.y() + L.c();
			const FT k_p = P.a() * K.x() + P.b() * K.y() + P.c() * K.z() + P.d();
			if (k_l * k_p < 0) L = CGAL_Line_2(-L.a(), -L.b(), -L.c());

			Intersection_Line* I = new Intersection_Line(id, L, it_l->first);
			lines[I->id_object] = I;
		}
	}
}



Intersection_Line* Support_Plane::get_line(const int id_plane)
{
	for (std::map<int, Intersection_Line*>::iterator it_l = lines.begin(); it_l != lines.end(); it_l++) {
		if (it_l->second->intersects(id_plane)) {
			return it_l->second;
		}
	}
	return nullptr;
}



void Support_Plane::get_concurrent_lines(const Intersection_Line* I_1, const Intersection_Line* I_2, std::list<Intersection_Line*> & L)
{
	int p = jmin(I_1->id_object, I_2->id_object), q = jmax(I_1->id_object, I_2->id_object);
	L.clear();

	// All triplets (t_0, t_1, t_2) are such that t_0 < t_1 < t_2, and we know that p < q.
	// We make use of this representation to select the interesting triplets and stop the search
	// once t_0 > p.

	for (std::set<Triplet, Triplet_Comparator>::const_iterator it_t = concurrences.begin() ; it_t != concurrences.end() ; it_t++) {
		const Triplet & T = (*it_t);

		const int t_0 = std::get<0>(T), t_1 = std::get<1>(T), t_2 = std::get<2>(T);

		if (t_0 < p) {
			if (t_1 == p && t_2 == q) {
				L.push_back(lines[t_0]);
			}
		} else if (t_0 == p) {
			if (t_1 < q) {
				if (t_2 == q) L.push_back(lines[t_1]);
			} else if (t_1 == q) {
				L.push_back(lines[t_2]);
			} else {
				// t_1 > q : Nothing else to search for
				return;
			}
		} else {
			return;
		}
	}
}



void Support_Plane::get_indices_of_intersecting_planes(const Intersection_Line* I_1, const Intersection_Line* I_2, std::list<int> & P)
{
	// We consider the point located at the intersection of this plane, I_1 and I_2,
	// and we determine the indices of all the planes that intersect in that location.

	P.clear();
	P.push_back(id);

	// If I_1 and I_2 represent more than one plane, then indices of such planes are added to the result.
	std::copy(I_1->planes.begin(), I_1->planes.end(), std::back_inserter(P));
	std::copy(I_2->planes.begin(), I_2->planes.end(), std::back_inserter(P));

	// However, if there exists another line I_3 which intersects I_1 at the same location as I_2,
	// then planes represented by I_3 are also added to the list P.

	std::list<Intersection_Line*> I_L;
	get_concurrent_lines(I_1, I_2, I_L);

	for (std::list<Intersection_Line*>::iterator it_l = I_L.begin() ; it_l != I_L.end() ; it_l++) {
		Intersection_Line* L = (*it_l);
		std::copy(L->planes.begin(), L->planes.end(), std::back_inserter(P));
	}

	// Sorts indices
	std::vector<int> V_P;
	V_P.reserve(P.size());

	std::copy(P.begin(), P.end(), std::back_inserter(V_P));
	std::sort(V_P.begin(), V_P.end());

	// Moves back sorted indices to P
	P.clear();
	std::copy(V_P.begin(), V_P.end(), std::back_inserter(P));
}



void Support_Plane::insert_triplets_of_concurrent_lines(const std::vector<Intersection_Line*> & I)
{
	size_t n = I.size();

	std::vector<int> J;
	J.reserve(n);

	// Gets a sorted vector of indices of concurrent lines

	for (int i = 0; i < n; i++) J.push_back(I[i]->id_object);
	std::sort(J.begin(), J.end());

	// From J, we extract a list of possible triplets
	// We test if the first triplet is already present in the list of triplets of concurrent lines
	// If so, then it's not necessary to insert any of the elements, since this intersection of n lines has already been processed before
	
	if (concurrences.find(std::make_tuple(J[0], J[1], J[2])) != concurrences.end()) return;

	for (int i = 0 ; i < n - 2 ; i++) {
		for (int j = i + 1 ; j < n - 1 ; j++) {
			for (int k = j + 1 ; k < n ; k++) {
				concurrences.insert(std::make_tuple(J[i], J[j], J[k]));
			}
		}
	}
}



void Support_Plane::init_bounding_polygon(const std::list<CGAL_Point_3> & bounding_polygon_3d, const std::vector<std::list<int> > & bounding_facets)
{
	// Converts the coordinates of the 3D polygon to the 2D local frame
	// and finds the coordinates of a squared bounding box for the support plane
	
	std::vector<CGAL_Point_2> points;

	FT ft_x_min = FT(FLT_MAX), ft_x_max = FT(-FLT_MAX);
	FT ft_y_min = ft_x_min, ft_y_max = ft_x_max;

	for (std::list<CGAL_Point_3>::const_iterator it_p = bounding_polygon_3d.begin() ; it_p != bounding_polygon_3d.end() ; it_p++) {
		CGAL_Point_2 pt = project(*it_p);

		FT x = pt.x(), y = pt.y();
		if (x < ft_x_min) ft_x_min = x;
		if (x > ft_x_max) ft_x_max = x;
		if (y < ft_y_min) ft_y_min = y;
		if (y > ft_y_max) ft_y_max = y;

		points.push_back(pt);
	}

	x_min = to_double(ft_x_min), x_max = to_double(ft_x_max);
	y_min = to_double(ft_y_min), y_max = to_double(ft_y_max);

	// Builds the polygon
	int n = int(points.size());

	for (int i = 0 ; i < n ; i++) {

		// We are going to define the Planar_Segment [P_i, P_{i + 1}]
		// It is attached to the facets whose indices are listed at BF[i].
		// In case BF[i] is a list of size greater than 1, segments are duplicated.
		for (std::list<int>::const_iterator it_f = bounding_facets[i].begin() ; it_f != bounding_facets[i].end() ; it_f++) {
			int f = (*it_f);

			Intersection_Line* L = get_line(f);
			Planar_Segment* s_f = new Planar_Segment(id, L, points[i], points[(i + 1) % n]);
		}
	}

	// Discards the lines that don't intersect the bounding polygon
	for (std::map<int, Intersection_Line*>::iterator it_l = lines.begin() ; it_l != lines.end() ; it_l++) {
		Intersection_Line* I = it_l->second;
		const CGAL_Line_2 & L = I->line;

		bool positive_vertices_exist = false;
		bool negative_vertices_exist = false;

		// We are going to test if the points that define the bounding polygon are located on both sides of L.
		// If not, it can be discarded, because no Polygon_Vertex will ever reach it

		for (int i = 0 ; i < n ; i++) {
			const CGAL_Point_2 & M = points[i];
			FT epsilon = L.a() * M.x() + L.b() * M.y() + L.c();

			if (epsilon == 0) {
				// L intersects the extrema of the bouding polygon so we already know it shouldn't be discarded
				positive_vertices_exist = negative_vertices_exist = true;
				break;
			}

			if (epsilon > 0) {
				positive_vertices_exist = true;
			} else {
				negative_vertices_exist = true;
			}

			if (positive_vertices_exist && negative_vertices_exist) break;
		}

		I->set_inside(positive_vertices_exist && negative_vertices_exist);
	}

	// Sets a vector with lines that are inside the bounding polygon
	lines_inside.clear();
	lines_inside.reserve(lines.size());
	for (std::map<int, Intersection_Line*>::const_iterator it_l = lines.begin() ; it_l != lines.end() ; ++it_l) {
		if (it_l->second->is_inside) {
			lines_inside.push_back(it_l->second);
		}
	}
	lines_inside.shrink_to_fit();
}



#if 0
void Support_Plane::init_rtree_lines()
{
	const double eps = 1e-3;
	uint b_index = -1;

	// In this function, we are going to insert N * rtree_boxes_per_line in a R-Tree,
	// where N represents the number of lines that intersect the bounding box.
	
	rtree_boxes_per_line = int(ceil(jmin((x_max - x_min) / rtree_boxes_width, (y_max - y_min) / rtree_boxes_width)));

	// The purpose of this operation comes later :
	// when we will schedule the intersections of a Polygon_Vertex v,
	// we can predict the collisions that may occur between [t_k, t_{k + 1}] by searching elements in this tree.

	// Boxes that are inserted in this tree have an index ('b_index', incremented at each insertion).
	// The conversion between box indices and lines is given by the table 'rtree_indices' :
	// rtree_indices[b_index / rtree_boxes_per_line] returns the associated line L.
	
	std::list<Intersection_Line*> discretized_lines;

	for (std::map<int, Intersection_Line*>::iterator it_l = lines.begin() ; it_l != lines.end() ; it_l++) {
		Intersection_Line* L = it_l->second;

		if (!L->is_inside) continue;
		const CGAL_Line_2 & l_eq = L->line;

		// Step 1.
		// For every line, we find its intersection points A and B with the bounding polygon of this support plane.
		bool A_exists = false, B_exists = false;
		CGAL_Point_2 A, B;

		for (std::map<int, Planar_Segment*>::iterator it_pl_s = borders.begin(); it_pl_s != borders.end(); it_pl_s++) {
			Planar_Segment* pl_s = it_pl_s->second;
			CGAL_Segment_2 S (pl_s->A, pl_s->B);

			CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Segment_2)>::type object = intersection(l_eq, S);
			if (object) {
				if (const CGAL_Point_2* ptr = boost::get<CGAL_Point_2>(&*object)) {
					if (!A_exists) {
						A = *ptr;
						A_exists = true;
					} else if (A != *ptr) {
						B = *ptr;
						B_exists = true;
					}
				} else if (const CGAL_Segment_2* ptr = boost::get<CGAL_Segment_2>(&*object)) {
					A = ptr->source();
					B = ptr->target();
					A_exists = B_exists = true;
				}
			}

			// Once we've found the intersection of the line l_eq and the bounding box, breaks the loop
			if (A_exists && B_exists) break;
		}
		
		if (!A_exists) {
			throw std::logic_error("Error : couldn't find the intersection of a line and a bounding polygon");
		} else if (!B_exists) {
			B = A;
		}
		
		// Step 2.
		// Discretizes the segment [AB] in thin boxes

		// We are going to consider N + 1 points (M_k) uniformly distributed on [AB],
		// each couple of points (M_k M_{k + 1}) defines directly or indirectly a box,
		// depending on the orientation of the segment [AB].

		CGAL_Vector_2 AB = B - A;
		const double x_ab = to_double(AB.x()), y_ab = to_double(AB.y());
		const double zero = 1e-12;

		if (fabs(x_ab) > zero && fabs(y_ab) > zero) {

			if (x_ab * y_ab > 0) {

				// The segment AB is oriented like this -> /
				// We consider the bottom left and top right corners and define the points M_k

				CGAL_Point_2 BL, TR;
				if (x_ab > 0) {
					BL = A, TR = B;
				} else {
					BL = B, TR = A;
				}

				const double x_bl = to_double(BL.x()), y_bl = to_double(BL.y());
				const double x_tr = to_double(TR.x()), y_tr = to_double(TR.y());
				const double dx = x_tr - x_bl, dy = y_tr - y_bl;
				assert(dx > 0 && dy > 0);

				double x_prev = x_bl, y_prev = y_bl, x_curr, y_curr;
				for (uint k = 1; k <= rtree_boxes_per_line - 1; k++) {
					x_curr = x_bl + k * dx / rtree_boxes_per_line;
					y_curr = y_bl + k * dy / rtree_boxes_per_line;

					Boost_Box K(Boost_Point(x_prev - eps, y_prev - eps), Boost_Point(x_curr + eps, y_curr + eps));
					rtree_lines.insert(std::make_pair(K, ++b_index));

					x_prev = x_curr;
					y_prev = y_curr;
				}

				Boost_Box K(Boost_Point(x_prev - eps, y_prev - eps), Boost_Point(x_tr + eps, y_tr + eps));
				rtree_lines.insert(std::make_pair(K, ++b_index));

			} else if (x_ab * y_ab < 0) {

				// The segment AB is oriented like this -> \
				// We find the bottom right and top left corners of the segment

				CGAL_Point_2 BR, TL;
				if (y_ab > 0) {
					BR = A, TL = B;
				} else {
					BR = B, TL = A;
				}

				const double x_br = to_double(BR.x()), y_br = to_double(BR.y());
				const double x_tl = to_double(TL.x()), y_tl = to_double(TL.y());
				const double dx = x_tl - x_br, dy = y_tl - y_br; // dx < 0, dy > 0
				assert(dx < 0 && dy > 0);

				double x_prev = x_br, y_prev = y_br, x_curr, y_curr;
				for (uint k = 1; k <= rtree_boxes_per_line - 1; k++) {
					x_curr = x_br + k * dx / rtree_boxes_per_line;
					y_curr = y_br + k * dy / rtree_boxes_per_line;

					Boost_Box K(Boost_Point(x_curr - eps, y_prev - eps), Boost_Point(x_prev + eps, y_curr + eps));
					rtree_lines.insert(std::make_pair(K, ++b_index));

					x_prev = x_curr;
					y_prev = y_curr;
				}

				Boost_Box K(Boost_Point(x_tl - eps, y_prev - eps), Boost_Point(x_prev + eps, y_tl + eps));
				rtree_lines.insert(std::make_pair(K, ++b_index));

			}

		} else if (fabs(x_ab) > zero && fabs(y_ab) < zero) {
		
			// The segment AB is oriented like this -> --
			// We find the left and right points
			
			CGAL_Point_2 L, R;
			if (x_ab > 0) {
				L = A, R = B;
			} else {
				L = B, R = A;
			}

			const double y = to_double(L.y());
			const double x_l = to_double(L.x()), x_r = to_double(R.x());
			const double dx = x_r - x_l;
			assert(dx > 0);

			double x_prev = x_l, x_curr;
			for (uint k = 1; k <= rtree_boxes_per_line - 1; k++) {
				x_curr = x_l + k * dx / rtree_boxes_per_line;
				Boost_Box K (Boost_Point(x_prev - eps, y - eps), Boost_Point(x_curr + eps, y + eps));
				rtree_lines.insert(std::make_pair(K, ++b_index));
				x_prev = x_curr;
			}

			Boost_Box K (Boost_Point(x_prev - eps, y - eps), Boost_Point(x_r + eps, y + eps));
			rtree_lines.insert(std::make_pair(K, ++b_index));

		} else if (fabs(y_ab) > zero && fabs(x_ab) < zero) {
		
			// The segment AB is oriented like this -> |
			// We identify the top and bottom points

			CGAL_Point_2 D, U;
			if (y_ab > 0) {
				D = A, U = B;
			} else {
				D = B, U = A;
			}

			const double x = to_double(D.x());
			const double y_d = to_double(D.y()), y_u = to_double(U.y());
			const double dy = y_u - y_d;
			assert(dy > 0);

			double y_prev = y_d, y_curr;
			for (uint k = 1; k <= rtree_boxes_per_line - 1; k++) {
				y_curr = y_d + k * dy / rtree_boxes_per_line;
				Boost_Box K (Boost_Point(x - eps, y_prev - eps), Boost_Point(x + eps, y_curr + eps));
				rtree_lines.insert(std::make_pair(K, ++b_index));
				y_prev = y_curr;
			}

			Boost_Box K (Boost_Point(x - eps, y_prev - eps), Boost_Point(x + eps, y_u + eps));
			rtree_lines.insert(std::make_pair(K, ++b_index));

		} else {
		
			// Case when the intersection of the line and the bounding polygon is just a point
			// Yet we are supposed to have N boxes along a line
			// As it would be dangerous to ignore this line for our further computations,
			// we have no choice but duplicating the same bounding box

			const double x = to_double(A.x()), y = to_double(A.y());
			for (uint k = 0; k < rtree_boxes_per_line ; k++) {
				Boost_Box K (Boost_Point(x - eps, y - eps), Boost_Point(x + eps, y + eps));
				rtree_lines.insert(std::make_pair(K, ++b_index));
			}
		}

		discretized_lines.push_back(L);
	}

	// Cf. supra
	rtree_boxes_to_lines = std::vector<Intersection_Line*>(discretized_lines.begin(), discretized_lines.end());
}



void Support_Plane::insert_in_rtree(uint & index, double & x_prev, double & y_prev, double & x_curr, double & y_curr, const double & eps)
{
	// Inserts a box in the rtree with coordinates (x_prev, y_prev, x_curr, y_curr)

	Boost_Box B (Boost_Point(x_prev - eps, y_prev - eps), Boost_Point(x_curr + eps, y_curr + eps));
	rtree_lines.insert(std::make_pair(B, ++index));

	// Iterates
	x_prev = x_curr;
	y_prev = y_curr;
}
#endif



void Support_Plane::search_lines_in_neighborhood(Polygon_Vertex* v, const FT & t_1, const FT & t_2, 
	std::list<std::tuple<Intersection_Line*, bool, FT> > & L)
{
	// As Polygon_Vertex v propagates, we want to answer the following question :
	// which Intersection_Lines are susceptible to get intersected by it between t_1 and t_2 ?

	// Step 1
	// Finds the coordinates of the bounding box for v between t_1 and t_2

	CGAL_Point_2 M_1 = v->pt(t_1);
	CGAL_Point_2 M_2 = v->pt(t_2);

	const double x_1 = to_double(M_1.x()), x_2 = to_double(M_2.x()), y_1 = to_double(M_1.y()), y_2 = to_double(M_2.y());


#if 1
	for (int j = 0 ; j < lines_inside.size() ; j++) {
		Intersection_Line* I_j = lines_inside[j];
		
		/*std::pair<FT, bool> R = v->get_intersection_time(I_j);
		if (R.second) {
			if (R.first >= t_1 && R.first < t_2) {
				L.push_back(std::make_tuple(I_j, true, R.first));
			}
		}*/

		const double & h_a = I_j->hint_a, & h_b = I_j->hint_b, & h_c = I_j->hint_c;
		const double dh_1 = h_a * x_1 + h_b * y_1 + h_c;
		const double dh_2 = h_a * x_2 + h_b * y_2 + h_c;

		if (dh_1 * dh_2 > 0.1) {
			continue;
		} else if (dh_1 * dh_2 < -0.1) {
			L.push_back(std::make_tuple(I_j, false, 0));
		} else {
			std::pair<FT, bool> R = v->get_intersection_time(I_j);
			if (R.second) {
				if (R.first >= t_1 && R.first < t_2) {
					L.push_back(std::make_tuple(I_j, true, R.first));
				}
			}
		}
	}
#else
	double bx_min, bx_max, by_min, by_max;

	if (x_1 < x_2) {
		bx_min = x_1; bx_max = x_2;
	} else {
		bx_min = x_2; bx_max = x_1;
	}

	if (y_1 < y_2) {
		by_min = y_1; by_max = y_2;
	} else {
		by_min = y_2; by_max = y_1;
	}

	// Step 2
	// Processes a query, only keeps the indices of the intersected boxes

	Boost_Box query(Boost_Point(bx_min, by_min), Boost_Point(bx_max, by_max));
	std::vector<Boost_Value> values;
	
	rtree_lines.query(bgi::intersects(query), std::back_inserter(values));
	std::vector<bool> lines_found = std::vector<bool>(rtree_boxes_to_lines.size(), false);

	for (size_t i = 0 ; i < values.size() ; i++) {
		uint b_index = values[i].second;
		lines_found[b_index / rtree_boxes_per_line] = true;
	}

	for (size_t i = 0 ; i < lines_found.size() ; i++) {
		if (lines_found[i]) L.push_back(rtree_boxes_to_lines[i]);
	}
#endif
}



void Support_Plane::init_polygon_set()
{
	polygon_set = new Polygon_Set(lines);
}



void Support_Plane::init_polygon(const std::vector<CGAL_Point_3> & polygon)
{
	// Gets a polygon and splits it hierarchically
	Polygon_Tree* T = decompose_and_init_segments(polygon);
	if (T == nullptr) return;

	std::list<Polygon*> LP;
	T->get_polygons(LP);

	// Loops on all the initial polygons
	// Each polygon is used to initialize a Polygon_Cell structure
	for (std::list<Polygon*>::iterator it_p = LP.begin() ; it_p != LP.end() ; it_p++) {
		Polygon* P = (*it_p);

		// Gets a signature for the polygon
		const CGAL_Point_2 M = P->get_barycenter(0);
		std::vector<bool> S = Polygon_Cell::make_signature(M, lines);

		// Inserts a couple (S, P) in the set of polygon cells
		polygon_set->insert(S, P);
	}

	T->remove_reference_to_polygons();
	delete T;
}



void Support_Plane::set_initial_propagation_directions(const std::vector<CGAL_Point_2> & R, CGAL_Point_2 & O, std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & D) const
{
	// Computes the barycenter of all points
	
	const int n = int(R.size());
	FT x_bar = 0, y_bar = 0;

	for (int i = 0 ; i < n ; i++) {
		const CGAL_Point_2 & pt = R[i];
		x_bar += pt.x(), y_bar += pt.y();
	}
	O = CGAL_Point_2(x_bar / n, y_bar / n);
	
	// The directions of all the points M[i] are given by the vector OM[i]
	// We finally get pairs of (M[i], OM[i]) that we store in the vector D

	D.reserve(n);
	for (int i = 0 ; i < n ; i++) {
		D.push_back(std::make_pair(R[i], R[i] - O));
	}
}



void Support_Plane::project_polygon(const std::vector<CGAL_Point_3> & P_0, std::vector<CGAL_Point_2> & P) const
{
	std::vector<CGAL_Point_2> Q;
	Q.reserve(P_0.size());

	for (int i = 0; i < int(P_0.size()); i++) {
		Q.push_back(project(P_0[i]));
	}

	CGAL::convex_hull_2(Q.begin(), Q.end(), std::back_inserter(P));
}



Polygon_Tree* Support_Plane::decompose_and_init_segments(const std::vector<CGAL_Point_3> & polygon)
{
	// Part 1.
	// Initializes a 2D polygon

	// Projects and regularizes polygon
	std::vector<CGAL_Point_2> P;
	std::vector<Intersection_Line*> reference_lines;
	project_polygon(polygon, P);

	// Gets the locations and directions of all initial vertices
	CGAL_Point_2 initial_barycenter;
	std::vector< std::pair<CGAL_Point_2, CGAL_Vector_2> > initial_directions;
	set_initial_propagation_directions(P, initial_barycenter, initial_directions);

	// Builds an initial polygon, itself used to initialize a structure of polygonal tree
	// This tree will be later hierarchically decomposed, depending on possible intersections with Intersection_Lines

	int seed = int(polygon_directions.size());
	Polygon_Directions* D = new Polygon_Directions(initial_barycenter, initial_directions, reference_lines);
	polygon_directions.push_back(D);

	Polygon_Tree* polygon_tree = new Polygon_Tree(id, seed, D);


	// Part 2.
	// Divides the polygon into subpolygons it is intersects with any of the intersection lines.
	// In case of intersection, we keep trace of the constrained vertices that delimit the intersection.
	// Later, we will initialize segments using them.

	typedef std::tuple<Polygon_Vertex*, Polygon_Vertex*, Polygon_Vertex*, Polygon_Vertex*> Quadruplet;

	std::map<int, Quadruplet> polyline_intersections;
	std::map<int, std::map<int, CGAL_Point_2> > mutual_intersections;

	for (std::map<int, Intersection_Line *>::const_iterator it_l = lines.begin(); it_l != lines.end(); it_l++) {
		Intersection_Line* I = it_l->second;
		
		int i = I->id_object;
		std::list<Polygon_Vertex *> intersection_pts;

		if (polygon_tree->split(I, intersection_pts, 0, Universe::params->K, NO_SCHEDULE, true)) {

			// If the tree had not been divided by another line during the initialization process,
			// or if only one leaf of the tree is split by the line, then 'intersection' only contains 4 active vertices.
			// (2 intersections points, that are duplicated since there are 2 sides to take into account.)

			// Otherwise it contains several 4-uplets of vertices that correspond to the intersections of the line and each of the intersected subpolygons of the tree. 
			// In the latter case, only 4 vertices are active (others correspond to the intersections of two lines, that's why they are disabled).

			Polygon_Vertex *v1_p = nullptr, *v1_n = nullptr, *v2_p = nullptr, *v2_n = nullptr;

			for (std::list<Polygon_Vertex *>::iterator it_v = intersection_pts.begin(); it_v != intersection_pts.end(); it_v++) {
				Polygon_Vertex* v = (*it_v);

				if (!v->is_active) {
					// If the vertex is not active, then it is at the intersection of I and another line H processed before
					// If we haven't met it yet, then we mark down it location in map crossed_lines_I
					// and we notify the line H that an intersection with I exists
					int h = v->get_constraint().first->id_object;
					if (h == i) h = v->get_second_constraint().first->id_object;
					if (mutual_intersections[i].find(h) == mutual_intersections[i].end()) {
						mutual_intersections[i][h] = v->M;
						mutual_intersections[h][i] = v->M;
					}
					continue;
				}

				// We find two pairs of opposite vertices
				if (v->get_constraint().second == PLUS) {
					set_element_in_quadruplet(v1_n, v2_n, v1_p, v2_p, v);
				} else {
					set_element_in_quadruplet(v1_p, v2_p, v1_n, v2_n, v);
				}
			}

			polyline_intersections[i] = std::make_tuple(v1_p, v1_n, v2_p, v2_n);
		}
	}

	// Part 3.
	// Initializes segments

	for (std::map<int, Quadruplet>::const_iterator it_pli = polyline_intersections.cbegin() ; it_pli != polyline_intersections.cend() ; it_pli++) {
		int i = it_pli->first;

		// Exhibits the 4-uplet of vertices

		const Quadruplet & Q = it_pli->second;
		Polygon_Vertex *v1_p = std::get<0>(Q), *v1_n = std::get<1>(Q), *v2_p = std::get<2>(Q), *v2_n = std::get<3>(Q);
		Intersection_Line* I = v1_p->get_constraint().first;

		bool v1_paired = Polygon_Vertex::are_paired_vertices(v1_p, v1_n);
		bool v2_paired = Polygon_Vertex::are_paired_vertices(v2_p, v2_n);

		const CGAL_Point_2 & V1 = v1_p->M, &V2 = v2_p->M;
		const CGAL_Point_3 W1 = backproject(V1), W2 = backproject(V2);

		const std::map<int, CGAL_Point_2> & L = mutual_intersections[i];

		if (L.empty()) {
			// If line I doesn't intersect any other line,
			// then we initialize a 4-uplet of opposite bidirectional segments without dilimitation
			// w.r.t. the other lines of the support plane.
			if (v1_paired && v2_paired) {
				init_4_uplets_of_bidirectional_segments(I, V1, W1, V2, W2, v1_p, v1_n, v2_p, v2_n, 0);
			} else {
				init_2_uplets_of_bidirectional_segments(I, V1, W1, V2, W2, v1_p, v2_p, 0);
				init_2_uplets_of_bidirectional_segments(I, V1, W1, V2, W2, v1_n, v2_n, 0);
			}

		} else {
			// I intersects (I_1 .. I_n).
			// We use I_1 to initialize unidirectional segments with an initial constraint,
			// and we assign the other lines (I_2 .. I_n) to one of the 2 2-uplets of segments thus created,
			// by marking them as intersected by the segments.

			// Sets J, which plays the role of I_1.
			// It intersects the line I in M_ij.
			std::map<int, CGAL_Point_2>::const_iterator it_l = L.cbegin();
			Intersection_Line* J = lines[it_l->first];
			const CGAL_Point_2 & M_ij = it_l->second;
			const CGAL_Point_3 W_ij = backproject(M_ij);

			Sign eps_1 = J->sign(v1_p->M);
			Sign eps_2 = (eps_1 == PLUS ? MINUS : PLUS);
			Constraint C_1(J, eps_1), C_2(J, eps_2);

			// We initialize two 2-uplets of segments [M_ij v1], [M_ij v2]
			// J is used to set the initial constraint of these segments.

			if (v1_paired) {
				init_2_uplets_of_unidirectional_segments(I, M_ij, W_ij, V1, W1, v1_p, v1_n, C_1, 0);
			} else {
				init_1_uplets_of_unidirectional_segments(I, M_ij, W_ij, V1, W1, v1_p, C_1, 0);
				init_1_uplets_of_unidirectional_segments(I, M_ij, W_ij, V1, W1, v1_n, C_1, 0);
			}

			if (v2_paired) {
				init_2_uplets_of_unidirectional_segments(I, M_ij, W_ij, V2, W2, v2_p, v2_n, C_2, 0);
			} else {
				init_1_uplets_of_unidirectional_segments(I, M_ij, W_ij, V2, W2, v2_p, C_2, 0);
				init_1_uplets_of_unidirectional_segments(I, M_ij, W_ij, V2, W2, v2_n, C_2, 0);
			}

			while (++it_l != L.cend()) {
				// Sets sucessive lines K, which play the roles of I_2 .. I_n.
				// They also intersect I in a set of points M_ik, 
				// we need to find if its intersects [M_ij v1], [M_ij v2], or both.
				Intersection_Line* K = lines[it_l->first];
				const CGAL_Point_2 & M_ik = it_l->second;

				int r = locate_intersection(M_ij, V1, V2, M_ik);
				if (r == 0) {
					// Intersects both segments [M_ij v1] and [M_ij v2]
					// Nothing to do ?
					std::cerr << "Warning : not implemented case in the initialization process -- no guarantee of result" << std::endl;
				} else if (r == 1) {
					// Intersects segment [M_ij v1]
					v1_p->indicate_line_initially_crossed_by_segments(K);
				} else {
					// Intersects segment [M_ij v2]
					v2_p->indicate_line_initially_crossed_by_segments(K);
					
				}
			}

			v1_p->copy_crossed_lines(v1_n);
			v2_p->copy_crossed_lines(v2_n);
		}
	}

	return polygon_tree;
}



int Support_Plane::locate_intersection(const CGAL_Point_2 & M_0, const CGAL_Point_2 & M_1, const CGAL_Point_2 & M_2, const CGAL_Point_2 & P) const
{
	// By calling this function, we have two segments s_1 = [M_1 M_0], s_2 = [M_0 M_2]
	// and we want to determine if the point P belongs to s_1, s_2 or both.
	// We suppose that P belongs to [M_1 M_2].

	const FT x_0 = M_0.x(), x_1 = M_1.x(), x_2 = M_2.x();
	if (x_1 != x_2) {
		// [M_1 M_2] is not vertical : the x coordinate can be used.
		// Three possibilities :
		// - P.x() = x_0 => P = M_0 => return 0
		// - P.x() is in ]x_0 x_1]  => return 1
		// - P.x() is in ]x_0 x_2]  => return 2
		const FT x = P.x();
		if (x_0 == x) {
			return 0;
		} else if ((x_0 <= x && x <= x_1) || (x_1 <= x && x <= x_0)) {
			return 1;
		} else {
			return 2;
		}
	} else {
		// [M_1 M_2] is vertical.
		// We should use the y coordinate, but with a similar reasoning.
		const FT y_0 = M_0.y(), y_1 = M_1.y(), y = P.y();
		if (y_0 == y) {
			return 0;
		} else if ((y_0 <= y && y <= y_1) || (y_1 <= y && y <= y_0)) {
			return 1;
		} else {
			return 2;
		}
	}
}



void Support_Plane::set_element_in_quadruplet(Polygon_Vertex* & v1_os, Polygon_Vertex* & v2_os, Polygon_Vertex* & v1_ts, Polygon_Vertex* & v2_ts, Polygon_Vertex* v) const
{
	if (v1_os != nullptr && v2_os == nullptr) {
		if ((v1_os->M - v->M) == CGAL::NULL_VECTOR) {
			v1_ts = v;
			return;
		}
	} else if (v1_os == nullptr && v2_os != nullptr) {
		if ((v2_os->M - v->M) == CGAL::NULL_VECTOR) {
			v2_ts = v;
			return;
		}
	} else if (v1_os != nullptr && v2_os != nullptr) {
		if ((v1_os->M - v->M).squared_length() < (v2_os->M - v->M).squared_length()) {
			v1_ts = v;
			return;
		} else {
			v2_ts = v;
			return;
		}
	}
	if (v1_ts == nullptr) {
		v1_ts = v;
		return;
	} else if (v2_ts == nullptr && (v1_ts->dM * v->dM < 0)) {
		v2_ts = v;
		return;
	}
}



void Support_Plane::init_schedule()
{
	// Loops on the list of all active vertices and schedule events

	for (std::map<int, Polygon_Vertex*>::iterator it_v = vertices.begin() ; it_v != vertices.end() ; it_v++) {
		Polygon_Vertex* v = it_v->second;

		if (v->is_active) {
			// Two possible cases :
			// - either the vertex is independant, and we compute its events normally,
			// - or it is paired to another vertex, and in this case, we duplicate the events of its 'master' vertex

			if (v->is_independent()) {
				if (v->unconstrained()) {
					v->schedule_events();
				} else {
					const Constraint & C = v->get_constraint();
					Intersection_Line* I = C.first;
					v->schedule_events(I);
				}
			}

			else {
				Polygon_Vertex* v_ts = v->get_master_vertex();
				v->schedule_events(v_ts);
			}
		}
	}
}



void Support_Plane::process_event(Event_Vertex_Line* e)
{
	clock_t t_begin = clock();

	int intersectant = e->intersectant;
	int intersected = e->intersected;
	const FT t_intersectant = e->t_intersectant;

	Polygon_Vertex* v = vertices[intersectant];
	Intersection_Line* I_0 = (v->unconstrained() ? nullptr : v->get_constraint().first);

	CGAL_Point_2 V_t = v->pt(t_intersectant);

	bool verbose = false;

	// Step 1.
	// In this function, we would like to associate the correct handle to the event e that has been popped from the queue.
	// We first factorize events, by obtaining the list of lines I intersected by vertex v at time t.

	std::list<Intersection_Line*> I;
	std::list<Event_Vertex_Line*> E;
	get_simultaneous_events(v, t_intersectant, V_t, E, e, I);
	v->decrement_queued_events(E.size());

	// Step 2.
	// Performs different tests, to associate the current handle to E_VL.

	try {
		if (Polygon_Vertex* v_n = v->get_constrained_neighbor(I_0, I)) {
			if (I_0 == nullptr) {
				// Case B.
				// An unconstrained vertex v intersects a set of lines I_k,
				// and meets at the same time one of its neighbors v_n, constrained by one of these lines.
				// The vertex is going to be passed to the adjacent polygon, if it exists.
				if (verbose) std::cout << "[" << id << "] " << e->intersectant << " " << e->intersected << " " << e->t_intersectant << " B" << std::endl;
				unconstrained_vertex_intersects_line_kth_time(E, v, V_t);

			} else {
				// Case C1.
				// A constrained vertex intersects a set of lines I_k,
				// and meets at the same time one of its neighbors, constrained by one of these lines.
				// The two constrained vertices should be merged into a double-constrained and still vertex.
				if (verbose) std::cout << "[" << id << "] " << e->intersectant << " " << e->intersected << " " << e->t_intersectant << " C1" << std::endl;
				constrained_vertices_intersect(E, v, v_n, V_t);
			}

		} else {
			if (Polygon_Vertex* v_p = v->get_neighbor_intersecting_identical_line(I_0, I, t_intersectant)) {
				// Case D.
				// There exists a line intersected, at the same time, by two vertices v and v_p.
				// We should actually process a collision between an edge and a line.
				if (verbose) {
					std::cout << "[" << id << "] " << e->intersectant << " " << e->intersected << " " << e->t_intersectant << " D" << std::endl;
					std::cout << "[" << id << "] " << v_p->id_object  << " " << e->intersected << " " << e->t_intersectant << " D *" << std::endl;
				}
				edge_intersects_line(E, v, v_p, V_t);

			} else if (I_0 == nullptr) {
				// Case A.
				// There exists no constrained neighbor, 
				// and no vertex intersecting one of the lines intersected by v at time t.
				// There is only a unconstrained vertex, intersecting a set of lines.
				if (verbose) std::cout << "[" << id << "] " << e->intersectant << " " << e->intersected << " " << e->t_intersectant << " A" << std::endl;
				unconstrained_vertex_intersects_line(E, v, V_t);

			} else {
				// Case C2.
				// There exists no neighbor constrained by one of the lines of I,
				// and there exists no vertex intersecting of the lines intersected by v at time t.
				if (verbose) std::cout << "[" << id << "] " << e->intersectant << " " << e->intersected << " " << e->t_intersectant << " C2" << std::endl;
				constrained_vertex_intersects_line(E, v, V_t);
			}
		}

	} catch (const std::exception & except) {
		// Prints the support plane
		std::cout << "Error raised while processing [" << id << "] " << intersectant << " " << intersected << " " << t_intersectant << std::endl;
		draw(to_double(t_intersectant), 0.5, 20000, 0.5, 0.5, 10);
		throw except;
	}

	if (Universe::params->rt_check) {
		real_time_check(t_intersectant);
	}

	clock_t t_end = clock();
	KP_Stats::process_events_time += (t_end - t_begin);
	++KP_Stats::process_events_calls;
}



void Support_Plane::real_time_check(const FT & t) 
{
	// For each cell
	for (std::map<Signature, Polygon_Cell*, Vector_Bool_Comparator>::iterator it_c = polygon_set->cells_begin() ; it_c != polygon_set->cells_end() ; it_c++) {
		const Signature & S = it_c->first;
		Polygon_Cell* C = it_c->second;

		// For each polygon contained in the cell
		for (std::list<Polygon*>::iterator it_p = C->polygons_begin() ; it_p != C->polygons_end() ; it_p++) {
			Polygon* P = (*it_p);
			for (std::list<Polygon_Vertex*>::iterator it_v = P->vertices.begin() ; it_v != P->vertices.end() ; it_v++) {
				Polygon_Vertex* v = (*it_v);
				if (!v->is_active) continue;

				// Checks that a vertex hasn't silently crossed a line
				CGAL_Point_2 V_t = v->pt(t);
				for (std::map<Intersection_Line*, int>::iterator it_l = polygon_set->dictionary.begin() ; it_l != polygon_set->dictionary.end() ; it_l++) {
					Intersection_Line* I = it_l->first;
					if (S[it_l->second]) {
						if (I->a() * V_t.x() + I->b() * V_t.y() + I->c() < 0) {
							std::cout << "Error with vertex : " << (*it_v)->id_object << " : missed line at t = " << t << std::endl;
							draw(to_double(t), 0.5, 20000, 0.5, 0.5, 10);
							exit(0);
						}
					} else {
						if (I->a() * V_t.x() + I->b() * V_t.y() + I->c() > 0) {
							std::cout << "Error with vertex : " << (*it_v)->id_object << " : missed line at t = " << t << std::endl;
							draw(to_double(t), 0.5, 20000, 0.5, 0.5, 10);
							exit(0);
						}
					}
				}
			}
		}
	}
}



void Support_Plane::get_simultaneous_events(Polygon_Vertex* v, const FT & t, const CGAL_Point_2 & V_t, 
	std::list<Event_Vertex_Line*> & E, Event_Vertex_Line* e, std::list<Intersection_Line*> & I)
{
	Event_Queue* Q = Universe::event_queue;

	I.clear();
	E.clear();

	// We have popped an event e from the queue
	// Now, we determine if some other events, involving the same vertex, occur simultaneously

	std::list<Event_Vertex_Line*> E_0;
	Q->get_simultaneous_events_for_this_vertex(e->intersectant, e->t_intersectant, E_0);

	for (std::list<Event_Vertex_Line*>::iterator it_e = E_0.begin() ; it_e != E_0.end() ; it_e++) {
		Event_Vertex_Line* e_p = (*it_e);
		Intersection_Line* I_p = lines[e_p->intersected];

		if (I_p->includes(V_t)) {
			// If the intersection of I_ref and v is the same as the intersection of I and v,
			// then e_ref and e occur simultaneously. We remove e from the queue and add it to E.
			// Q->erase(e_p);
			E.push_back(e_p);
			I.push_back(I_p);
		}
	}

	// Finally adds e_ref to E and I
	if (e != nullptr) {
		E.push_back(e);
		I.push_back(lines[e->intersected]);
	}
}



void Support_Plane::get_simultaneous_events_as_edge_intersects_line(Polygon_Vertex* v, const FT & t, const CGAL_Point_2 & V_t, 
	Intersection_Line* I_L, std::list<Event_Vertex_Line*> & E, Event_Vertex_Line* & e_vl)
{
	e_vl = nullptr;
	E.clear();

	Event_Queue* Q = Universe::event_queue;

	// We support that an edge e = (v0 v) intersects the line I, at time t.
	// We've already popped the event (v0 intersects I_L), we should now pop (v intersects I_L)

	Q->get_simultaneous_events_for_this_vertex(v->id_object, t, E);
	for (std::list<Event_Vertex_Line *>::iterator it_e = E.begin() ; it_e != E.end() ; it_e++) {
		if ((*it_e)->intersected == I_L->id_object) {
			e_vl = (*it_e);
			break;
		}
	}
}




void Support_Plane::append(Event_Vertex_Line* e_vl, std::string type)
{
	FILE* file = fopen("events.txt", "a+");
	if (file != NULL) {
		fprintf(file, "[%i] %i %i %lf %s\n", id, e_vl->intersectant, e_vl->intersected, to_double(e_vl->t_intersectant), type.c_str());
		fclose(file);
	}
}



void Support_Plane::get_polygon_description(std::list<std::list<CGAL_Point_3> > & P, std::list<CGAL_Color> & C, const double t)
{
	polygon_set->get_polygon_description(P, C, t);
}



void Support_Plane::group_final_polygons()
{
	// Step 1.
	// Gets a list containing all the polygons listed in the tree
	std::list<Polygon*> P;
	polygon_set->get_polygons(P);

	// Step 2.
	// Groups adjacent polygons, i.e. polygons separated by a single edge without segment on it

	int dummy = -1;

	for (std::list<Polygon*>::iterator it_p = P.begin(); it_p != P.end(); it_p++) {
		Polygon* P_curr = (*it_p);

		// We define a tuple with 3 elements :
		// - iterator to the group G compatible with P_curr
		// - list of edges of G adjacent to P_curr
		// - list of edges of P_curr adjacent to G
		typedef std::tuple<std::list<Polygon_Group*>::iterator, std::list<Polygon_Edge*>, std::list<Polygon_Edge*> > Compatibility_Info;
		std::list<Compatibility_Info> G_adj;

		for (std::list<Polygon_Group*>::iterator it_g = groups.begin() ; it_g != groups.end() ; it_g++) {
			std::list<Polygon_Edge*> e_adj_G, e_adj_P;
			if ((*it_g)->is_adjacent_to(P_curr, e_adj_G, e_adj_P)) {
				G_adj.push_back(std::make_tuple(it_g, e_adj_G, e_adj_P));
			}
		}

		// If no adjacent polygon has been found, then we create a new group of adjacent polygons with a single element
		if (G_adj.empty()) {
			groups.push_back(new Polygon_Group(P_curr));
		}

		// If a group of adjacent polygons is found, then we add P[i] to that group
		// In case there are more than one compatible groups, then we add P[i] to the first one 
		// and merge all other groups to the first one 
		else {
			std::list<Compatibility_Info>::iterator it_1 = G_adj.begin();
			Polygon_Group* G_1 = *std::get<0>(*it_1);
			G_1->append(P_curr, std::get<1>(*it_1), std::get<2>(*it_1));

			std::list<Compatibility_Info>::iterator it_k = it_1;
			while (++it_k != G_adj.end()) {
				// Case when P_curr is adjacent to more than one group
				// P_curr has previously been appended to G_1 and all its edges have been appended
				// So, as we merge G_1 and G_k, we remove from G_1 the edges of P_curr adjacent to G_k,
				// and add all borders of G_k except the edges of G_k adjacent to P_curr

				std::list<Polygon_Group*>::iterator it_gk = std::get<0>(*it_k);
				Polygon_Group* G_k = (*it_gk);

				G_1->merge(G_k, std::get<2>(*it_k), std::get<1>(*it_k));
				groups.erase(it_gk);
			}
		}
	}

	// std::cout << "id : " << id << ", |G| = " << groups.size() << std::endl;
}

}