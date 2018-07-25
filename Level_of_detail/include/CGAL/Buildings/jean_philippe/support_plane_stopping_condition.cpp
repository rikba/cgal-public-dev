#include "support_plane.h"
#include "universe.h"
#include "parameters.h"
#include "octree_base_vertex.h"

namespace JPTD {

bool Support_Plane::K_based_stopping_condition(bool exist_segments, Polygon_Vertex* v) const
{
	// In this function, we determine if vertex v should keep propagating.

	// If v = nullptr, then it means that this function is called in a context of hole prevention,
	// and therefore only the existence of a segment at v->pt(t) matters.
	// If not, then we check if v->K == 0 and decrement this variable.

	// We return true if the propagation continues, false otherwise.

	if (v == nullptr) {
		return !exist_segments;
	} else {
		if (v->K == 0) {
			return false;
		} else {
			--v->K;
			return true;
		}
	}
}



bool Support_Plane::density_based_stopping_condition(const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t, const CGAL_Vector_2 & dM) const
{
	// At time t, a vertex v intersects a line I in V_t with a propagation speed dM.
	// We would like to define a rectangular cuboid and count the number of input points
	// that belong to this cuboid. Their faces should be expressed as planes.

	// Part 1.
	// We are interested in defining a local 3D frame in W_t whose two vectors
	// correspond to dM and the orthogonal vector to this support plane.

	const CGAL_Point_3 & O = W_t;
	CGAL_Vector_3 j = backproject(V_t + dM) - W_t;
	CGAL_Vector_3 k = plane.orthogonal_vector();
	CGAL_Vector_3 i = CGAL::cross_product(j, k);

	const FT & a = Universe::params->density_box_width;
	const FT & b = Universe::params->density_box_length;
	const FT & c = Universe::params->density_box_height;

	// Part 2.
	// In the local 3D frame (O ni nj nk) where ni, nj and nk are the normalized versions of i, j and k,
	// we consider the 8 points with the following coordinates :
	// P[0 .. 3] : A- = (-a/2, 0, -c/2), D- = (a/2, 0, -c/2), B- = (-a/2, b, -c/2), C- = (a/2, b, -c/2),
	// P[4 .. 7] : A+ = (-a/2, 0,  c/2), D+ = (a/2, 0,  c/2), B+ = (-a/2, b,  c/2), C+ = (a/2, b,  c/2).

	// As we work with i, j and k which are not normalized,
	// we should find approximate coordinates of these 8 points with the exact kernel.

	// We want |t_a| such that O + t_a * i is at distance 'a' from O and belongs to a plane parallel to (O j k),
	//         |t_b| such that O + t_b * j is at distance 'b' from O and belongs to a plane parallel to (O i k),
	//         |t_c| such that O + t_c * k is at distance 'c' from O and belongs to a plane parallel to (O i j).

	FT half_t_a = project_to_parallel_plane(a, i) / FT(2);
	FT t_b = project_to_parallel_plane(b, j);
	FT half_t_c = project_to_parallel_plane(c, k) / FT(2);

	std::vector<CGAL_Point_3> P;
	for (int r = 0; r < 8; r++) {
		FT x = (r % 2 == 0 ? -half_t_a : half_t_a);
		FT y = ((r / 2) % 2 == 0 ? 0 : t_b);
		FT z = (r > 3 ? -half_t_c : half_t_c);
		P.push_back(O + x * i + y * j + z * k);
	}

	// Step 3.
	// Reasons on the minimal/maximal coordinates of the 8 points along all axes
	// in order to formulate a query for the octree

	double _query[6] = { FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX };
	std::vector<double> query(_query, _query + sizeof(_query) / sizeof(double));

	for (int r = 0; r < 8; r++) {
		double x = CGAL::to_double(P[r].x()), y = CGAL::to_double(P[r].y()), z = CGAL::to_double(P[r].z());
		if (x < query[0]) query[0] = x;
		if (x > query[1]) query[1] = x;
		if (y < query[2]) query[2] = y;
		if (y > query[3]) query[3] = y;
		if (z < query[4]) query[4] = z;
		if (z > query[5]) query[5] = z;
	}

	std::list< Octree_Base_Vertex<CGAL_Point_3, CGAL_Inexact_Point_3>* > C;
	Universe::point_cloud_octree->search(query, C);

	if (int(C.size()) < Universe::params->density_pts) {
		return false;
	}

	// Step 4.
	// Filters the results returned by the previous query,
	// by discarding the points that are not inside the six planes that delimit the previous cuboid.
	// First we set the plane equations. For each plane the three first coefficients are provided by [ijk].xyz(),
	// the last one can be deduced from a point included in the plane.

	FT I_inf = -(i.x() * P[0].x() + i.y() * P[0].y(), i.z() * P[0].z()), I_sup = -(i.x() * P[7].x() + i.y() * P[7].y(), i.z() * P[7].z());
	FT J_inf = -(j.x() * P[0].x() + j.y() * P[0].y(), j.z() * P[0].z()), J_sup = -(j.x() * P[7].x() + j.y() * P[7].y(), j.z() * P[7].z());
	FT K_inf = -(k.x() * P[0].x() + k.y() * P[0].y(), k.z() * P[0].z()), K_sup = -(k.x() * P[7].x() + k.y() * P[7].y(), k.z() * P[7].z());

	// Now loops on the points.
	// Initially, we consider that all points are between each couple of plane (I_inf, I_sup), (J_inf, J_sup), (K_inf, K_sup).
	// We decrease the number of points if it appears that one of these relations is not satisfied.

	int inside = C.size();

	for (auto it_p = C.begin(); it_p != C.end(); ++it_p) {
		if (inside < Universe::params->density_pts) return false;

		const CGAL_Point_3 & M = (*it_p)->M;
		FT d, d_inf, d_sup;

		// First test : I_inf, I_sup
		d = i.x() * M.x() + i.y() * M.y() + i.z() * M.z();
		d_inf = d + I_inf, d_sup = d + I_sup;
		if (d_inf * d_sup > 0) {
			--inside;
			continue;
		}

		// Second test : J_inf, J_sup
		d = j.x() * M.x() + j.y() * M.y() + j.z() * M.z();
		d_inf = d + J_inf, d_sup = d + J_sup;
		if (d_inf * d_sup > 0) {
			--inside;
			continue;
		}

		d = k.x() * M.x() + k.y() * M.y() + k.z() * M.z();
		d_inf = d + K_inf, d_sup = d + K_sup;
		if (d_inf * d_sup > 0) {
			--inside;
			continue;
		}
	}

	return (inside >= Universe::params->density_pts);
}



FT Support_Plane::project_to_parallel_plane(const FT & D, const CGAL_Vector_3 & n) const
{
	// For a point M (x, y, z) in a plane (P) ax + by + cz + d = 0,
	// We want N in a plane (P') ax + by + cz + d' = 0 at a distance D from M.
	// Such point is obtained by computing N = M + t * n (= vec(a, b, c)) where t = D / sqrt(a2 + b2 + c2).

	FT s = n.x() * n.x() + n.y() * n.y() + n.z() * n.z();

	double s_d = CGAL::to_double(s);
	double sqrt_s_d = sqrt(s_d);

	std::stringstream stream;
	stream << std::setprecision(15);
	stream << sqrt_s_d;

	std::string str_sqrt_s = stream.str();
	FT sqrt_s = FT(str_sqrt_s);

	return D / sqrt_s;
}



bool Support_Plane::propagation_continues_outside_intersections(Intersection_Line* I, Polygon_Vertex* v, const CGAL_Point_2 & V_t, const CGAL_Point_3 & W_t, const FT & t, const std::vector<bool> & S, const int seed) const
{
	// The propagation of a vertex v, such that V_t = v(t) at time t, continues beyond line I
	// when the three following conditions are respected :

	// a) there is no polygon already existing in the cell identified by S ;
	// b) I is not a border.
	// c) there is no couple of segments, of opposite signs, supported by I containing V_t ;

	/*bool keep_propagating = !(polygon_set->exists(S, seed) || I->is_border || I->exist_segments_including_point_outside_intersections(V_t, t));
	return keep_propagating;*/

	if (polygon_set->exists(S, seed) || I->is_border) {
		return false;
	}

	bool exist_segments = I->exist_segments_including_point_outside_intersections(V_t, t);
	if (!exist_segments) return true;

	if (Universe::params->stopping_condition == 0) {
		return K_based_stopping_condition(exist_segments, v);
	} else {
		return density_based_stopping_condition(V_t, W_t, v->dM);
	}
}



bool Support_Plane::propagation_continues_at_intersection(Intersection_Line* I, Polygon_Vertex* v, const CGAL_Point_2 & V_t, 
	const CGAL_Point_3 & W_t, const Constraint & C, const FT & t, const std::vector<bool> & S, const int seed) const
{
	// Same function as above, except that we make use of the fact that V_t is the intersection of I and C.first

	/*bool keep_propagating = !(polygon_set->exists(S, seed) || I->is_border || I->exist_segments_including_point_at_intersection(V_t, C, t));
	return keep_propagating;*/

	if (polygon_set->exists(S, seed) || I->is_border) {
		return false;
	}

	bool exist_segments = I->exist_segments_including_point_at_intersection(V_t, C, t);
	if (!exist_segments) return true;

	if (Universe::params->stopping_condition == 0) {
		return K_based_stopping_condition(exist_segments, v);
	} else {
		return density_based_stopping_condition(V_t, W_t, v->dM);
	}
}



bool Support_Plane::propagation_continues_at_intersection(Intersection_Line* I, Polygon_Vertex* v, const CGAL_Point_2 & V_t, 
	const CGAL_Point_3 & W_t, const std::list<Constraint> & C_limits, const FT & t, const std::vector<bool> & S, const int seed) const
{
	// Same function as above, in the case of a multi-line intersection

	/*bool keep_propagating = !(polygon_set->exists(S, seed) || I->is_border || I->exist_segments_including_point_at_intersection(V_t, C_limits, t));
	return keep_propagating;*/

	if (polygon_set->exists(S, seed) || I->is_border) {
		return false;
	}

	bool exist_segments = I->exist_segments_including_point_at_intersection(V_t, C_limits, t);
	if (!exist_segments) return true;

	if (Universe::params->stopping_condition == 0) {
		return K_based_stopping_condition(exist_segments, v);
	} else {
		return density_based_stopping_condition(V_t, W_t, v->dM);
	}
}

}