#include "propagation.h"
#include "universe.h"
#include "support_plane.h"
#include "event.h"
#include "event_queue.h"
#include "partition.h"
#include "ply_in.h"
#include "ply_out.h"
#include "stats.h"
#include "vars.h"

#include <fstream>
#include <iomanip>
#include <CGAL/intersections.h>
#include <boost/filesystem.hpp>
#include <boost/algorithm/string.hpp>
#include <CGAL/linear_least_squares_fitting_2.h>
#include <CGAL/linear_least_squares_fitting_3.h>
#include <CGAL/convex_hull_2.h>

#include <CGAL/point_generators_2.h>
#include <CGAL/random_convex_set_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_with_holes_2.h>
#include <CGAL/Boolean_set_operations_2.h>
#include <CGAL/Polygon_2_algorithms.h>

namespace JPTD {

using namespace boost::filesystem;
using CGAL::to_double;



Kinetic_Propagation::Kinetic_Propagation(const int N)
{
	generate_polygons_randomly(N);
}



Kinetic_Propagation::Kinetic_Propagation(const std::string & filename, Preprocess process)
{
	try {
		get_polygons_from_input_file(filename);
	} catch (std::exception & e) {
		polygons.clear();
		std::cout << e.what() << std::endl;
	}

	if (process != NONE) {
		apply(process);
		set_option("--input", filename);
	}
}



Kinetic_Propagation::Kinetic_Propagation(int argc, char *argv[], Preprocess process)
{
	set_options(argc, argv);

	try {
		const std::string & filename = Universe::params->path;
		const int N = Universe::params->rand_n;
		const int P = Universe::params->rand_p;
		const double D = Universe::params->rand_d;

		if (!filename.empty()) {
			get_polygons_from_input_file(filename);
		} else {
			generate_polygons_randomly(N, P, D);
		}

		const std::string & filename_point_cloud = Universe::params->path_point_cloud;
		if (!filename_point_cloud.empty()) {
			init_point_cloud_octree(filename_point_cloud);
		}

	} catch (std::exception & e) {
		polygons.clear();
		std::cout << e.what() << std::endl;
	} 

	if (process != NONE) {
		apply(process);
	}
}



Kinetic_Propagation::Kinetic_Propagation(const std::vector<std::vector<CGAL_Inexact_Point_3> > & primitives, Preprocess process)
{
	const size_t n = primitives.size();

	polygons.clear();
	polygons.reserve(n);
	for (size_t i = 0 ; i < n ; i++) {

		std::vector<CGAL_Point_3> P;
		const size_t p = primitives[i].size();
		P.reserve(p);

		for (size_t j = 0 ; j < p ; j++) {
			const CGAL_Inexact_Point_3 & p_ij = primitives[i][j];
			FT x = p_ij.x(), y = p_ij.y(), z = p_ij.z();
			P.push_back(CGAL_Point_3(x, y, z));
		}
		polygons.push_back(P);
	}

	if (process != NONE) {
		apply(process);
	}
}



Kinetic_Propagation::Kinetic_Propagation(const std::vector<std::vector<CGAL_Point_3> > & primitives, Preprocess process)
{
	const size_t n = primitives.size();

	polygons.clear();
	polygons.reserve(n);
	for (size_t i = 0 ; i < n ; i++) {
		polygons.push_back(primitives[i]);
	}

	if (process != NONE) {
		apply(process);
	}
}



void Kinetic_Propagation::set_options(int argc, char *argv[]) const
{
	Universe::params->set_options(argc, argv);
}



void Kinetic_Propagation::set_option(const std::string & option) const
{
	Universe::params->set_option(option);
}



void Kinetic_Propagation::set_option(const std::string & option, const std::string & argument_1) const
{
	Universe::params->set_option(option, argument_1);
}



void Kinetic_Propagation::set_option(const std::string & option, const std::string & argument_1, const std::string & argument_2) const
{
	Universe::params->set_option(option, argument_1, argument_2);
}



Kinetic_Propagation::~Kinetic_Propagation()
{
	if (partition != nullptr) {
		delete partition;
	}
}


bool Kinetic_Propagation::data() const
{
	return (!polygons.empty());
}



void Kinetic_Propagation::run()
{
	if (polygons.empty()) return;

	try {
		// Initializes a kinetic data structure
		init_kinetic_data_structure();

		// Processes all events
		unstack();

		// Builds the partition
		build_partition();

		// Cleans memory
		// delete_kinetic_data_structure();

		/*
		{
			std::cout << "** schedule_events() : " << std::endl;
			std::cout << "calls      : " << KP_Stats::schedule_events_vertices << std::endl;
			std::cout << "lines      : " << KP_Stats::schedule_events_lines << std::endl;
			std::cout << "average    : " << double(KP_Stats::schedule_events_lines) / KP_Stats::schedule_events_vertices << std::endl;
			std::cout << "searching  : " << double(KP_Stats::schedule_events_search_time) / CLOCKS_PER_SEC << " s." << std::endl;
			std::cout << "scheduling : " << double(KP_Stats::schedule_events_computation_time) / CLOCKS_PER_SEC << " s." << std::endl;

			std::cout << "** average vertex life expectancy : " << std::endl;
			std::cout << "lines      : " << double(KP_Stats::life_expectancy_lines) / KP_Stats::life_expectancy_terms << std::endl;
			std::cout << "distance   : " << double(KP_Stats::life_expectancy_distance) / KP_Stats::life_expectancy_terms << std::endl;

			std::cout << "** process_events() : " << std::endl;
			std::cout << "time       : " << double(KP_Stats::process_events_time) / CLOCKS_PER_SEC << " s." << std::endl;
			std::cout << "events     : " << KP_Stats::process_events_calls << std::endl;
			std::cout << "average    : " << double(KP_Stats::process_events_time) / (CLOCKS_PER_SEC * KP_Stats::process_events_calls) << std::endl;
		} */

	} catch (std::exception & except) {
		std::cout << except.what() << std::endl;
		throw except;
	}
}



void Kinetic_Propagation::get_polygons_from_input_file(const std::string & filename)
{
	std::string extension = path(filename).extension().string();

	try {
		if (extension == ".txt") {
			get_polygons_from_txt_file(filename);
		} else if (extension == ".ply") {
			get_polygons_from_ply_file(filename);
		} else {
			throw std::ios_base::failure("Error : unexpected input format.");
		}
	} catch (std::exception & except) {
		throw except;
	}
}



void Kinetic_Propagation::get_polygons_from_txt_file(const std::string & filename)
{
	// Opens file
	std::ifstream file (filename);
	if (!file.is_open()) {
		throw std::ios_base::failure("Error : input file cannot be opened.");
	}

	try {
		std::string line;
		if (!std::getline(file, line)) {
			throw std::logic_error("Parsing error in TXT file : unexpected end of file.");
		}

		int n, p;
		std::istringstream stream (line);

		// The first element to read is the number of polygons
		if (!(stream >> n)) {
			throw std::logic_error("Parsing error in TXT file : a number of polygons was expected.");
		}
		polygons.reserve(n);

		// Once the number of polygons is known, we start enumerating them
		// A polygon consists in a number of vertices, then a list of vertices

		for (int i = 0 ; i < n ; i++) {

			// Gets a number of vertices
			if (!std::getline(file, line)) {
				throw std::logic_error("Parsing error in TXT file : unexpected end of file.");
			}
			stream = std::istringstream (line);
			if (!(stream >> p)) {
				throw std::logic_error("Parsing error in TXT file : a number of points was expected.");
			}

			// Reads points
			std::vector<CGAL_Point_3> polygon;
			polygon.reserve(p);
			for (int j = 0 ; j < p ; j++) {

				if (!std::getline(file, line)) {
					throw std::logic_error("Parsing error in TXT file : unexpected end of file.");
				}

				CGAL_Point_3 pt;
				stream = std::istringstream (line);
				if (!(stream >> pt)) {
					throw std::logic_error("Error : a triplet of coordinates was expected.");
				}
				polygon.push_back(pt);
			}

			// Inserts polygon
			polygons.push_back(polygon);
		}

		file.close();

	} catch (std::exception & e) {
		polygons.clear();
		file.close();
		throw e;
	}
}



void Kinetic_Propagation::get_polygons_from_ply_file(const std::string & filename)
{
	try {
		Ply_In::read(filename, polygons);
	} catch (std::exception & except) {
		throw except;
	}
}



void Kinetic_Propagation::generate_polygons_randomly(const int N, const int P, const double D)
{
	std::string basename = "rand_" + std::to_string(N) + "_" + std::to_string(P);
	set_option("--basename", basename);

	std::default_random_engine & generator = Universe::generator;
	std::uniform_real_distribution<double> R (-1, 1);

	typedef CGAL::Random_points_in_square_2<CGAL_Point_2, CGAL::Creator_uniform_2<double, CGAL_Point_2> > Point_generator;

	// Part 1.
	// Although this is not an optimal strategy, given the fact that it will be computed again later,
	// we have to compute the bounding polygon of a random plane inside the 3D cube [-1, 1].
	// This operation has to be repeated N times (N is the number of polygons to generate).
	// Here, we define several variables that will be used as we compute that bounding polygon.

	CGAL_Point_3 pt_min (-1, -1, -1);
	CGAL_Point_3 pt_max (1, 1, 1);
	std::vector<CGAL_Point_3> box_corners;
	std::vector<std::pair<int, int> > box_edges;

	int vertices[8][3] = { {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1} };
	int edges[12][2] = { {0, 2}, {1, 3}, {4, 6}, {5, 7}, {0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 4}, {2, 6}, {1, 5}, {3, 7} };

	for (int i = 0 ; i < 8 ; i++) {
		FT x = (vertices[i][0] == 0 ? pt_min.x() : pt_max.x());
		FT y = (vertices[i][1] == 0 ? pt_min.y() : pt_max.y());
		FT z = (vertices[i][2] == 0 ? pt_min.z() : pt_max.z());
		box_corners.push_back(CGAL_Point_3(x, y, z));
	}

	for (int i = 0 ; i < 12 ; i++) {
		box_edges.push_back(std::make_pair(edges[i][0], edges[i][1]));
	}

	// Part 2.
	// Generates N random polygons

	polygons.reserve(N);

	while (polygons.size() < N) {

		// Part 2.1
		// Generates a random plane equation that intersects the cube
		// by picking a random point inside the cube and a random normal

		const double x_0 = R(generator), y_0 = R(generator), z_0 = R(generator);
		const CGAL_Point_3 O_ref(x_0, y_0, z_0);

		const double a = R(generator), b = R(generator), c = R(generator);
		const double d = -(a * x_0 + b * y_0 + c * z_0);
		const CGAL_Plane P_ref (a, b, c, d);

		// Part 2.2
		// Generates the 3D intersection of the plane with the bounding box
		// Maps the points to the plane to obtain a 2D polygon

		std::list<CGAL_Point_3> BP_ref_3d;
		std::vector<CGAL_Point_2> BP_ref_2d;

		Support_Plane::construct_bounding_polygon_of_support_plane(pt_min, pt_max, box_corners, box_edges, P_ref, BP_ref_3d);

		BP_ref_2d.reserve(BP_ref_3d.size());
		for (std::list<CGAL_Point_3>::iterator it_p = BP_ref_3d.begin() ; it_p != BP_ref_3d.end() ; it_p++) {
			BP_ref_2d.push_back(P_ref.to_2d(*it_p));
		}

		CGAL_Polygon_2 BP_ref (BP_ref_2d.begin(), BP_ref_2d.end());
		if (BP_ref.orientation() == CGAL::Orientation::CLOCKWISE) {
			BP_ref.reverse_orientation();
		}

		// Part 2.3
		// Generates a convex hull of P points

		std::vector<CGAL_Point_2> CH;
		CGAL::random_convex_set_2(P, std::back_inserter(CH), Point_generator(D));
		
		std::vector<CGAL_Point_2> K_pts;
		for (size_t j = 0 ; j < CH.size() ; j++) {
			FT lambda = CH[j].x(), mu = CH[j].y();
			CGAL_Point_3 M = O_ref + lambda * P_ref.base1() + mu * P_ref.base2();
			K_pts.push_back(P_ref.to_2d(M));
		}

		CGAL_Polygon_2 K (K_pts.begin(), K_pts.end());

		// Part 2.4
		// Computes the intersection of the plane and the convex hull
		// This will discards points that are not inside the cube [-1, 1]
		// Finally backprojects the points in 3D

		std::list<CGAL_Polygon_with_holes_2> BP_K_intersection;
		CGAL::intersection(BP_ref, K, std::back_inserter(BP_K_intersection));

		if (!BP_K_intersection.empty()) {

			// This intersection may be empty :
			// Even if the idea is to generate points all around O_ref,
			// it is not garanteed that all points will be inside the bounding cube.

			// However, if the intersection exists,
			// it consists of a unique polygon because K and BP are convex
			// (CGAL imposes to use Polygon_with_holes_2)

			for (std::list<CGAL_Polygon_with_holes_2>::iterator it_p = BP_K_intersection.begin(); it_p != BP_K_intersection.end(); it_p++) {
				CGAL_Polygon_2 S = it_p->outer_boundary();

				std::vector<CGAL_Point_3> P_generated;
				P_generated.reserve(S.size());

				for (std::vector<CGAL_Point_2>::iterator it_v = S.vertices_begin(); it_v != S.vertices_end(); it_v++) {
					P_generated.push_back(P_ref.to_3d(*it_v));
				}

				// Finally inserts the polygon
				polygons.push_back(P_generated);
			}
		}
	}

	// Part 3.
	// Prints the vertices

	std::ofstream stream(basename + ".txt", std::ofstream::out);
	if (stream.is_open()) {
		stream << polygons.size() << std::endl;

		for (int i = 0 ; i < polygons.size() ; i++) {
			stream << polygons[i].size() << std::endl;
		
			for (int j = 0 ; j < int(polygons[i].size()) ; j++) {
				const CGAL_Point_3 & P_ij = polygons[i][j];
				stream << P_ij.x() << " " << P_ij.y() << " " << P_ij.z() << std::endl;
			}
		}

		stream.close();
	}
}



void Kinetic_Propagation::generate_polygons_randomly(const int N)
{
	std::string basename = "rand_" + std::to_string(N);
	set_option("--basename", basename);

	std::default_random_engine & generator = Universe::generator;
	std::uniform_real_distribution<double> R (-1, 1);
	std::exponential_distribution<double> E (0.8);

	for (int i = 0 ; i < N ; i++) {

		// Generates a random 3D vector n
		const double a = R(generator), b = R(generator), c = R(generator);

		// Generates a random 3D point, and gets a equation of plane
		const double x_0 = R(generator), y_0 = R(generator), z_0 = R(generator);
		const CGAL_Point_3 O_ref(x_0, y_0, z_0);

		const double d = -(a * x_0 + b * y_0 + c * z_0);
		const CGAL_Plane P_i(a, b, c, d);

		// Generates a number of polygons
		int M = 1; // int(ceil(E(generator)));

		for (int j = 0; j < M; j++) {

			// Generates 4 couples of coordinates (mu, nu) corresponding 
			// to a linear combination of the base vectors of the plane
			const CGAL_Vector_3 u = P_i.base1(), v = P_i.base2();
			const CGAL_Point_3 O (O_ref + R(generator) * u + R(generator) * v);

			// mu and nu belongs to [0, 0.3]
			const double mu_1 = 0.15 * (R(generator) + 1);
			const double mu_2 = 0.15 * (R(generator) + 1);
			const double nu_1 = 0.15 * (R(generator) + 1);
			const double nu_2 = 0.15 * (R(generator) + 1);

			std::vector<CGAL_Point_3> P;
			P.reserve(4);

			P.push_back(O + mu_1 * u);
			P.push_back(O + nu_1 * v);
			P.push_back(O - mu_2 * u);
			P.push_back(O - nu_2 * v);

			// Inserts P
			polygons.push_back(P);
		}
	}

	// Prints the vertices

	std::ofstream stream(basename + ".txt", std::ofstream::out);
	if (stream.is_open()) {
		stream << polygons.size() << std::endl;

		for (int i = 0 ; i < polygons.size() ; i++) {
			stream << polygons[i].size() << std::endl;
		
			for (int j = 0 ; j < int(polygons[i].size()) ; j++) {
				const CGAL_Point_3 & P_ij = polygons[i][j];
				stream << P_ij.x() << " " << P_ij.y() << " " << P_ij.z() << std::endl;
			}
		}

		stream.close();
	}
}



void Kinetic_Propagation::init_point_cloud_octree(const std::string & filename)
{
	// We suppose that filename points to a list of points.
	// We call our simplified parser to read it.
	std::vector<CGAL_Point_3> C_points;
	Ply_In::read(filename, C_points);
 
	// We convert the CGAL_Point_3 objects into MPIR_Point_3 points.
	// We simultaneously identify the extrema along each axis
	double _extrema[6] = { FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX, FLT_MAX, -FLT_MAX };
	std::vector<double> extrema(_extrema, _extrema + sizeof(_extrema) / sizeof(double));

	std::vector<Octree_Base_Vertex<CGAL_Point_3, CGAL_Inexact_Point_3>*> points;
	points.reserve(C_points.size());
	for (size_t i = 0 ; i < C_points.size() ; i++) {
		const CGAL_Point_3 & C = C_points[i];

		double xd = CGAL::to_double(C.x());
		double yd = CGAL::to_double(C.x());
		double zd = CGAL::to_double(C.x());
		CGAL_Inexact_Point_3 hint_C(xd, yd, zd);

		if (xd < extrema[0]) extrema[0] = xd;
		if (xd > extrema[1]) extrema[1] = xd;
		if (yd < extrema[2]) extrema[2] = yd;
		if (yd > extrema[3]) extrema[3] = yd;
		if (zd < extrema[4]) extrema[4] = zd;
		if (zd > extrema[5]) extrema[5] = zd;

		points.push_back(new Octree_Base_Vertex<CGAL_Point_3, CGAL_Inexact_Point_3>(C, hint_C));
	}

	// Initialize the octree
	Universe::point_cloud_octree = new Octree_Base<CGAL_Point_3, CGAL_Inexact_Point_3>(extrema);
	for (size_t i = 0 ; i < points.size() ; i++) {
		Universe::point_cloud_octree->add(points[i]);
	}
}



void Kinetic_Propagation::apply(Preprocess process)
{
	if (process == CONVEX_ENVELOPS) {
		make_convex_envelops();
	} else if (process == OPTIMAL_RECTANGLES) {
		make_fit_rectangles();
	}
}



void Kinetic_Propagation::make_convex_envelops()
{
	typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 Inexact_Point_3;
	typedef CGAL::Exact_predicates_inexact_constructions_kernel::Plane_3 Inexact_Plane_3;

	for (size_t i = 0 ; i < polygons.size() ; i++) {

		// Switches from points with exact precision,
		// to points with double coordinates
		std::vector<Inexact_Point_3> Q;
		for (size_t j = 0 ; j < polygons[i].size() ; j++) {
			const CGAL_Point_3 & p_ij = polygons[i][j];
			FT x = p_ij.x(), y = p_ij.y(), z = p_ij.z(); 
			Q.push_back(Inexact_Point_3(to_double(x), to_double(y), to_double(z)));
		}

		// Fits a plane
		Inexact_Plane_3 Q_fit;
		linear_least_squares_fitting_3(Q.begin(), Q.end(), Q_fit, CGAL::Dimension_tag<0>());

		// Maps points of polygons[i] as 2D coordinates

		CGAL_Plane P_fit (Q_fit.a(), Q_fit.b(), Q_fit.c(), Q_fit.d());
		std::vector<CGAL_Point_2> map_of_points, convex_map;

		CGAL_Point_3 O = P_fit.point();
		CGAL_Vector_3 b_1 = P_fit.base1(), b_2 = P_fit.base2();
		FT det = (b_1.x() * b_2.y() - b_1.y() * b_2.x());
		
		for (size_t j = 0 ; j < polygons[i].size() ; j++) {
			CGAL_Point_3 M = P_fit.projection(polygons[i][j]);
			CGAL_Vector_3 OM = M - O;

			// p_ij = O + lambda * u + mu * v : determines lambda and mu
			FT lambda = (OM.x() * b_2.y() - OM.y() * b_2.x()) / det;
			FT mu = (b_1.x() * OM.y() - b_1.y() * OM.x()) / det;
			map_of_points.push_back(CGAL_Point_2(lambda, mu));
		}

		// Computes the convex envelop
		convex_hull_2(map_of_points.begin(), map_of_points.end(), std::back_inserter(convex_map));

		// Backprojects it and replaces polygons[i]
		std::vector<CGAL_Point_3> R;
		R.reserve(convex_map.size());
		for (size_t j = 0 ; j < convex_map.size() ; j++) {
			R.push_back(O + convex_map[i].x() * b_1 + convex_map[i].y() * b_2);
		}

		polygons[i].clear();
		std::copy(R.begin(), R.end(), std::back_inserter(polygons[i]));
	}
}



void Kinetic_Propagation::make_fit_rectangles()
{
	typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_2 Inexact_Point_2;
	typedef CGAL::Exact_predicates_inexact_constructions_kernel::Line_2 Inexact_Line_2;

	typedef CGAL::Exact_predicates_inexact_constructions_kernel::Point_3 Inexact_Point_3;
	typedef CGAL::Exact_predicates_inexact_constructions_kernel::Plane_3 Inexact_Plane_3;

	for (size_t i = 0 ; i < polygons.size() ; i++) {

		// Switches from points with exact precision,
		// to points with double coordinates
		std::vector<Inexact_Point_3> Q;
		for (size_t j = 0 ; j < polygons[i].size() ; j++) {
			const CGAL_Point_3 & p_ij = polygons[i][j];
			FT x = p_ij.x(), y = p_ij.y(), z = p_ij.z(); 
			Q.push_back(Inexact_Point_3(to_double(x), to_double(y), to_double(z)));
		}

		// Fits a plane
		Inexact_Plane_3 Q_fit;
		linear_least_squares_fitting_3(Q.begin(), Q.end(), Q_fit, CGAL::Dimension_tag<0>());

		// Maps points of polygons[i] as 2D coordinates

		CGAL_Plane P_fit (Q_fit.a(), Q_fit.b(), Q_fit.c(), Q_fit.d());
		std::vector<Inexact_Point_2> map_of_points;

		CGAL_Point_3 O = P_fit.point();
		CGAL_Vector_3 b_1 = P_fit.base1(), b_2 = P_fit.base2();
		FT det = (b_1.x() * b_2.y() - b_1.y() * b_2.x());
		
		for (size_t j = 0 ; j < polygons[i].size() ; j++) {
			CGAL_Point_3 M = P_fit.projection(polygons[i][j]);
			CGAL_Vector_3 OM = M - O;

			// p_ij = O + lambda * u + mu * v : determines lambda and mu
			FT lambda = (OM.x() * b_2.y() - OM.y() * b_2.x()) / det;
			FT mu = (b_1.x() * OM.y() - b_1.y() * OM.x()) / det;

			map_of_points.push_back(Inexact_Point_2(to_double(lambda), to_double(mu)));
		}

		// Fits a line to the map of points
		Inexact_Line_2 L_fit;
		linear_least_squares_fitting_2(map_of_points.begin(), map_of_points.end(), L_fit, CGAL::Dimension_tag<0>());

		// Computes the coordinates of the 
		double lambda_min = FLT_MAX, lambda_max = -lambda_min;
		double mu_min = FLT_MAX, mu_max = -mu_min;
		for (size_t j = 0 ; j < map_of_points.size() ; j++) {
			double l = map_of_points[i].x(), m = map_of_points[i].y();
			if (l < lambda_min) lambda_min = l;
			if (l > lambda_max) lambda_max = l;
			if (m < mu_min) mu_min = m;
			if (m > mu_max) mu_max = m;
		}

		// We have obtained the coordinates of the four vertices of the best fit rectangle
		// We backproject it and replace polygons[i]
		std::vector<CGAL_Point_3> R;
		R.push_back(O + FT(lambda_min) * b_1 + FT(mu_min) * b_2);
		R.push_back(O + FT(lambda_min) * b_1 + FT(mu_max) * b_2);
		R.push_back(O + FT(lambda_max) * b_1 + FT(mu_max) * b_2);
		R.push_back(O + FT(lambda_max) * b_1 + FT(mu_min) * b_2);

		polygons[i].clear();
		std::copy(R.begin(), R.end(), std::back_inserter(polygons[i]));
	}
}



void Kinetic_Propagation::init_kinetic_data_structure()
{
	CGAL_Point_3 pt_min, pt_max;
	std::vector<CGAL_Point_3> box_corners;
	std::vector<std::pair<int, int> > box_edges;
	std::map<int, std::map<int, CGAL_Line_3> > lines;

	try {
		// Initialization, part 1 : geometric data structures
		init_bounding_box(pt_min, pt_max, box_corners, box_edges);
		init_plane_equations(90);
		init_intersection_lines(lines);

		init_supporting_planes(pt_min, pt_max, box_corners, box_edges, lines);
		build_polygons();

		// Initialization, part 2 : optimized structure of queue
		init_schedule();

	} catch (std::exception & except) {
		throw except;
	}

	// std::cout << "** Initialized structure" << std::endl;
}



void Kinetic_Propagation::delete_kinetic_data_structure()
{
	// Deletes the support planes and their objects
	for (size_t i = 0 ; i < Universe::map_of_planes.size() ; i++) {
		delete Universe::map_of_planes[i];
	}
	Universe::map_of_planes.clear();

	// Deletes the queue
	delete Universe::event_queue;

	// Delete counters
	Counters::id_planes = -1;
	Counters::id_objects = -1;

	Counters::id_partition_vertex = -1;
	Counters::id_partition_edge = -1;
	Counters::id_partition_facet = -1;
	Counters::id_partition_polyhedron = -1;

	Counters::par_v_local_ids = std::vector<int>();
	Counters::par_e_local_ids = std::vector<int>();
}



void Kinetic_Propagation::init_bounding_box(CGAL_Point_3 & pt_min,
	CGAL_Point_3 & pt_max,
	std::vector<CGAL_Point_3> & box_corners,
	std::vector<std::pair<int, int> > & box_edges)
{
	// Step 1.
	// First of all, we set the dimensions of the bounding box, by reading the coordinates of the different polygons.

	FT x_min = FLT_MAX, x_max = -FLT_MAX;
	FT y_min = FLT_MAX, y_max = -FLT_MAX;
	FT z_min = FLT_MAX, z_max = -FLT_MAX;

	for (int i = 0 ; i < polygons.size() ; i++) {
		for (int j = 0 ; j < polygons[i].size(); j++) {

			const CGAL_Point_3 & p_ij = polygons[i][j];
			const FT x = p_ij.x(), y = p_ij.y(), z = p_ij.z();

			// Updates extrema
			if (x < x_min) x_min = x;
			if (x > x_max) x_max = x;
			if (y < y_min) y_min = y;
			if (y > y_max) y_max = y;
			if (z < z_min) z_min = z;
			if (z > z_max) z_max = z;
		}
	}

	// Extra margins : x_min and x_max are shifted by 5% of (x_max - x_min),
	// and the same operation is applied to the other dimensions of the bounding box

	FT m_x = (x_max - x_min) / FT(10);
	FT m_y = (y_max - y_min) / FT(10);
	FT m_z = (z_max - z_min) / FT(10);

	x_min -= m_x, y_min -= m_y, z_min -= m_z;
	x_max += m_x, y_max += m_y, z_max += m_z;

	pt_min = CGAL_Point_3(x_min, y_min, z_min);
	pt_max = CGAL_Point_3(x_max, y_max, z_max);

	// Step 2.
	// Defines the corners of the bounding box.

	int vertices[8][3] = { {0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {1, 1, 0}, {0, 0, 1}, {1, 0, 1}, {0, 1, 1}, {1, 1, 1} };
	int edges[12][2] = { {0, 2}, {1, 3}, {4, 6}, {5, 7}, {0, 1}, {2, 3}, {4, 5}, {6, 7}, {0, 4}, {2, 6}, {1, 5}, {3, 7} };

	for (int i = 0 ; i < 8 ; i++) {
		FT x = (vertices[i][0] == 0 ? pt_min.x() : pt_max.x());
		FT y = (vertices[i][1] == 0 ? pt_min.y() : pt_max.y());
		FT z = (vertices[i][2] == 0 ? pt_min.z() : pt_max.z());
		box_corners.push_back(CGAL_Point_3(x, y, z));
	}

	for (int i = 0 ; i < 12 ; i++) {
		box_edges.push_back(std::make_pair(edges[i][0], edges[i][1]));
	}

	// Step 3.
	// We define the six first planes of the universe, that correspond to the six facets of the bounding box
	// We also compute the centers of the facet

	planes.push_back(CGAL_Plane(1, 0, 0, -x_min));
	planes.push_back(CGAL_Plane(1, 0, 0, -x_max));
	planes.push_back(CGAL_Plane(0, 1, 0, -y_min));
	planes.push_back(CGAL_Plane(0, 1, 0, -y_max));
	planes.push_back(CGAL_Plane(0, 0, 1, -z_min));
	planes.push_back(CGAL_Plane(0, 0, 1, -z_max));
}



void Kinetic_Propagation::init_plane_equations(const int subdivs)
{
	std::map<std::pair<int, int>, std::list<int> > atlas;

	for (int i = 0 ; i < polygons.size() ; i++) {

		std::vector<CGAL_Point_3> & polygon = polygons[i];
		int n = int(polygon.size());

		// We compute the barycenter of the polygon.

		FT bar_x = 0, bar_y = 0, bar_z = 0;
		for (int j = 0; j < n; j++) {
			const CGAL_Point_3 & p_ij = polygon[j];
			bar_x += p_ij.x();
			bar_y += p_ij.y();
			bar_z += p_ij.z();

		}

		bar_x /= n, bar_y /= n, bar_z /= n;
		CGAL_Point_3 O(bar_x, bar_y, bar_z);

		// We estimate the normal to the polygon.
		// We discretize that normal as a couple (longitude, latitude).

		CGAL_Vector_3 N = CGAL::NULL_VECTOR;

		for (int j = 0 ; j < n / 3 ; j++) {
			CGAL_Point_3 P_0 = polygon[j], P_1 = polygon[n / 3 + j], P_2 = polygon[2 * n / 3 + j];
			N = CGAL::cross_product(P_1 - P_0, P_2 - P_0);
			if (N != CGAL::NULL_VECTOR) {
				if (N.z() < 0 || (N.z() == 0 && N.y() < 0) || (N.z() == 0 && N.y() == 0 && N.x() < 0)) N = -N;
				break;
			}
		}
		
		if (N == CGAL::NULL_VECTOR) {
			for (int j = 0 ; j < polygon.size() ; j++) {
				std::cout << polygon[j] << std::endl;
			}
			throw std::logic_error("Error : read polygon is degenerate.");
		}

		int k_longitude, k_latitude;
		discretize_normal(subdivs, N, k_longitude, k_latitude);
		
		FT a = 0, b = 0, c = 0;
		bool compute_plane_coefficients = true;
		int plane_index = -1;

		// Now, two possible solutions :
		// either there already exists a parallel plane sufficiently close to O in the list of planes,
		// and if so we consider that the current polygon actually belongs to that plane,
		// or such a plane doesn't exist and we register P.

		std::pair<int, int> key = std::make_pair(k_longitude, k_latitude);
		std::map<std::pair<int, int>, std::list<int> >::iterator it_key = atlas.find(key);

		if (it_key != atlas.end()) {
			const std::list<int> & parallel_planes = it_key->second;
			const CGAL_Plane & P_parallel = planes[parallel_planes.front()];

			a = P_parallel.a(), b = P_parallel.b(), c = P_parallel.c();
			compute_plane_coefficients = false;

			FT dist_min = 1e-4;
			for (std::list<int>::const_iterator it_p = parallel_planes.begin() ; it_p != parallel_planes.end() ; it_p++) {
				const CGAL_Plane & P = planes[*it_p];
				CGAL_Vector_3 dP(P.projection(O) - O);

				if (dP.squared_length() < dist_min) {
					dist_min = dP.squared_length();
					plane_index = (*it_p);
				}
			}
		}

		if (plane_index != -1) {
			// Assimilates the polygon to the plane pointed by the previous variable
			polygons_to_planes.push_back(plane_index);

		} else {
			// Registers a new plane
			if (compute_plane_coefficients) {

				if (k_latitude == 0) {
					a = 0, b = 0, c = 1;
				} else if (k_latitude == 90) {
					if (k_longitude == 0 || k_longitude == 360) {
						a = 1;
						b = 0;
					} else if (k_longitude == 90) {
						a = 0;
						b = 1;
					} else if (k_longitude == 180) {
						a = -1;
						b = 0;
					} else if (k_longitude == 270) {
						a = 0;
						b = -1;
					} else {
						a = cos(k_longitude * PI / 180);
						b = sin(k_longitude * PI / 180);
					}
					c = 0;
				} else if (k_longitude % 90 == 0) {
					if (k_longitude == 0 || k_longitude == 360) {
						a = sin(k_latitude * PI / 180);
						b = 0;
					} else if (k_longitude == 90) {
						a = 0;
						b = sin(k_latitude * PI / 180);
					} else if (k_longitude == 180) {
						a = -sin(k_latitude * PI / 180);
						b = 0;
					} else if (k_longitude == 270) {
						a = 0;
						b = -sin(k_latitude * PI / 180);
					}
					c = cos(k_latitude * PI / 180);
				} else {
					a = sin(k_latitude * PI / 180) * cos(k_longitude * PI / 180);
					b = sin(k_latitude * PI / 180) * sin(k_longitude * PI / 180);
					c = cos(k_latitude * PI / 180);
				}
				
				if (c < 0 || (c == 0 && b < 0) || (c == 0 && b == 0 && a < 0)) {
					if (a != 0) a = -a;
					if (b != 0) b = -b;
					if (c != 0) c = -c;
				}
			}

			FT d = -(a * O.x() + b * O.y() + c * O.z());
			CGAL_Plane P (a, b, c, d);
			plane_index = int(planes.size());

			planes.push_back(P);
			polygons_to_planes.push_back(plane_index);
			atlas[key].push_back(plane_index);
		}
 		
		// Projects points on their assigned plane
		for (int j = 0; j < n; j++) {
			polygon[j] = planes[plane_index].projection(polygon[j]);
		}
	}
}
	


void Kinetic_Propagation::discretize_normal(const int n, const CGAL_Vector_3 & v, int & k_longitude, int & k_latitude) const
{
	// We receive a unit vector v.
	// We assume that it is oriented towards the top of the unit sphere.

	// We compute the longitude and latitude values associated to n.
	// These values belong to the [0, 360] and [0, 90] intervals,
	// which are discretized in (4 * n) and (n + 1) values.

	double f_latitude, f_longitude;

	if (v.x() == 0 && v.y() == 0) {
		k_latitude = k_longitude = 0;
		return;
	} else {
		double x = to_double(v.x()), y = to_double(v.y()), z = to_double(v.z());
		double r = sqrt(x * x + y * y + z * z);

		f_latitude = acos(z / r);
		f_longitude = atan2(y / (r * sin(f_latitude)), x / (r * sin(f_latitude)));
		if (f_longitude < 0) f_longitude += 2 * PI;
	}

	f_latitude = f_latitude * 180 / PI;   // in [0, 90]
	f_longitude = f_longitude * 180 / PI; // in [0, 360]

	// Discretizes
	
	int i_inf = f_latitude * n / 90, i_sup = i_inf + 1;
	k_latitude = ((90 * i_sup / n - f_latitude) < (f_latitude - 90 * i_inf / n) ? i_sup : i_inf);

	int j_inf = 4 * f_longitude * n / 360, j_sup = j_inf + 1;
	k_longitude = ((90 * j_sup / n - f_longitude) < (f_longitude - 90 * j_inf / n) ? j_sup : j_inf);
	if (k_longitude == 4 * n) k_longitude = 0;
}



void Kinetic_Propagation::discretize_normal(const int n, const CGAL_Vector_3 & v, CGAL_Vector_3 & v_disc) const
{
	// Same as before but compute the new normal v_disc

	int i_latitude, i_longitude;
	discretize_normal(n, v, i_longitude, i_latitude);

	double latitude = 90 * i_latitude / n;
	double longitude = 90 * i_longitude / n;

	double x = sin(i_latitude * PI / 180) * cos(i_longitude * PI / 180);
	double y = sin(i_latitude * PI / 180) * sin(i_longitude * PI / 180);
	double z = cos(i_latitude * PI / 180);

	v_disc = CGAL_Vector_3(x, y, z);
}



void Kinetic_Propagation::init_intersection_lines(std::map<int, std::map<int, CGAL_Line_3> > & lines) const
{
	int n = int(planes.size());

	// We compute the intersection of each couple of planes (P_i, P_j).
	// If this intersection exists, and is a line, then it is added to the map of intersection lines.

	for (int i = 0 ; i < n ; i++) {
		const CGAL_Plane & P_i = planes[i];
		const CGAL_Vector_3 N_i = P_i.orthogonal_vector();

		for (int j = i + 1 ; j < n ; j++) {
			const CGAL_Plane & P_j = planes[j];
			const CGAL_Vector_3 N_j = P_j.orthogonal_vector();

			CGAL::cpp11::result_of<K::Intersect_3(CGAL_Plane, CGAL_Plane)>::type object = intersection(P_i, P_j);
			if (object) {
				if (const CGAL_Line_3* ptr_L_ij = boost::get<CGAL_Line_3>(&*object)) {
					CGAL_Line_3 L_ij = (*ptr_L_ij);
					lines[i][j] = L_ij;
					lines[j][i] = L_ij;
				}
			}
		}
	}
}



void Kinetic_Propagation::init_supporting_planes(const CGAL_Point_3 & pt_min,
	const CGAL_Point_3 & pt_max,
	const std::vector<CGAL_Point_3> & box_corners,
	const std::vector<std::pair<int, int> > & box_edges,
	const std::map<int, std::map<int, CGAL_Line_3> > & lines) const
{
	int id_vertices[6][4] = { {0, 4, 6, 2}, {1, 5, 7, 3}, {1, 5, 4, 0}, {3, 7, 6, 2}, {1, 0, 2, 3}, {5, 4, 6, 7} };
	int id_facets[6][4] = { {2, 5, 3, 4}, {2, 5, 3, 4}, {1, 5, 0, 4}, {1, 5, 0, 4}, {2, 0, 3, 1}, {2, 0, 3, 1} };

	// There are two different types of support planes to build :
	// - The first category, from planes[0] to planes[5], corresponds to the facets of the bounding box;
	// - The second category, corresponds to the planes for each input planar polygon.

	// This is a two-step initialization process.
	// Over a first phase, all the Support_Plane objects are constructed using their plane equations.
	
	int n = int(planes.size());
	for (int i = 0 ; i < n ; i++) {

		// No need to care about the return value : the constructor adds the address of the new'ed object
		// to a global variable in 'Universe' namespace.
		Support_Plane* SP = new Support_Plane(planes[i]);
	}

	// Once all the Support_Planes exist, and their respective frames are initialized, 
	// we construct the Intersection_Lines and the bounding polygons.
	for (int i = 0 ; i < n ; i++) {

		Support_Plane* SP = Universe::map_of_planes[i];

		// As for the Intersection_Lines, their 3D versions can be obtained very easily
		std::map<int, std::map<int, CGAL_Line_3> >::const_iterator it_lines = lines.find(i);
		const std::map<int, CGAL_Line_3> & L = it_lines->second;
		SP->init_intersection_lines(L);
		SP->init_polygon_set();

		// As for the bounding polygons, there are different ways to obtain them.
		// planes[0 .. 5] correspond to facets of the bounding box.
		// planes[6 ..] are planes associated to polygons, we run a specific algorithm to obtain BP.
		std::list<CGAL_Point_3> BP;
		std::vector<std::list<int> > BF;

		if (i < 6) {
			for (int j = 0 ; j < 4 ; j++) {
				BP.push_back(box_corners[id_vertices[i][j]]);
				BF.push_back(std::list<int>(1, id_facets[i][j]));
			}
		} else {
			const CGAL_Plane & P = planes[i];
			Support_Plane::construct_bounding_polygon_of_support_plane(pt_min, pt_max, box_corners, box_edges, P, BP, BF);
		}
	
		SP->init_bounding_polygon(BP, BF);
		/*if (i >= 6) {
			SP->init_rtree_lines();
		}*/
	}

	// Now that the number of planes is known, initializes the queue of events
	Universe::event_queue = new Event_Queue();
}



void Kinetic_Propagation::build_polygons() const
{
	// Over a first phase, polygons are simply bound to their support plane
	// Segments corresponding to this polygon are created on the other planes
	for (size_t i = 0 ; i < polygons.size() ; i++) {
		
		const std::vector<CGAL_Point_3> & P = polygons[i];
		
		Support_Plane* S = Universe::map_of_planes[polygons_to_planes[i]];
		S->init_polygon(P);
	}
}



void Kinetic_Propagation::init_schedule() const
{
	// Detects events once all vertices are known,
	// and sorts them to initialize the queue

	//clock_t t_begin = clock();

//#pragma omp parallel for
	for (int i = 6 ; i < Universe::map_of_planes.size() ; i++) {
		Universe::map_of_planes[i]->init_schedule();
	}

	/*clock_t t_end = clock();
	std::cout << "Init schedule : " << double(t_end - t_begin) / CLOCKS_PER_SEC << std::endl;
	exit(0);*/
}



void Kinetic_Propagation::unstack()
{
	if (Universe::params->print_drawings) {
		write_polygons(0);
		for (int i = 6; i < Universe::map_of_planes.size(); i++) {
			Universe::map_of_planes[i]->draw(0, 0.5, 3000, 0.5, 0.5, 10);
		}
	}

	Event_Queue* Q = Universe::event_queue;

	while (Universe::moving_objects > 0) {

		// As long as there remain moving objects in the universe, we should expect to pop an event.
		Event_Vertex_Line* e_vl = Q->pop();

		// Finds the Support_Plane to which this event corresponds
		const int intersectant = e_vl->intersectant;
		const FT & t_intersectant = e_vl->t_intersectant;

		Polygon_Vertex* intersectant_object = Universe::map_of_objects[intersectant];
		Support_Plane* SP = Universe::map_of_planes[intersectant_object->id_plane];
		
		if (Universe::params->print_schedule) {
			std::cout << "* [" << SP->id << "] " << e_vl->intersectant << " " << e_vl->intersected << " " << e_vl->t_intersectant << std::endl;
		}

		if (Universe::params->print_drawings && SP->id == 6) {
			SP->draw(t_intersectant, 0.5, 3000, 0.5, 0.5, 10);
		}

		// Requests the plane to process the event
		SP->process_event(e_vl);
	}

	// std::cout << "** Processed all events" << std::endl;

	if (Universe::params->print_drawings) {
		write_polygons(FLT_MAX);
		for (int i = 0; i < Universe::map_of_planes.size(); i++) {
			Universe::map_of_planes[i]->draw(5000, 0.5, 5000, 0.5, 0.5, 10);
		}
	}

	if (Universe::params->check) {
		clock_t t_check_begin = clock();
		for (int i = 0; i < Universe::map_of_planes.size(); i++) {
			Support_Plane* SP = Universe::map_of_planes[i];
			for (auto it_s = SP->segments.begin() ; it_s != SP->segments.end() ; it_s++) {
				it_s->second->check();
			}
		}
		clock_t t_check_end = clock();
		std::cout << "** Checked data in " << double(t_check_end - t_check_begin) / CLOCKS_PER_SEC << " s." << std::endl;
	}
}



#if 0
void Kinetic_Propagation::unstack()
{
	write_polygons(0);
	/*for (int i = 6 ; i < Universe::map_of_planes.size() ; i++) {
		Universe::map_of_planes[i]->draw(0, 0.5, 15000, 0.5, 0.5, 10);
	}*/

	Event_Queue* Q = Universe::event_queue;
	while (Event_Vertex_Line* e_vl = Q->pop()) {

		// Finds the Support_Plane to which this event corresponds
		const int intersectant = e_vl->intersectant;
		const FT t_intersectant = e_vl->t_intersectant;

		Support_Plane_Object* intersectant_object = Universe::map_of_objects[intersectant];
		Support_Plane* SP = Universe::map_of_planes[intersectant_object->id_plane];

		//std::cout << "* [" << SP->id << "] " << e_vl->intersectant << " " << e_vl->intersected << " " << e_vl->t_intersectant << std::endl;
		//SP->draw(t_intersectant, 0.5, 15000, 0.5, 0.5, 10);

		// Requests this plane to process the event
		SP->process_event(e_vl);
		Q->sort();
	}

	//for (int i = 6 ; i < Universe::map_of_planes.size() ; i++) Universe::map_of_planes[i]->draw(10000, 0.5, 5000, 0.5, 0.5, 10);

	write_polygons(FLT_MAX);
	std::cout << "** Processed all events" << std::endl;
}
#endif



void Kinetic_Propagation::build_partition()
{
	const std::string & basename = Universe::params->basename;

	partition = new Partition();
	partition->build();

	// std::cout << "** Built a partition with " << partition->polyhedrons_size() << " polyhedrons" << std::endl;

	if (Universe::params->output_facets) {
		const std::string filename = basename + "_facets.ply";
		partition->ply_facets(filename);
	}

	if (Universe::params->output_polyhedrons) {
		for (std::list<Partition_Polyhedron*>::const_iterator it_p = partition->polyhedrons_begin(); it_p != partition->polyhedrons_end(); it_p++) {
			int i = (*it_p)->id;
			const std::string filename_poly = basename + "_poly_" + std::to_string(i) + ".ply";
			partition->ply_individual_polyhedron(filename_poly, i);
		}
	}
}



void Kinetic_Propagation::write_polygons(const double t) const
{
	const std::string & basename = Universe::params->basename;

	std::string filename = "/Users/danisimo/Documents/pipeline/logs/tmp/bad_input/ply/" + basename;
	if (t == 0) {
		filename += "_in.ply";
	} else if (t == FLT_MAX) {
		filename += "_inf.ply";
	} else {
		std::string stamp = std::to_string(t);
		for (size_t c = 0 ; c < stamp.size() ; c++) {
			if (stamp[c] == '.') stamp[c] = '_'; break;
		}
		filename += ("_" + stamp + ".ply");
	}

	// Gets the definitions of all polygons, once their points no longer propagate
	std::list<std::list<CGAL_Point_3> > polygons;
	std::list<CGAL_Color> colors;
	for (int i = 6 ; i < int(Universe::map_of_planes.size()) ; i++) {
		Universe::map_of_planes[i]->get_polygon_description(polygons, colors, t);
	}

	Ply_Out::print_plain_colorful_facets(filename, polygons, colors, 25);
}

}