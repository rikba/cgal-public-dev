#include "support_plane.h"
#include "svg.h"
#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <CGAL/convex_hull_2.h>

namespace JPTD {
using CGAL::to_double;


void Support_Plane::init_bounding_square(double & x_min, double & y_min, double & x_max, double & y_max, double & dx, double & dy)
{
	/*std::cout << "Plane : " << plane << std::endl;
	for (std::map<int, Planar_Segment *>::iterator it_s = borders.begin() ; it_s != borders.end() ; it_s++) {
		Planar_Segment* s = it_s->second;
		FT x_a = s->A.x(), y_a = s->A.y();
		std::cout << to_double(x_a) << " " << to_double(y_a) << std::endl;
	}*/

	/*if (id == 18) {
		std::cout << "break" << std::endl;
	}*/

	FT ft_x_min = FLT_MAX, ft_x_max = -FLT_MAX;
	FT ft_y_min = FLT_MAX, ft_y_max = -FLT_MAX;

	for (std::map<int, Planar_Segment *>::iterator it_s = borders.begin() ; it_s != borders.end() ; it_s++) {
		Planar_Segment* s = it_s->second;
		FT x_a = s->A.x(), y_a = s->A.y();

		if (x_a < ft_x_min) ft_x_min = x_a;
		if (x_a > ft_x_max) ft_x_max = x_a;
		if (y_a < ft_y_min) ft_y_min = y_a;
		if (y_a > ft_y_max) ft_y_max = y_a;

		FT x_b = s->B.x(), y_b = s->B.y();
		if (x_b < ft_x_min) ft_x_min = x_b;
		if (x_b > ft_x_max) ft_x_max = x_b;
		if (y_b < ft_y_min) ft_y_min = y_b;
		if (y_b > ft_y_max) ft_y_max = y_b;
	}

	x_min = to_double(ft_x_min);
	x_max = to_double(ft_x_max);
	y_min = to_double(ft_y_min);
	y_max = to_double(ft_y_max);

	dx = x_max - x_min;
	dy = y_max - y_min;
}



void Support_Plane::init_colormap(std::map<Intersection_Line *, CGAL_Color> & colors)
{
	std::default_random_engine generator;
	std::uniform_int_distribution<int> distribution(100, 255);

	for (std::map<int, Intersection_Line*>::iterator it_il = lines.begin() ; it_il != lines.end() ; it_il++) {
		int r = distribution(generator);
		int g = distribution(generator);
		int b = distribution(generator);
		colors[it_il->second] = CGAL_Color(r, g, b);
	}
}



void Support_Plane::draw_grid(std::ostream & os, int rows, int cols, double x_min, double y_min, double dx, double dy, int font_size, double x_s, double y_s)
{
	CGAL_Color gray = CGAL_Color(200, 200, 200);

	double x_max = x_min + dx, y_max = y_min + dy;
	int k_x_min = int(ceil(x_min / x_s)), k_x_max = int(floor(x_max / x_s));
	int k_y_min = int(ceil(y_min / y_s)), k_y_max = int(floor(y_max / y_s));

	// Prints a set of vertical lines
	for (int k = k_x_min ; k <= k_x_max ; k++) {
		double x = cols * (k * x_s - x_min) / dx;
		SVG::markup_dashed_line(os, 10, 5, x, x, 0, rows, gray);
		SVG::markup_text(os, x + 2, rows - 2, k * x_s, font_size, gray);
	}

	// Prints a set of horizontal lines
	for (int k = k_y_min ; k <= k_y_max ; k++) {
		double y = rows * (1 - (k * y_s - y_min) / dy);
		SVG::markup_dashed_line(os, 10, 5, 0, cols, y, y, gray);
		SVG::markup_text(os, 0 + 2, y - 2, k * y_s, font_size, gray);
	}
}



void Support_Plane::draw_bounding_polygon(std::ostream & os, int rows, int cols, double x_min, double y_min, double dx, double dy, int font_size, std::map<Intersection_Line*, CGAL_Color> & colors)
{
	for (std::map<int, Planar_Segment *>::iterator it_s = borders.begin() ; it_s != borders.end() ; it_s++) {
		
		Planar_Segment* s = it_s->second;
		double x1 = cols * (to_double(s->A.x()) - x_min) / dx;
		double x2 = cols * (to_double(s->B.x()) - x_min) / dx;
		double y1 = rows * (1 - (to_double(s->A.y()) - y_min) / dy);
		double y2 = rows * (1 - (to_double(s->B.y()) - y_min) / dy);

		CGAL_Color color = colors[s->support];
		SVG::markup_line(os, x1, x2, y1, y2, color);
	}
}



void Support_Plane::draw_intersection_lines(std::ostream & os, int rows, int cols, double x_min, double y_min, double dx, double dy, int font_size, std::map<Intersection_Line*, CGAL_Color> & colors)
{
	double y_max = y_min + dy;

	for (std::map<int, Intersection_Line *>::iterator it_il = lines.begin() ; it_il != lines.end() ; it_il++) {
	
		Intersection_Line* I = it_il->second;
		const CGAL_Line_2 & I_line = I->line;
		CGAL_Color & I_col = colors[I];

		std::string I_legend = "L(" + std::to_string(I->id_object) + ") [";
		std::list<int>::iterator it_p = I->planes.begin();
		while (true) {
			I_legend += std::to_string(*it_p);
			if (++it_p != I->planes.end()) {
				I_legend += ", ";
			} else {
				I_legend += "] ";
				break;
			}
		}

		// Computes the restriction to the domain D = [x_min, x_max] x [y_min, y_max]
		// We compute the intersection I_line with the four lines defining the domain,
		// this should normally result in a segment [AB]

		bool A_exists = false, B_exists = false;
		CGAL_Point_2 A, B;

		for (std::map<int, Planar_Segment*>::iterator it_pl_s = borders.begin(); it_pl_s != borders.end(); it_pl_s++) {
			Planar_Segment* pl_s = it_pl_s->second;
			CGAL_Segment_2 S (pl_s->A, pl_s->B);

			CGAL::cpp11::result_of<K::Intersect_2(CGAL_Line_2, CGAL_Segment_2)>::type object = intersection(I_line, S);
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
		
		if (A_exists && B_exists) {
			double x_a = to_double(A.x()), y_a = to_double(A.y());
			double x_b = to_double(B.x()), y_b = to_double(B.y());

			// Obtains the final coordinates
			double x1 = cols * (x_a - x_min) / dx;
			double x2 = cols * (x_b - x_min) / dx;
			double y1 = rows * (1 - (y_a - y_min) / dy);
			double y2 = rows * (1 - (y_b - y_min) / dy);
			SVG::markup_dashed_line(os, 5, 5, x1, x2, y1, y2, I_col);
			SVG::markup_text(os, (x1 + x2) / 2 + 10, (y1 + y2) / 2 + 10, I_legend, font_size, I_col);
		}
	}
}



void Support_Plane::draw_polygon_segments(std::ostream & os, int rows, int cols, double x_min, double y_min, double dx, double dy, int font_size, std::map<Intersection_Line*, CGAL_Color> & colors, FT t, double dt)
{
	for (std::map<int, Polygon_Segment *>::iterator it_s = segments.begin() ; it_s != segments.end() ; it_s++) {
		Polygon_Segment* s = it_s->second;
		if (t < s->t_init) continue;

		FT u = jmin(t, s->t_stop);

		CGAL_Point_2 O = s->origin();
		double x1 = cols * (to_double(O.x()) - x_min) / dx;
		double y1 = rows * (1 - (to_double(O.y()) - y_min) / dy);

		CGAL_Point_2 P = s->pt(u);
		double x2 = cols * (to_double(P.x()) - x_min) / dx;
		double y2 = rows * (1 - (to_double(P.y()) - y_min) / dy);

		CGAL_Color color = colors[s->support];
		SVG::markup_line(os, x1, x2, y1, y2, color);

		std::string s_legend = "S(" + std::to_string(s->id_object) + ")";
		SVG::markup_text(os, (x1 + x2) / 2 + 10, (y1 + y2) / 2 + 10, s_legend, font_size, color);

		if (s->t_stop == FLT_MAX) {

			// If the segment is running, we also draw its speed
			Segment_Translation* Tr = nullptr;
			for (std::list<Segment_Translation*>::iterator it_t = s->Tr_previous.begin() ; it_t != s->Tr_previous.end() ; it_t++) {
				if ((*it_t)->t_int_start <= t && t <= (*it_t)->t_int_end) {
					Tr = (*it_t);
					break;
				}
			}
			if (Tr == nullptr) Tr = s->Tr;

			CGAL_Point_2 Q = P + dt * Tr->dA;
			double x3 = cols * (to_double(Q.x()) - x_min) / dx;
			double y3 = rows * (1 - (to_double(Q.y()) - y_min) / dy);

			std::string arrow_id = "arrow-" + std::to_string(s->support->id_object);
			SVG::markup_line(os, x2, x3, y2, y3, color, arrow_id);
		}
	}
}



void Support_Plane::draw_polygon(std::ostream & os, int rows, int cols, double x_min, double y_min, double dx, double dy, 
	int font_size, std::map<Intersection_Line*, CGAL_Color> & colors, FT t, Polygon_Set* T, double dt)
{
	for (std::map<Signature, Polygon_Cell*, Vector_Bool_Comparator>::iterator it_c = T->cells_begin() ; it_c != T->cells_end() ; it_c++) {
		Polygon_Cell* C = it_c->second;
		for (std::list<Polygon*>::iterator it_p = C->polygons_begin() ; it_p != C->polygons_end() ; it_p++) {

			Polygon* P = (*it_p);
			bool has_active_vertices = false;
			for (std::list<Polygon_Vertex*>::iterator it_v = P->vertices.begin(); it_v != P->vertices.end(); it_v++) {
				Polygon_Vertex* v = (*it_v);
				if (v->is_active) {
					has_active_vertices = true;
					break;
				}
			}

			//if (!has_active_vertices) continue;

			for (std::list<Polygon_Edge *>::iterator it_e = P->edges.begin(); it_e != P->edges.end(); it_e++) {
				Polygon_Edge* e = (*it_e);
				std::string e_legend = "E(" + std::to_string(e->id_object) + ")";

				CGAL_Point_2 pt_1 = e->v1->pt(t);
				CGAL_Point_2 pt_2 = e->v2->pt(t);
				double x1 = cols * (to_double(pt_1.x()) - x_min) / dx;
				double x2 = cols * (to_double(pt_2.x()) - x_min) / dx;
				double y1 = rows * (1 - (to_double(pt_1.y()) - y_min) / dy);
				double y2 = rows * (1 - (to_double(pt_2.y()) - y_min) / dy);
				CGAL_Color color = (e->is_constrained ? colors[e->constraint.first] : CGAL::BLACK);

				// SVG::markup_line(os, x1, x2, y1, y2, color);
				SVG::markup_text(os, (x1 + x2) / 2 + 10, (y1 + y2) / 2 + 10, e_legend, font_size, CGAL::BLACK);
			}

			std::vector<CGAL_Point_2> V, CH_V;
			for (std::list<Polygon_Vertex*>::iterator it_v = P->vertices.begin(); it_v != P->vertices.end(); it_v++) {
				Polygon_Vertex* v = (*it_v);
				std::string v_legend = "V(" + std::to_string(v->id_object) + ")";

				CGAL_Point_2 pt_1 = v->pt(t);
				double x1 = cols * (to_double(pt_1.x()) - x_min) / dx;
				double y1 = rows * (1 - (to_double(pt_1.y()) - y_min) / dy);
				SVG::markup_point(os, x1, y1, 2, "black", 0.5, (v->is_active ? "black" : "white"));
				V.push_back(CGAL_Point_2(x1, y1));
				if (v->is_active) {
					// Also draws the direction of the point
					CGAL_Point_2 pt_2 = v->pt(t + dt);
					double x2 = cols * (to_double(pt_2.x()) - x_min) / dx;
					double y2 = rows * (1 - (to_double(pt_2.y()) - y_min) / dy);
					SVG::markup_line(os, x1, x2, y1, y2, CGAL::BLACK, "arrow-black");
				}

				SVG::markup_text(os, x1, y1, v_legend, font_size, CGAL::BLACK);
			}

			CGAL::convex_hull_2(V.begin(), V.end(), std::back_inserter(CH_V));
			CGAL_Color fill = CGAL_Color(200, 191, 231);
			SVG::markup_polygon(os, CH_V, fill, 1, 0.8);
		}
	}
}



void Support_Plane::init_markers(std::ostream & os, std::map<Intersection_Line *, CGAL_Color> & colors)
{
	SVG::markup_defs_header(os);
	SVG::markup_marker(os, "arrow-black", CGAL::BLACK);

	for (std::map<Intersection_Line *, CGAL_Color>::iterator it_c = colors.begin() ; it_c != colors.end() ; it_c++) {
		std::string arrow_id = "arrow-" + std::to_string(it_c->first->id_object);
		SVG::markup_marker(os, arrow_id, it_c->second);
	}

	SVG::markup_defs_footer(os);
}



void Support_Plane::draw(FT t, double dt, int size, double x_step, double y_step, int font_size)
{
	// Obtains the dimensions of the plane
	double x_min, x_max, y_min, y_max, dx, dy;
	init_bounding_square(x_min, y_min, x_max, y_max, dx, dy);

	// Randomly associates a color to each intersection line
	std::map<Intersection_Line*, CGAL_Color> colors;
	init_colormap(colors);

	// Opens a SVG file
	int int_t = int(1000000 * to_double(t));
	std::string filename = "T_" + std::to_string(int_t) + "_P_" + std::to_string(id) + ".svg";

	// Don't override existing files
	FILE* existing_file = fopen(filename.c_str(), "r");
	if (existing_file != NULL) {
		fclose(existing_file);
		return;
	}

	std::filebuf fb;
	fb.open(filename, std::ios::out);
	std::ostream os(&fb);

	// Computes dimensions of the image
	int lx = int(size * dx / jmin(dx, dy));
	int ly = int(size * dy / jmin(dx, dy));
	SVG::markup_header(os, ly, lx);

	// Defines markers
	init_markers(os, colors);

	// Step 1 : prints a background grid
	draw_grid(os, ly, lx, x_min, y_min, dx, dy, font_size, x_step, y_step);

	// Step 2 : draws the intersection lines
	draw_intersection_lines(os, ly, lx, x_min, y_min, dx, dy, font_size, colors);

	// Step 3 : draws the bounding polygon
	draw_bounding_polygon(os, ly, lx, x_min, y_min, dx, dy, font_size, colors);

	// Step 4 : draws the segments
	draw_polygon_segments(os, ly, lx, x_min, y_min, dx, dy, font_size, colors, t, dt);

	// Step 5 : draws the polygon
	draw_polygon(os, ly, lx, x_min, y_min, dx, dy, font_size, colors, t, polygon_set, dt);

	// Closes the file
	SVG::markup_footer(os);
	fb.close();
}
}