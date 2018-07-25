#include "polygon.h"

namespace JPTD {

Polygon_Directions::Polygon_Directions(const CGAL_Point_2 & _O, 
	const std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & _vertices,
	const std::vector<Intersection_Line*> & _reference_lines)
	: O (_O)
{
	// Copies vectors
	std::copy(_vertices.begin(), _vertices.end(), std::back_inserter(vertices));
	std::copy(_reference_lines.begin(), _reference_lines.end(), std::back_inserter(reference_lines));
}


#if 0
Polygon_Directions::Polygon_Directions(const CGAL_Point_2 & _O, 
	const std::vector<std::pair<CGAL_Point_2, CGAL_Vector_2> > & unsorted_vertices)
	: O (_O)
{
	// Inserts all elements of 'unsorted_vertices' in the vector 'vertices' such that :
	// - they are sorted by increasing value of theta_i, where theta_i is the angle made by the directions of V[i]
	// - the first element of 'vertices' makes the closest angle with 0 while being positive.

	typedef std::pair<CGAL_Point_2, CGAL_Vector_2> Point_Vector;
	typedef std::pair<double, Point_Vector> Orientation;

	size_t n = unsorted_vertices.size();
	std::vector<Orientation> orientations;
	orientations.reserve(n);

	for (size_t i = 0 ; i < n ; i++) {
		const CGAL_Vector_2 & u_i = unsorted_vertices[i].second;
		double theta_i = atan2(u_i.y(), u_i.x());
		if (theta_i < 0) theta_i += 2 * PI;
		
		orientations.push_back(std::make_pair(theta_i, unsorted_vertices[i]));
	}

	struct _Orientation_Comparator {
		bool operator() (const Orientation & OL, const Orientation & OR) {
			return OL.first < OR.first;
		}
	} Orientation_Comparator;

	// std::sort(orientations.begin(), orientations.end(), Orientation_Comparator);

	for (size_t i = 0 ; i < n ; i++) {
		vertices.push_back(orientations[i].second);
	}
}
#endif


Polygon_Directions::~Polygon_Directions()
{
	vertices.clear();
}


size_t Polygon_Directions::size() const
{
	return vertices.size();
}


std::pair<CGAL_Point_2, CGAL_Vector_2> & Polygon_Directions::get_vertex(int i) 
{
	size_t n = vertices.size();
	return vertices[i % n];
}


Intersection_Line* Polygon_Directions::get_reference_line(int i)
{
	return reference_lines[i];
}


const CGAL_Point_2 & Polygon_Directions::get_barycenter() const
{
	return O;
}

}