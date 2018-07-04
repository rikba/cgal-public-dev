#pragma once
#include "defs.h"
#include <fstream>
#include <iomanip>

namespace JPTD {
namespace Ply_Out 
{
	void print_colorless_vertices(const std::string & name, const std::list<CGAL_Point_3> & vertices, int precision = 6);

	void print_colorful_vertices(const std::string & name, const std::list<CGAL_Point_3> & vertices, const std::list<CGAL_Color> & colors, int precision = 6);

	void print_colorful_facets(const std::string & name, const std::list<CGAL_Point_3> & vertices, 
		const std::list<std::list<int> > & facets, const std::list<CGAL_Color> & colors, int precision = 6);

	void print_plain_colorful_facets(const std::string & name, const std::list<CGAL_Point_3> & vertices, 
		const std::list<std::list<int> > & facets, const std::list<CGAL_Color> & colors, int precision = 6);

	void print_plain_colorful_facets(const std::string & name, const std::list<std::list<CGAL_Point_3> > & polygons,
		const std::list<CGAL_Color> & colors, int precision = 6);

	void print_colorless_facets(const std::string & name, const std::list<CGAL_Point_3> & vertices, const std::list<std::list<int> > & facets, int precision = 6);
};
}