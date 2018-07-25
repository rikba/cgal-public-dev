#pragma once
#include "defs_cgal.h"

namespace JPTD {

namespace Ply_In
{
	
	void read(const std::string & filename, std::vector<std::vector<CGAL_Point_3> > & polygons);

	void read(const std::string & filename, std::vector<CGAL_Point_3> & points);

	bool get_words(std::ifstream & file, std::vector<std::string> & words);

	void get_number_of_vertices_and_facets(std::ifstream & file, int & V, int & F);

	void get_vertices(std::ifstream & file, const int V, std::vector<CGAL_Point_3> & vertices);

	void get_facets(std::ifstream & file, const int F, const std::vector<CGAL_Point_3> & vertices, 
		std::vector<std::vector<CGAL_Point_3> > & polygons);
};

}