#include "ply_out.h"

namespace JPTD {
void Ply_Out::print_colorless_vertices(const std::string & name, const std::list<CGAL_Point_3> & vertices, int precision)
{
	std::ofstream stream(name, std::ofstream::out);
	if (stream.is_open()) {

		// Writes the header
		stream << "ply" << std::endl
			<< "format ascii 1.0" << std::endl
			<< "element vertex " << vertices.size() << std::endl
			<< "property float x" << std::endl
			<< "property float y" << std::endl
			<< "property float z" << std::endl
			<< "end_header" << std::endl;

		// Writes the list of vertices
		stream << std::setprecision(precision);
		for (std::list<CGAL_Point_3>::const_iterator it_d = vertices.begin(); it_d != vertices.end(); it_d++) {
			stream << it_d->x() << " " << it_d->y() << " " << it_d->z() << std::endl;
		}

		// Closes stream
		stream.close();
	}
}



void Ply_Out::print_colorful_vertices(const std::string & name, const std::list<CGAL_Point_3> & vertices, const std::list<CGAL_Color> & colors, int precision)
{
	assert(vertices.size() == colors.size());

	std::ofstream stream(name, std::ofstream::out);
	if (stream.is_open()) {

		// Writes the header
		stream << "ply" << std::endl
			<< "format ascii 1.0" << std::endl
			<< "element vertex " << vertices.size() << std::endl
			<< "property float x" << std::endl
			<< "property float y" << std::endl
			<< "property float z" << std::endl
			<< "property uchar red" << std::endl
			<< "property uchar green" << std::endl
			<< "property uchar blue" << std::endl
			<< "end_header" << std::endl;

		// Writes the list of vertices
		std::list<CGAL_Color>::const_iterator it_c = colors.begin();
		std::list<CGAL_Point_3>::const_iterator it_d = vertices.begin();
		while (it_d != vertices.end()) {
			stream << std::setprecision(precision) << it_d->x() << " " << it_d->y() << " " << it_d->z();
			stream << std::setprecision(3) << " " << int(it_c->red()) << " " << int(it_c->green()) << " " << int(it_c->blue()) << std::endl;
			++it_c, ++it_d;
		}

		// Closes stream
		stream.close();
	}
}



void Ply_Out::print_colorful_facets(const std::string & name, const std::list<CGAL_Point_3> & vertices,
	const std::list<std::list<int> > & facets, const std::list<CGAL_Color> & colors, int precision)
{
	// Computes the number of vertices and facets
	int nb_facets = int(facets.size());
	int nb_vertices = int(vertices.size());
	int nb_colors = int(colors.size());

	std::ofstream stream(name, std::ofstream::out);
	if (stream.is_open()) {

		// Header
		stream << "ply" << std::endl
			<< "format ascii 1.0" << std::endl;

		// Case when there is a color per vertex
		// The color of the facet is then obtained by interpolation

		stream << "element vertex " << nb_vertices << std::endl
			<< "property float x" << std::endl
			<< "property float y" << std::endl
			<< "property float z" << std::endl
			<< "property uchar red" << std::endl
			<< "property uchar green" << std::endl
			<< "property uchar blue" << std::endl;

		stream << "element face " << nb_facets << std::endl
			<< "property list uchar int vertex_index" << std::endl
			<< "end_header" << std::endl;

		// Writes the list of vertices
		std::list<CGAL_Point_3>::const_iterator it_v = vertices.begin();
		std::list<CGAL_Color>::const_iterator it_c = colors.begin();
		while (it_v != vertices.end() && it_c != colors.end()) {
			stream << std::setprecision(12) << it_v->x() << " " << it_v->y() << " " << std::setprecision(6) << it_v->z() << " ";
			stream << std::setprecision(3) << int(it_c->red()) << " " << int(it_c->green()) << " " << int(it_c->blue()) << std::endl;
			++it_v; ++it_c;
		}

		// Writes the list of facets
		for (std::list<std::list<int> >::const_iterator it_i1 = facets.begin(); it_i1 != facets.end(); it_i1++) {
			const std::list<int> & facet = (*it_i1);
			std::list<int>::const_iterator it_i2 = facet.begin();
			stream << facet.size() << " " << (*it_i2);
			while (++it_i2 != facet.end()) stream << " " << (*it_i2);
			stream << std::endl;
		}

		// Closes file
		stream.close();
	}
}



void Ply_Out::print_plain_colorful_facets(const std::string & name, const std::list<CGAL_Point_3> & vertices, 
	const std::list<std::list<int> > & facets, const std::list<CGAL_Color> & colors, int precision)
{
	// Computes the number of vertices and facets
	int nb_facets = int(facets.size());
	int nb_vertices = int(vertices.size());
	int nb_colors = int(colors.size());

	std::ofstream stream(name, std::ofstream::out);
	if (stream.is_open()) {

		// Header
		stream << "ply" << std::endl
			<< "format ascii 1.0" << std::endl;

		// Facets are plain, to the color attribute goes to the facet
		stream << "element vertex " << nb_vertices << std::endl
			<< "property float x" << std::endl
			<< "property float y" << std::endl
			<< "property float z" << std::endl;
		stream << "element face " << nb_facets << std::endl
			<< "property list uchar int vertex_index" << std::endl
			<< "property uchar red" << std::endl
			<< "property uchar green" << std::endl
			<< "property uchar blue" << std::endl
			<< "end_header" << std::endl;

		// Prints the list of vertices, which only have coordinates
		stream << std::setprecision(precision);
		for (std::list<CGAL_Point_3>::const_iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
			stream << it_v->x() << " " << it_v->y() << " " << it_v->z() << " " << std::endl;
		}

		// Prints the list of facets
		std::list<std::list<int> >::const_iterator it_i1 = facets.begin();
		std::list<CGAL_Color>::const_iterator it_c = colors.begin();

		while (it_i1 != facets.end() && it_c != colors.end()) 
		{
			// Writes the definition of the facet
			const std::list<int> & facet = (*it_i1);
			std::list<int>::const_iterator it_i2 = facet.begin();

			stream << facet.size() << " " << (*it_i2);
			while (++it_i2 != facet.end()) stream << " " << (*it_i2);

			// Writes its color
			stream << " " << int(it_c->red()) << " " << int(it_c->green()) << " " << int(it_c->blue()) << std::endl;

			// Iterates
			++it_i1; ++it_c;
		}

		// Closes file
		stream.close();
	}
}



void Ply_Out::print_plain_colorful_facets(const std::string & name, const std::list<std::list<CGAL_Point_3> > & polygons,
	const std::list<CGAL_Color> & colors, int precision)
{
	int nb_facets = int(polygons.size()), nb_vertices = 0;
	for (auto it_p = polygons.begin(); it_p != polygons.end(); it_p++) nb_vertices += int(it_p->size());

	// Opens the file
	std::ofstream stream(name, std::ofstream::out);
	if (stream.is_open()) {

		// Writes the header
		stream << "ply" << std::endl
			<< "format ascii 1.0" << std::endl
			<< "element vertex " << nb_vertices << std::endl
			<< "property float x" << std::endl
			<< "property float y" << std::endl
			<< "property float z" << std::endl
			<< "property uchar red" << std::endl
			<< "property uchar green" << std::endl
			<< "property uchar blue" << std::endl;
		stream << "element face " << nb_facets << std::endl
			<< "property list uchar int vertex_index" << std::endl
			<< "end_header" << std::endl;
		
		std::list<std::list<CGAL_Point_3> >::const_iterator it_p;
		std::list<CGAL_Color>::const_iterator it_c;

		// Writes the list of vertices
		stream << std::setprecision(precision);
		for (it_p = polygons.begin(), it_c = colors.begin(); it_p != polygons.end() && it_c != colors.end(); it_p++, it_c++) {
			for (std::list<CGAL_Point_3>::const_iterator it_v = it_p->begin(); it_v != it_p->end(); it_v++) {
				stream << it_v->x() << " " << it_v->y() << " " << it_v->z() << " " << int(it_c->red()) << " " << int(it_c->green()) << " " << int(it_c->blue()) << std::endl;
			}
		}

		// Writes the list of facets, defined by a sequence of indices
		int offset = 0;
		for (it_p = polygons.begin() ; it_p != polygons.end() ; it_p++) {
			int n = (int)it_p->size();
			stream << n << " ";
			for (int i = 0 ; i < n ; i++) {
				stream << offset + i;
				stream << (i != n - 1 ? ' ' : '\n');
			}
			offset += n;
		}

		// Closes stream
		stream.close();
	}
}


void Ply_Out::print_colorless_facets(const std::string & name, const std::list<CGAL_Point_3> & vertices, const std::list<std::list<int> > & facets, int precision)
{
	// Computes the number of vertices and facets
	int nb_facets = int(facets.size());
	int nb_vertices = int(vertices.size());

	std::ofstream stream(name, std::ofstream::out);
	if (stream.is_open()) {

		// Header
		stream << "ply" << std::endl
			<< "format ascii 1.0" << std::endl;
		stream << "element vertex " << nb_vertices << std::endl
			<< "property float x" << std::endl
			<< "property float y" << std::endl
			<< "property float z" << std::endl;
		stream << "element face " << nb_facets << std::endl
			<< "property list uchar int vertex_index" << std::endl
			<< "end_header" << std::endl;

		// Writes the list of vertices
		stream << std::setprecision(precision);
		for (std::list<CGAL_Point_3>::const_iterator it_v = vertices.begin(); it_v != vertices.end(); it_v++) {
			 stream << it_v->x() << " " << it_v->y() << " " << it_v->z() << std::endl;
		}

		// Writes the list of facets
		for (std::list<std::list<int> >::const_iterator it_i1 = facets.begin(); it_i1 != facets.end(); it_i1++) {
			const std::list<int> & facet = (*it_i1);

			std::list<int>::const_iterator it_i2 = facet.begin();

			stream << facet.size() << " " << (*it_i2);
			while (++it_i2 != facet.end()) stream << " " << (*it_i2);
			stream << std::endl;
		}

		// Closes file
		stream.close();
	}
}
}