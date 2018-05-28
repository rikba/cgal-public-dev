#include <iostream>
#include <fstream>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_2 Point_2;

int main() {
    std::ifstream stream("../data/points.xy");
    CGAL::set_ascii_mode(stream);
    Point_2 a, b, c;
    stream >> a >> b >> c;
    std::cout << b << std::endl;
    std::cout << (CGAL::collinear(a, b, c) ? "Collinear" : "Non-collinear") << std::endl;
}
