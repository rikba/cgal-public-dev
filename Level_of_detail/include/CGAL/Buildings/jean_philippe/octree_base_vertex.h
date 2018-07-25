#pragma once

namespace JPTD {

template <typename Point_3, typename Inexact_Point_3>
class Octree_Base_Vertex
{
public:
	Octree_Base_Vertex(const Point_3 & _M, const Inexact_Point_3 & _hint_M);

	virtual ~Octree_Base_Vertex();

public:
	Point_3 M;
	Inexact_Point_3 hint_M;
};


template <typename Point_3, typename Inexact_Point_3>
Octree_Base_Vertex<Point_3, Inexact_Point_3>::Octree_Base_Vertex(const Point_3 & _M, const Inexact_Point_3 & _hint_M)
{
	M = _M;
	hint_M = _hint_M;
}


template <typename Point_3, typename Inexact_Point_3>
Octree_Base_Vertex<Point_3, Inexact_Point_3>::~Octree_Base_Vertex()
{
}

}