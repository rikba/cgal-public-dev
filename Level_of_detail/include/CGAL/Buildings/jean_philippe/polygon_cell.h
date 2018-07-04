#pragma once
#include "polygon.h"

namespace JPTD {
typedef std::vector<bool> Signature;


class Polygon_Cell
{
public:
	Polygon_Cell(const Signature & S, Polygon* P);

	~Polygon_Cell();

	static Signature make_signature(const CGAL_Point_2 & O, const std::map<int, Intersection_Line*> & L);

	inline Polygon* get_unique_polygon() {
		assert(polygons.size() == 1);
		return polygons.front();
	}

	Polygon* get_one_polygon() const;

	void push(Polygon* P);

	inline Signature & get_signature() { return signature; }

	size_t size() const;
	const std::list<Polygon*>::iterator polygons_begin();
	const std::list<Polygon*>::iterator polygons_end();

protected:
	Signature signature;
	std::list<Polygon*> polygons;
};


struct Vector_Bool_Comparator
{
	bool operator() (const std::vector<bool> & SL, const std::vector<bool> & SR) const
	{
		// Loops on all bits of SL and SR
		// If bits differ, SL < SR iif SL[i] == 0 and SR[i] == 1
		int n = int(SL.size());
		for (int i = 0 ; i < n ; ++i) {
			if (SL[i] != SR[i]) {
				return (!SL[i] && SR[i]);
			}
		}

		// SL == SR
		return false;
	}
};
}