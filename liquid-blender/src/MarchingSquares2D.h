#include "stdafx.h"
#include "Animation.h"

template<int DIM>
struct SDF {
	SDF() {}

	void computeSDF(const Animation<DIM>& anim, int time, Float edgeLength); 

	bool isInside(int x, int y) const {
		return (field[numPointsX*y + x] < 0);
	}

	bool isInside(int x, int y, int z) const {
		return (field[(numPointsX*numPointsX)*z + numPointsX*y + x] < 0);
	}

	Float value(int x, int y) const {
		return field[ numPointsX * y + x];
	}

	vector < Float > field;
	int numPointsX; 
	Float edgeLength;
};

template<int DIM>
class MarchingSquares {
public:
	MarchingSquares(const SDF<DIM>& sdf): sdf(sdf) {
	}

	void extract(string filename);
	
	Vector2 computeIntersection(int x0, int y0, int x1, int y1) {
		Float a = sdf.value(x0,y0);
		Float b = sdf.value(x1,y1);
		Float t = a / (a-b);
		Vector2 v0 = Vector2(x0 * sdf.edgeLength,y0 * sdf.edgeLength);
		Vector2 v1 = Vector2(x1 * sdf.edgeLength,y1 * sdf.edgeLength);
		return (1-t)*v0 + (t)*v1;
	}

protected:
	const SDF<DIM>& sdf;
};