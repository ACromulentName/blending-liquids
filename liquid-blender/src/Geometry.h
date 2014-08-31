#pragma once
#include "stdafx.h"

using namespace  std;
namespace Geometry
{


	struct Vertex4 {

		Vertex4() {}
		Vertex4(Vector4 p, int _index) {
			pos = p;
			index = _index;
		}

		Vector4 pos;
		int index;
	};

	struct Intersection4 {
		Vector4 pos;
		int index0;
		int index1;
	};

	struct Triangle {
		Intersection4 a,b,c;
		Vector3 normal;
		int origTriIndex;
		int type;

		Triangle() {
			origTriIndex = -1;
		}

	};

	bool projectOntoEdge(const Vector3& p, const Vector3&a, const Vector3&b, Vector3& projPt, Float& projLength);

	template<class T> 
	T lerp(const T& a, const T& b, const Float t) {
		return (1-t)*a + t*b;
	}

	Vector3 nlerp(const Vector3&, const Vector3&, const Float t);
	inline bool isSamePoint(const Intersection4& a, const Intersection4& b);

	inline bool onSameSide(const Vector3& p1, const Vector3& p2, const Vector3& a, const Vector3& b);
	inline bool isInsideTriangle(const Vector3& p, const Vector3& a, const Vector3& b, const Vector3& c);
	inline void projectOntoPlane(const Vector3& p, const Vector3& origin, const Vector3& normal, Vector3& projPt);
	bool projectOntoFace(const Vector3& p, const Vector3& a, const Vector3& b, const Vector3& c, Vector3& projPt, Vector3& normal);

	bool doesPlaneIntersectEdge( Float time, const Vector3& startEdge, const Vector3& endEdge, Vector3& intersection );

	bool doesPlaneIntersectEdge(Float time, const Vertex4& a, const Vertex4& b, Intersection4& intersection);

	inline bool doesPlaneIntersectTriangle(Float time, const Vertex4& a, const Vertex4& b, const Vertex4& c, vector<Intersection4>& intersections);

	bool doesPlaneIntersectTetrahedron(Float time, const Vertex4& a, const Vertex4& b, const Vertex4& c, const Vertex4& d, vector<Triangle>& triangles, int origTri=-1, int type=-1);

	Vector4 computeOrthogonalVector(const Vector4& st0, const Vector4& st1,const Vector4& st2);
	bool doesRayIntersectTriangle(const Vector3& o, const Vector3& d, const Vector3&a, const Vector3&b, const Vector3&c, Vector3&);
}