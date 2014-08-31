#include "stdafx.h"
#include "Geometry.h"
#include "Utility.hpp"


namespace Geometry
{

/*
Computes projection of point p onto edge ab.
Returns true if the projected point lies on AB (closed). 
Stores projected point and length in the last two arguments
*/
bool projectOntoEdge(const Vector3& p, const Vector3&a, const Vector3&b, Vector3& projPt, Float& interpolant) {
	Float lengthAB = (b-a).norm();
	Vector3 unitAB = (b-a) / lengthAB;
	Float projLength = (p - a).dot(unitAB);
	projPt = a + projLength * unitAB;
	interpolant = projLength / lengthAB;
	return (projLength >= 0 && projLength <= lengthAB);
}

/*
Interpolates between two vectors, normalizes the result
*/
Vector3 nlerp (const Vector3& a, const Vector3& b, Float t) {
	Vector3 n = Geometry::lerp<Vector3> (a, b, t);
	n.normalize();
	return n;
}

/*
Returns true if p1 and p2 lie on the same side of line segment ab
*/
bool onSameSide(const Vector3& p1, const Vector3& p2, const Vector3& a, const Vector3& b) {
	return (  ((b-a).cross(p1-a)).dot((b-a).cross(p2-a)) >= 0 );
}

/*
Checks if p lies inside triangle abc. Note: p needs to be in the plane of abc for this to be meaningful.
*/
bool isInsideTriangle(const Vector3& p, const Vector3& a, const Vector3& b, const Vector3& c) {
	return ( onSameSide(p,a,b,c) && onSameSide(p,b,a,c) && onSameSide(p,c,a,b) );
}

/*
Projects p onto plane defined by (origin,normal). Stores point in projPt.
*/
void projectOntoPlane(const Vector3& p, const Vector3& origin, const Vector3& normal, Vector3& projPt) {
	projPt = p - (p-origin).dot(normal) * normal;
}

/* 
Projects p onto face defined by abc. Returns true if p lies inside triangle abc. Stores projected point and normal in projPt and normal respectively.
*/
bool projectOntoFace(const Vector3& p, const Vector3& a, const Vector3& b, const Vector3& c, Vector3& projPt, Vector3& normal) {
	normal = (b-a).cross(c-a);
	normal.normalize();
	projectOntoPlane(p, a, normal, projPt);

	return isInsideTriangle(projPt, a, b, c);
}

/*
Tests if 3D plane parallel to time axis, intersects the edge ab at time t
*/
bool doesPlaneIntersectEdge( Float time, const Vector3& startEdge, const Vector3& endEdge, Vector3& intersection )
{
	if (endEdge[2] == startEdge[2])
		return false;

	Float alpha = (time - startEdge[2]) / (endEdge[2] - startEdge[2]);

	if (alpha >= 0 && alpha <= 1) {
		intersection = (1-alpha) * startEdge + (alpha) * endEdge;
		return true;
	}

	return false;
}


/* 
Tests if the 4D hyperplane at time intersects edge ab 
*/
bool doesPlaneIntersectEdge(Float time, const Vertex4& a, const Vertex4& b, Intersection4& intersection) {
	Float eps = 1e-10;
	if (std::abs(a.pos[3]-b.pos[3])<eps && std::abs(time-a.pos[3])<eps) {
		intersection.pos = a.pos;
		intersection.index0 = min(a.index, b.index);
		intersection.index1 = max(a.index, b.index);
		return true;
	}

	Float alpha = (time - a.pos[3]) / (b.pos[3] - a.pos[3]);
	if (alpha >= 0 && alpha < 1) {
		intersection.pos = (1-alpha)*a.pos + alpha * b.pos;
		intersection.index0 = min(a.index, b.index);
		intersection.index1 = max(a.index, b.index);

		return true;
	}

	return false;
}

/*
Tests if the 4D hyperplane at time t intersects triangle defined by abc. 
*/
bool doesPlaneIntersectTriangle(Float time, const Vertex4& a, const Vertex4& b, const Vertex4& c, vector<Intersection4>& intersections) {
	bool hasIntersection = false;
	intersections.clear();

	Intersection4 intersection;
	if (doesPlaneIntersectEdge(time, a, b, intersection)) {
		intersections.emplace_back(intersection);
		hasIntersection = true;
	}
	
	if (doesPlaneIntersectEdge(time, b, c, intersection)) {
		intersections.emplace_back(intersection);		
		hasIntersection = true;
	}
	
	if (doesPlaneIntersectEdge(time, c, a, intersection)) {
		intersections.emplace_back(intersection);
		hasIntersection = true;
	}

	return hasIntersection;
}

bool isSamePoint(const Intersection4& a, const Intersection4& b) {
	if (a.index0 == b.index0 && a.index1 == b.index1)
		return true;
	return false;
}

/*
Tests if the 4D hyperplane intersects the tetrahedron defined by a,b,c,d. 
Returns triangulated output
*/
bool doesPlaneIntersectTetrahedron(Float time, const Vertex4& a, const Vertex4& b, const Vertex4& c, const Vertex4& d, vector<Triangle>& triangles, int origTri, int type) {
	vector< vector<Intersection4> > allIntersections;
	vector<Intersection4> triIntersections;
	triangles.reserve(4);
	allIntersections.reserve(4);

	//LOG(INFO) << " Tet: ";
	//LOG(INFO) << " A: " << a;
	//LOG(INFO) << " B: " << b;
	//LOG(INFO) << " C: " << c;
	//LOG(INFO) << " D: " << d;

	if (doesPlaneIntersectTriangle(time, a, c, b, triIntersections)) {	
		allIntersections.emplace_back(triIntersections);
	}


	if (doesPlaneIntersectTriangle(time, a, b, d, triIntersections)) {
		allIntersections.emplace_back(triIntersections);
	}
	//

	if (doesPlaneIntersectTriangle(time, b, c, d, triIntersections)) {		
		allIntersections.emplace_back(triIntersections);		
	}
	////

	if (doesPlaneIntersectTriangle(time, a, d, c, triIntersections)) {
		allIntersections.emplace_back(triIntersections);		
	}


	vector<Intersection4> vertices;


	
	// Create list of vertices that need to be triangulated
	for (int i = 0; i < allIntersections.size(); i++) {
		if (allIntersections[i].size() == 0)
			continue;

		if (vertices.size() == 0) {
			for (int j = 0; j < allIntersections[i].size(); j++)
				vertices.emplace_back (allIntersections[i][j]);
		}
		else {
			if (allIntersections[i].size() != 2) 
				continue;

			Intersection4 start = allIntersections[i][0];
			Intersection4 end = allIntersections[i][1];

			if (isSamePoint(start, vertices[vertices.size()-1])) {
				if (!isSamePoint(end, vertices[0]))
					vertices.emplace_back(end);
			}
			else if (isSamePoint(end, vertices[vertices.size()-1])) {
				if (!isSamePoint(start, vertices[0]))
					vertices.emplace_back(start);
			}
			else if (isSamePoint(start, vertices[0])) {
				if (!isSamePoint(end, vertices[vertices.size()-1]))
					vertices.insert(vertices.begin(), end);
			}
			else if (isSamePoint(end, vertices[0])) {
				if (!isSamePoint(start, vertices[vertices.size()-1]))
					vertices.insert(vertices.begin(), start);
			}

			//else {
			//	vertices.push_back(start);
			//	vertices.push_back(end);
			//}
		}
	}

	//LOG(INFO) << "# intersections:" << allIntersections.size();
	////if (vertices.size() == 4 || vertices.size() == 3) 
	//
	//{
	//	
	//	for (int i = 0; i < allIntersections.size(); i++) {			
	//		LOG(INFO) << "Intersection: "<< i;
	//		for (int j = 0; j < allIntersections[i].size(); j++) {
	//			LOG(INFO) << allIntersections[i][j];
	//		}
	//	}

	//	LOG(INFO) << "# vertices:" << vertices.size();
	//	for (int i =0; i< vertices.size(); i++) {
	//		LOG(INFO) << "Vertex:" << vertices[i];
	//	}
	//}
	

	// Make triangles
	if (vertices.size() == 3) {
		Triangle t;
		t.a = vertices[0];
		t.b = vertices[1];
		t.c = vertices[2];
		t.origTriIndex = origTri;
		t.type = type;
		triangles.emplace_back(t);
		return true;
	}
	else if (vertices.size() == 4) {
		Triangle t;
		t.a = vertices[0];
		t.b = vertices[1];
		t.origTriIndex = origTri;
		t.c = vertices[2];
		t.type = type;
		triangles.emplace_back(t);

		t.a = vertices[0];
		t.b = vertices[2];
		t.c = vertices[3];
		t.origTriIndex = origTri;
		t.type = type;
		triangles.emplace_back(t);
		return true;
	}

	return false;

}

// Computes *(st0^st1^st2) and returns normalized vector
Vector4 computeOrthogonalVector(const Vector4& st0, const Vector4& st1,const Vector4& st2) {


	Matrix3 mI,mJ,mK,mL;
	mI << st0[1], st0[2], st0[3],
		st1[1], st1[2], st1[3],
		st2[1], st2[2], st2[3];

	mJ << st0[0], st0[2], st0[3],
		st1[0], st1[2], st1[3],
		st2[0], st2[2], st2[3];

	mK << st0[0], st0[1], st0[3],
		st1[0], st1[1], st1[3],
		st2[0], st2[1], st2[3];

	mL << st0[0], st0[1], st0[2],
		st1[0], st1[1], st1[2],
		st2[0], st2[1], st2[2];

	Vector4 n = Vector4(determinant(mI), -determinant(mJ), determinant(mK), -determinant(mL));

	if (n.norm() < 1e-15)  {
		//LOG(INFO) << "Zero normal! " << n << ";" << st0 << ";" << st1 << ";" << st2;
	}
	n.normalize();
	return n;
}

bool doesRayIntersectTriangle(const Vector3& o, const Vector3& d, const Vector3&a, const Vector3&b, const Vector3&c, Vector3& intersectionPoint) {
	Vector3 normal = (b-a).cross(c-a);
	float dr = normal.dot(d);
	if (fabs(dr) == 1e-6) 
		return false;

	Float t = - (o-a).dot(normal) / dr;
	intersectionPoint = o + t * d;
	
	return isInsideTriangle(intersectionPoint, a, b, c);
}

}
