#include "stdafx.h"
#include "Geometry.h"

#include "gtest/gtest.h"


using namespace std;
using namespace Geometry;

TEST(GeometryTest, isInsideTriangle) {
	Vector3 p (0,0,1);
	Vector3 a (1,0,1);
	Vector3 b (1,1,1);
	Vector3 c (0,1,1);

	EXPECT_FALSE(isInsideTriangle(p,a,b,c));
	EXPECT_TRUE(isInsideTriangle(a,a,b,c));
	EXPECT_TRUE(isInsideTriangle(b,a,b,c));
	EXPECT_TRUE(isInsideTriangle(c,a,b,c));

	EXPECT_FALSE(isInsideTriangle(a-Vector3(1e-9,0,0),a,b,c));
}

TEST(GeometryTest, projectOntoFace) {
	Vector3 p (0,0,1);
	

	Vector3 a (1,0,1);
	Vector3 b (1,1,1);
	Vector3 c (0,1,1);

	Vector3 projPt;
	Vector3 projNormal;
	EXPECT_FALSE(projectOntoFace(p,a,b,c,projPt,projNormal));

	EXPECT_TRUE(projectOntoFace(Vector3(1,0,1),a,b,c,projPt,projNormal));
	DLOG(INFO)<< projPt << projNormal ;

	EXPECT_TRUE(projectOntoFace(Vector3(1,0,2),a,b,c,projPt,projNormal));
	DLOG(INFO)<< projPt << projNormal ;
}

TEST(GeometryTest, projectOntoEdge) {
	Vector3 projPt;
	Float projLength;
	EXPECT_TRUE(projectOntoEdge(Vector3(0.5,0,0), Vector3(0,0,0), Vector3(1,0,0), projPt, projLength));
	DLOG(INFO)<< projPt << projLength;

	EXPECT_FALSE(projectOntoEdge(Vector3(1+1e-6,0,0), Vector3(0,0,0), Vector3(1,0,0), projPt, projLength));
	DLOG(INFO)<< projPt << projLength;

	EXPECT_TRUE(projectOntoEdge(Vector3(0.25,0,1), Vector3(0,0,0), Vector3(1,0,0), projPt, projLength));
	DLOG(INFO)<< projPt << projLength;

}

TEST(GeometryTest, PlaneTriangle4D) {
	Vertex4 a(Vector4(0,0,0,0), 0);
	Vertex4 b(Vector4(1,0,0,0), 1);
	Vertex4 c(Vector4(0,1,0,0), 2);

	vector<Intersection4> intersections;
	doesPlaneIntersectTriangle(0.0, a, b, c, intersections);

	EXPECT_EQ(3, intersections.size());

}


TEST(GeometryTest, PlaneTet4D) {
	Vertex4 a(Vector4(0,0,0,0), 0);
	Vertex4 b(Vector4(1,0,0,0), 1);
	Vertex4 c(Vector4(0,1,0,0), 2);
	Vertex4 d(Vector4(1,0,0,1), 3);

	vector<Triangle> triangles;
	
	doesPlaneIntersectTetrahedron(0.5, a, b, c, d, triangles);

	EXPECT_EQ(1, triangles.size());
	//for (int i = 0; i < triangles.size(); i++) {
	//	cout << triangles[i].a << " " << triangles[i].b << " " << triangles[i].c << endl;
	//}
	

}