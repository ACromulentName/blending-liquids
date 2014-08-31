#include "stdafx.h"
#include "Mesh3D.h"

#include "gtest/gtest.h"
#include <array>
using namespace std;

TEST(Mesh3DTest, LoadObj) {
	Mesh3D mesh;
	EXPECT_TRUE(mesh.loadObj("tests/cube.obj"));

	EXPECT_EQ(mesh.vertices.size(), 8);
	EXPECT_EQ(mesh.faces.size(), 12);

	array<int,5> expected = {4,6,2,3,1};
	vector<int> verts = mesh.getOneRing(0);
	for (int i = 0; i < verts.size(); i++) {
		EXPECT_EQ(expected[i], verts[i]);
	}
	
}

TEST(Mesh3DTest, LoadBinaryObj) {
	Mesh3D mesh;
	EXPECT_TRUE(mesh.loadBinaryObj("tests/basic3d_0000.bobj.gz"));

	EXPECT_EQ(mesh.vertices.size(), 20600);
	EXPECT_EQ(mesh.faces.size(), 41196);

}

TEST(Mesh3DTest, WriteBinaryObj) {
	Mesh3D mesh;
	EXPECT_TRUE(mesh.loadBinaryObj("tests/basic3d_0000.bobj.gz"));

	EXPECT_EQ(mesh.vertices.size(), 20600);
	EXPECT_EQ(mesh.faces.size(), 41196);


	mesh.writeBinaryObj("objTest.bobj.gz");

	Mesh3D newMesh;
	EXPECT_TRUE(newMesh.loadBinaryObj("objTest.bobj.gz"));
	EXPECT_EQ(newMesh.vertices.size(), 20600);
	EXPECT_EQ(newMesh.faces.size(), 41196);
}

TEST(Mesh3DTest, Visit) {
	Mesh3D mesh;
	EXPECT_TRUE(mesh.loadObj("tests/cube.obj"));

	EXPECT_EQ(mesh.vertices.size(), 8);
	EXPECT_EQ(mesh.faces.size(), 12);
	
	int face = mesh.visit(true);
	while (face != -1) {
		cout << "Visited face:" << face << endl;
		face = mesh.visit();
	}

}
