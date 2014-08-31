#include "stdafx.h"
#include "Animation.h"
#include <fstream>
#include "gtest/gtest.h"

Float tol = 1e-6;

TEST(Animation2DTest, Load) {		
	Animation2D anim;
	EXPECT_TRUE(anim.load("tests/drop0/drop0",10, 30));
	GlobalIndex gIndex = anim.frameToGlobalMap[FrameIndex(0,0)];
	EXPECT_EQ(anim.frameToGlobalMap[FrameIndex(0,0)], 0);
}

TEST(Animation2DTest, FindClosestVertex) {
	Animation2D anim;
	anim.load("tests/drop0/drop0",10, 12);

	Vector3 queryPoint(0.5,0.5,0);
	int index = anim.findClosestVertex(queryPoint);	
	FrameIndex fIndex = anim.globalToFrameMap[index];

	EXPECT_EQ(0, fIndex.first);
	EXPECT_EQ(53, fIndex.second);
	
}

TEST(Animation2DTest, FindNeighbours) {	
	Animation2D anim;
	anim.load("tests/drop0/drop0",10, 12);
	int index = anim.frameToGlobalMap[FrameIndex(1,85)];
	EXPECT_EQ(188, anim.neighbours[index][PrevSpace].index);
	EXPECT_EQ(190, anim.neighbours[index][NextSpace].index);
	EXPECT_EQ(291, anim.neighbours[index][NextTime].index);
	EXPECT_EQ(85, anim.neighbours[index][PrevTime].index);
}

TEST(Animation2DTest, ComputeNormals) {	
	Animation2D anim;
	anim.load("tests/drop0/drop0",10, 12);
	EXPECT_NEAR(0.0, (Vector3(1,0,0) - anim.normals[197]).norm(), tol);
	EXPECT_NEAR(0.0, (Vector3(-0.01674694, -0.73940759, 0.67304975) - anim.normals[235]).norm(), tol);
}

TEST(Animation2DTest, ComputeAllNormals) {
	Animation2D animA, animB;
	animA.load("tests/drop0/drop0",0, 10);
	animB.load("tests/drop1/drop1",0, 10);

	fstream f;

	f.open("tests/drop0/normalsDrop0.txt",ios::in);
	for (int i = 0; i < animA.vertices.size(); i++) {
		Vector3 tmpNormal;
		f >> tmpNormal[0] >> tmpNormal[1] >> tmpNormal[2];

		EXPECT_NEAR(0.0, (tmpNormal - animA.normals[i]).norm(), 1e-3);
	}
	f.close();

	f.open("tests/drop1/normalsDrop1.txt",ios::in);
	for (int i = 0; i < animB.vertices.size(); i++) {
		Vector3 tmpNormal;
		f >> tmpNormal[0] >> tmpNormal[1] >> tmpNormal[2];
		EXPECT_NEAR(0.0, (tmpNormal - animB.normals[i]).norm(), 1e-3);
	}
	f.close();
}

TEST(Animation2DTest, FindClosestPoint) {
	Animation2D anim;
	anim.load("tests/drop0/drop0",10, 12);

	Vector3 queryPoint;
	
	{
		queryPoint = Vector3(0.5,0.5,0);
		auto c = anim.findClosestPoint(queryPoint);	
		EXPECT_NEAR(0.0, (Vector3(0.49991639,  0.42535209,  0.) - c.p).norm(), tol);
		LOG(INFO)<<c.p<<c.n;
	}


	{
		queryPoint = Vector3(0.5,0.5, anim.getDepth(1));
		auto c = anim.findClosestPoint(queryPoint);	
		EXPECT_NEAR(0.0, (Vector3(0.49991639,  0.42535209,  0.) - c.p).norm(), tol);
		LOG(INFO)<<c.p<<c.n;
	}

	{
		queryPoint = Vector3(0.1,0.2, anim.getDepth(1));
		auto c = anim.findClosestPoint(queryPoint);	
		EXPECT_NEAR(0.0, (Vector3(0.317971  ,  0.409987  , -0.02308143) - c.p).norm(), tol);
		LOG(INFO)<<c.p<<c.n;
	}

	{
		queryPoint = Vector3(0.3,0.8, anim.getDepth(2));
		auto c = anim.findClosestPoint(queryPoint);	
		EXPECT_NEAR(0.0, (Vector3( 0.315625,  0.734376, -0.      ) - c.p).norm(), tol);
		LOG(INFO)<<c.p<<c.n;
	}

	{
		queryPoint = Vector3(0.8,0.4, anim.getDepth(2));
		auto c = anim.findClosestPoint(queryPoint);	
		EXPECT_NEAR(0.0, (Vector3(0.61183421,  0.42287672, -0.02308143) - c.p).norm(), tol);
		LOG(INFO)<<c.p<<c.n;
	}
}



TEST(Animation2DTest, Assignment) {
	Animation2D a;
	a.load("tests/drop0/drop0",10, 12);
	Animation2D b;	
	b = a;
	EXPECT_EQ(a.vertices.size(), b.vertices.size());
}

TEST(Animation2DTest, CopyConstructor) {
	Animation2D a;
	a.load("tests/drop0/drop0",10, 12);
	Animation2D b(a);	
	EXPECT_EQ(a.vertices.size(), b.vertices.size());
}