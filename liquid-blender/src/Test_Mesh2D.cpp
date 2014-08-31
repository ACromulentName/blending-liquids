#include "stdafx.h"
#include "Mesh2D.h"

#include "gtest/gtest.h"

TEST(Mesh2DTest, Load) {
	Mesh2D mesh;

	LOG(INFO) << "Running mesh load test";
	EXPECT_TRUE(mesh.load("tests/drop0/drop0_0000.mesh"));
	EXPECT_EQ(102, mesh.vertices.size());
	EXPECT_EQ(102, mesh.edges.size());

	EXPECT_TRUE(mesh.loadVelocities("tests/drop0/drop0_0000.vvel"));
	EXPECT_EQ(102, mesh.velocities.size());

	EXPECT_EQ(97, mesh.nextVertexIndex(0));
	EXPECT_EQ(95, mesh.prevVertexIndex(0));

	EXPECT_EQ(207, mesh.nextVertex(0).id);
	EXPECT_EQ(203, mesh.prevVertex(0).id);
}
