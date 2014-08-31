#include "stdafx.h"
#include "Utility.hpp"
#include "Animation.h"
#include "gtest/gtest.h"

TEST(DeterminantTest, Identity) {
	Matrix3 mat;
	mat << 
		1, 0, 0,
		0, 1, 0,
		0, 0, 1;

	EXPECT_EQ(determinant(mat), 1.0);
}

TEST(DeterminantTest, Diag) {
	Matrix3 mat;
	mat << 
		2, 0, 0,
		0, 3, 0,
		0, 0, 4;

	EXPECT_EQ(determinant(mat), 24.0);
}

TEST(DeterminantTest, Skew) {
	Matrix3 mat;
	mat << 
		0, 2, -1,
		-2, 0, -4,
		1, 4, 0;

	EXPECT_EQ(determinant(mat), 0.0);
}