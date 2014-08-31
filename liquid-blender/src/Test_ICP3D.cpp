#include "stdafx.h"
#include "ICP.hpp"
#include "Animation.h"
#include "gtest/gtest.h"
#include "Utility.hpp"
#include "Configuration.hpp"


TEST(ICP3DTest, Constructor)  {
	Animation3D a;
	Animation3D b;

	

	a.load("tests/basic3d/basic3d",0, 10);
	b.load("tests/basic3d/basic3d",0, 10);

	PointCorrs pointCorrs;
	LineCorrs lineCorrs;
	int maxIterations = 10;
	int saveEvery = 10;
	ICP<3> icp(a,b, 10, 1,maxIterations, saveEvery, pointCorrs, lineCorrs );

	icp.initialize();
}