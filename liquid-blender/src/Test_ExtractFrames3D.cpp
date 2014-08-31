#include "stdafx.h"
#include "ICP.hpp"
#include "Animation.h"
#include "gtest/gtest.h"
#include "Utility.hpp"
#include "Configuration.hpp"

TEST(ExtractFrames3DTest, BoxResampled)  {
	Animation3D a;
	Animation3D b;

	a.load("tests/boxResampled/boxResampled",0, 10);
	b.load("tests/boxResampled/boxResampled",0, 10);
	int maxIterations = 10;
	int saveEvery = 10;
	PointCorrs pointCorrs;
	LineCorrs lineCorrs;
	ICP<3> icp(a,b, 10, 1, maxIterations, saveEvery, pointCorrs, lineCorrs );

	icp.initialize();

	Animation3D interpolatedAnim = a;
	interpolatedAnim.extractSlice("boxResampled", interpolatedAnim.getDepth(1), true);
}

//TEST(ExtractFrames3DTest, Box)  {
//	Animation3D a;
//	Animation3D b;
//
//	a.load("tests/boxGravity/boxGravity",0, 10);
//	b.load("tests/boxGravity/boxGravity",0, 10);
//
//	PointCorrs pointCorrs;
//	LineCorrs lineCorrs;
//	ICP<3> icp(a,b, 10, pointCorrs, lineCorrs );
//
//	icp.initialize();
//
//	Animation3D interpolatedAnim = a;
//	interpolatedAnim.extractSlice("boxGravity", interpolatedAnim.getDepth(0.99));
//}

//TEST(ExtractFrames3DTest, Box)  {
//	Animation3D a;
//	Animation3D b;
//
//	a.load("tests/box/box",0, 10);
//	b.load("tests/box/box",0, 10);
//
//	PointCorrs pointCorrs;
//	LineCorrs lineCorrs;
//	ICP<3> icp(a,b, 10, pointCorrs, lineCorrs );
//
//	icp.initialize();
//
//	Animation3D interpolatedAnim = a;
//	interpolatedAnim.extractSlice("box", interpolatedAnim.getDepth(0.1));
//}


//TEST(ExtractFrames3DTest, Ball)  {
//	Animation3D a;
//	Animation3D b;
//
//	a.load("tests/dropGravity/dropGravity",0, 30);
//	b.load("tests/dropGravity/dropGravity",0, 30);
//
//	PointCorrs pointCorrs;
//	LineCorrs lineCorrs;
//	ICP<3> icp(a,b, 10, pointCorrs, lineCorrs );
//
//	icp.initialize();
//
//	Animation3D interpolatedAnim = a;
//	interpolatedAnim.extractSlice("dropGravity", interpolatedAnim.getDepth(0.4));
//}