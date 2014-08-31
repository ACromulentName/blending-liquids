#include "stdafx.h"
#include "Animation.h"
#include "MarchingSquares2D.h"

#include "gtest/gtest.h"

TEST(MarchingSquaresTest, Basic) {
	Animation2D anim;
	anim.load("tests/dropPool0/dropPool0",0, 50);

	for (int i = 0; i < 50; i++) {
		SDF<2> sdf;
		sdf.computeSDF(anim, i, Globals::AvgEdgeLength);
		MarchingSquares<2> ms(sdf);
		char frameName[1024];
		sprintf_s(frameName, 1024, "test_%04d.ms", i);
		ms.extract(frameName);
	}
}