#include "stdafx.h"
#include "Mesh3D.h"
#include "Animation.h"
#include "Utility.hpp"
#include "gtest/gtest.h"

TEST(Animation3DTest, Load) {
	Animation3D anim;
	EXPECT_TRUE(anim.load("tests/basic3d",0,1));
}