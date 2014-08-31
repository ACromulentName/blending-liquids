#include "stdafx.h"
#include "ICP.hpp"
#include "Animation.h"
#include "gtest/gtest.h"
#include "Utility.hpp"
#include "Configuration.hpp"
#include "Blend.hpp"

//TEST(Blend3DTest, BoxGravity) {
//	wstring corrFilename(L"tests/configs/boxGravity.corr");
//	blend<3> (corrFilename, 0);
//}


//TEST(Blend3DTest, BoxResampled) {
//	wstring corrFilename(L"tests/configs/boxResampled.corr");
//	blend<3> (corrFilename, 0);
//}

TEST(Blend3DTest, Drop0) {
	string corrFilename("tests/configs/drop0.corr");
	blend<3> (corrFilename, 0);
}