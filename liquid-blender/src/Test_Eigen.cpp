#include "stdafx.h"
#include "Utility.hpp"
#include "gtest/gtest.h"

using namespace std;

struct A {
	Eigen::Vector2d v;
};

TEST(EigenTest, Output) {
	Vector4 v(1,2,3,4);

	
	LOG(INFO) << v.format(CommaFmt);
}

TEST(EigenTest, Alignment) {
	A *a = new A;
	a->v = Eigen::Vector2d(0,0);
	std::cout << a->v << endl;
	delete a;
}


TEST(EigenTest, AlignmentVector) {
	vector<A> a;
	a.resize(100);
	for (int i = 0 ; i < 10; i++) {
		A tmp;
		a.push_back(tmp);
	}
	cout << a.size() << endl;
}