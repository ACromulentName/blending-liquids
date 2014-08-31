#include "stdafx.h"
#include "Utility.hpp"
#include "Animation.h"
#include "gtest/gtest.h"

using namespace std;

TEST(BFSTest, 2D) {
	Animation2D a;	
	a.load("tests/drop0/drop0",0, 10);

	int startIndex = 0;
	vector<Neighbour> neighbours = doBFS<2>(a, startIndex, 3 * a.timeScaling);


	float expectedPairs[] = {
		0 , 0 ,
		97 , 0.00796526013386 ,
		2 , 0.0159344360234 ,
		3 , 0.0284344360234 ,
		197 , 0.0203372972116 ,
		102 , 0.0123706675604 ,
		199 , 0.0203359276942 ,
		104 , 0.0283045156095 ,
		129 , 0.0283051101492 ,
		298 , 0.0329115881217 ,
		203 , 0.0249453064773 ,
		300 , 0.0329074367083 ,
		27 , 0.0159344425888 ,
		29 , 0.0284344425888 ,
		95 , 0.00796682577944 
	};
	
	int numExpectedPairs = 15;
		
	/*cout << "Neighbours of: "<< startIndex << endl;*/
	for (auto n : neighbours) {
		//cout<< a.globalToFrameMap[n.index] << ":" << n.index <<" " << n.d << endl;
		bool found = false;
		for (int i = 0; i < numExpectedPairs; i++) {
			if (expectedPairs[2*i] == n.index) {
				EXPECT_NEAR(expectedPairs[2*i+1], n.d, 1e-3);
				found = true;
				break;
			}			
		}
		EXPECT_TRUE(found);

	}	
}