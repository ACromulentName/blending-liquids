#include "stdafx.h"
#include "ICP.hpp"
#include "Animation.h"
#include "gtest/gtest.h"
#include <cstdlib>
#include <fstream>
#include <array>
#include "Utility.hpp"
#include "Configuration.hpp"

typedef Eigen::Matrix<Float, 3, 3> Matrix3;
typedef ICP<2> ICP2D;

TEST(ICPTest, Constructor)  {
	Animation2D a;
	Animation2D b;

	a.load("tests/drop0/drop0",0, 10);
	b.load("tests/drop0/drop0",0, 10);

	PointCorrs pointCorrs;
	LineCorrs lineCorrs;
	ICP<2> icp(a,b, 10, pointCorrs, lineCorrs );
	
	

	EXPECT_EQ (0.0, icp.getWeights().weightRigid);
	EXPECT_EQ (1.0, icp.getWeights().weightSmooth);

	icp.initialize();
	
	array <int,22> expectedNodes = {
		0 ,
		43 ,
		63 ,
		77 ,
		87 ,
		94 ,
		56 ,
		36 ,
		22 ,
		12 ,
		2 ,
		1027 ,
		1070 ,
		1090 ,
		1104 ,
		1113 ,
		1129 ,
		1087 ,
		1067 ,
		1051 ,
		1041 ,
		1031
	};


	for (int i = 0; i < expectedNodes.size(); i++ ) {
		EXPECT_EQ(expectedNodes[i], icp.nodes[i]);
	}

	array<Float, 12> expectedSmoothnessWeights = {
		5 , 0.408248290464,
		4 , 0.408248290464,
		17, 0.408248290464,
		15, 0.408248290464,
		14, 0.408248290464,
		18, 0.408248290464
	};

	for (int i = 0; i < icp.nodeWeights[16].size(); i++) {
		for (int j = 0; j < 6; j++) {
			if (expectedSmoothnessWeights[2*j] == icp.nodeWeights[16][i].nodeIndex) {
				EXPECT_NEAR(expectedSmoothnessWeights[2*j+1], icp.nodeWeights[16][i].weight,1e-3);
			}
		}				
	}

}


TEST(ICPTest, DeformationTest) {
	Animation2D animA;
	Animation2D animB;

	animA.load("tests/drop0/drop0",0, 10);
	animB.load("tests/drop0/drop0",0, 10);

	PointCorrs pointCorrs;
	LineCorrs lineCorrs;
	ICP<2> icp(animA, animB, 10, pointCorrs, lineCorrs);


	icp.initialize();

	Animation2D animC;
	vector< Matrix3 > A;
	vector < Vector3 > b;

	for (int i = 0; i < icp.getNumNodes(); i++) {		
		Matrix3 tmpA;
		tmpA.setIdentity();
		Vector3 tmpB = Vector3::Zero();
		tmpB[2] = 1;

		A.push_back(tmpA);
		b.push_back(tmpB);

	}

	icp.deformAnimation(animC, A, b);


	for (int i = 0; i < animA.vertices.size(); i++) {
		EXPECT_NEAR(0.0, (animA.vertices[i] + Vector3(0,0,1) - animC.vertices[i]).norm(), 1e-6);
	}

}

TEST(ICPTest, DISABLED_Correspondences) {
	Animation2D animA;
	Animation2D animB;

	animA.load("tests/drop0/drop0",0, 10);
	animB.load("tests/drop1/drop1",0, 10);

	PointCorrs pointCorrs;
	LineCorrs lineCorrs;
	ICP2D icp(animA, animB, 16, pointCorrs, lineCorrs);
	icp.initialize();

	icp.computeCorrespondences();
	
	int count = 0;
	for (int i = 0; i < icp.corrs.size(); i++)
	{
		if (icp.corrs[i].valid)
			count++;
	}

	EXPECT_EQ(1131, count);

}



TEST(ICPTest, DISABLED_ESmooth) {
	Animation2D animA;
	Animation2D animB;

	animA.load("tests/drop0/drop0",0, 10);
	animB.load("tests/drop1/drop1",0, 10);

	PointCorrs pointCorrs;
	LineCorrs lineCorrs;
	ICP2D icp(animA, animB, 16, pointCorrs, lineCorrs);
	icp.initialize();
	icp.setupMatrices();

	for (int i = 0; i < icp.nodes.size(); i++) {
		for (int j = 0; j < icp.nodeWeights[i].size(); j++) {
			icp.computeESmooth(i, j);
		}		
	}
	
}

TEST(ICPTest, DISABLED_JSmooth) {
	Animation2D animA;
	Animation2D animB;

	animA.load("tests/drop0/drop0",0, 10);
	animB.load("tests/drop1/drop1",0, 10);

	PointCorrs pointCorrs;
	LineCorrs lineCorrs;
	ICP2D icp(animA, animB, 10, pointCorrs, lineCorrs);
	icp.initialize();
	icp.setupMatrices();
	
	vector< Matrix3 > A;
	vector< Vector3> b;
	icp.getNodes(A, b);

	for (int i = 0; i < icp.nodes.size(); i++) {
		for (int j = 0; j < 3; j++) {
			b[i][j] = (Float)rand()/(float)RAND_MAX;
		}
	}
	icp.setNodes(A,b);

	icp.curRow = 0;
	
	for (int i = 1; i < icp.nodes.size(); i++) {
		for (int j = 0; j < icp.nodeWeights[i].size(); j++) {				
			icp.computeJSmooth(i,j);
			Eigen::VectorXd jfd = icp.computeJSmoothFD(i,j);
			for (int k = 0; k < icp.numFreeVariables; k++) {				
				EXPECT_NEAR(0.0, jfd[k]-icp.jacobian.coeff(icp.curRow,k), 1e-5);			
			}
			icp.curRow++;
		}
	}
 }

TEST(ICPTest, JFitPlane) {
	Animation2D animA;
	Animation2D animB;

	animA.load("tests/drop0/drop0",0, 10);
	animB.load("tests/drop1/drop1",0, 10);

	PointCorrs pointCorrs;
	LineCorrs lineCorrs;
	ICP2D icp(animA, animB, 16, pointCorrs, lineCorrs);
	icp.initialize();


	vector< Matrix3 > A;
	vector< Vector3> b;
	icp.getNodes(A, b);

	// Note: Do NOT test near discontinuous point!
	for (int i = 0; i < icp.nodes.size(); i++) {
		for (int j = 0; j < 3; j++) {
			b[i][j] = 0.001;
		}
	}
	icp.setNodes(A,b);

	icp.curRow = 0;
	icp.deformAnimation(icp.curAnim, icp.A, icp.b);
	icp.computeCorrespondences();
	icp.setupMatrices();

	int count = 0;
	for (int i = 0; i < animA.vertices.size(); i++) {
		if (!icp.corrs[i].valid) 
			continue;

		icp.computeJFitPlane(i);

		cout << "Computing:" << i << " "<< icp.corrs[i].p << " " << icp.corrs[i].n << endl;
		Eigen::VectorXd jfd = icp.computeJFitPlaneFD(i);
		for (int k = 0; k < icp.numFreeVariables; k++) {	
			EXPECT_NEAR(jfd[k], icp.jacobian.coeff(icp.jFitPlaneMap[i],k), 1e-5);			
		}		
		
	}

 }

TEST(ICPTest, JFitPoint) {
	Animation2D animA;
	Animation2D animB;

	//animA.load("tests/drop0/drop0",0, 10);
	//animB.load("tests/drop1/drop1",0, 10);

	//PointCorrs pointCorrs;
	//LineCorrs lineCorrs;
	//ICP2D icp(animA, animB, 16, pointCorrs, lineCorrs);
	//icp.initialize();


	//vector< Matrix3 > A;
	//vector< Vector3> b;
	//icp.getNodes(A, b);

	//// Note: Do NOT test near discontinuous point!
	//for (int i = 0; i < icp.nodes.size(); i++) {
	//	for (int j = 0; j < 3; j++) {
	//		b[i][j] = 0.01;
	//	}
	//}
	//icp.setNodes(A,b);

	//icp.curRow = 0;
	//icp.deformAnimation(icp.curAnim, icp.A, icp.b);
	//icp.computeCorrespondences();
	//icp.setupMatrices();

	//int count = 0;
	//for (int i = 0; i < animA.vertices.size(); i++) {
	//	if (!icp.corrs[i].valid) 
	//		continue;

	//	icp.computeJFitPoint(i);

	//	//cout << "Computing:" << i << " "<< icp.corrs[i].p << " " << icp.corrs[i].n << endl;
	//	Eigen::VectorXd jfd = icp.computeJFitPointFD(i);
	//	for (int k = 0; k < icp.numFreeVariables; k++) {	
	//		EXPECT_NEAR(jfd[k], icp.jacobian.coeff(icp.jFitPointMap[i],k), 1e-5);			
	//	}				
	//}

}
//
//TEST(ICPTest, Drop) {
//	Animation2D animA;
//	Animation2D animB;
//
//	animA.load("tests/drop0/drop0",0, 10);
//	animB.load("tests/drop1/drop1",0, 10);
//
//	PointCorrs pointCorrs;
//	LineCorrs lineCorrs;
//	ICP2D icp(animA, animB, 16, pointCorrs, lineCorrs);
//	icp.initialize();
//
//	icp.doIteration();
//}
//
//TEST(ICPTest, LoadConfig0) {
//	Animation2D animA, animB;
//	ICPWeights initialWeights;
//	int samplingRate;
//	int maxIterations;
//	PointCorrs pointCorrs;
//	LineCorrs lineCorrs;
//	string configName;
//
//	EXPECT_TRUE(loadConfig("tests/configs/drop0.cfg", configName, animA, animB, initialWeights, maxIterations, samplingRate));
//
//	ICP2D icp(animA, animB, samplingRate, pointCorrs, lineCorrs);
//	icp.initialize();
//	icp.doIteration();
//}

//TEST(ICPTest, LoadConfig1) {
//	Configuration<2> cfg;
//
//	EXPECT_TRUE(loadConfig<2>("tests/configs/drop0.cfg", cfg));
//
//	ICP2D icp(cfg);
//	icp.initialize();
//	
//	int saveEvery = 5;
//	for (int i = 1; i <= cfg.maxIterations; i++) {
//		cout << "Iteration: " << i << endl;
//		icp.doIteration();		
//
//		if (i % saveEvery == 0)
//			saveRegistration<2> (cfg, icp, i);
//	}
//	
//}

TEST(ICPTest, LoadCorrs) {
	Configuration<2> cfg;

	EXPECT_TRUE(cfg.load("tests/configs/dropPoolLong-2-3.cfg"));
	cfg.maxIterations = 1;
	ICP2D icp(cfg);
	icp.initialize();

	int saveEvery = 5;
	for (int i = 1; i <= cfg.maxIterations; i++) {
		cout << "Iteration: " << i << endl;
		icp.doIteration();		

		if (i % saveEvery == 0)
			saveRegistration<2> (cfg, icp, i);
	}

}