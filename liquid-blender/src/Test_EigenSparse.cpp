#include "stdafx.h"
#include "gtest/gtest.h"

#include <Eigen/SparseCholesky>
#include <Eigen/PardisoSupport>
#include <Eigen/IterativeLinearSolvers>

#include <Eigen/OrderingMethods>
#include <Eigen/SparseQR>

using namespace Eigen;
TEST(EigenSparseTest, SparseMatrix) {
	int rows = 3;
	int cols = 3;
	SparseMatrix<Float> mat(rows, cols);
	
	mat.reserve(VectorXi::Constant(rows, 1));

	mat.insert(0,0) = 1;
	mat.insert(1,1) = 2;
	mat.insert(2,2) = 3;
	mat.makeCompressed();
	
	VectorX rhs(3);
	VectorX x(3);

	rhs[0] = 1;
	rhs[1] = 2;
	rhs[2] = 3;

	//SimplicialLDLT< SparseMatrix<Float> >  solver;	
	//PardisoLDLT< SparseMatrix<Float> >  solver;

	//solver.compute(mat);	
	//x = solver.solve(rhs);

	
	// CG Solve
	{
		ConjugateGradient< SparseMatrix<Float> >  solver;
		VectorX x0(3);
		x0[0] = x0[1] = x0[2] = 0.0;
		solver.compute(mat);	
		solver.setTolerance(1e-6);
		solver.setMaxIterations(100);
		x = solver.solveWithGuess(rhs,x0);
		std::cout << "#iterations: " << solver.iterations() << std::endl;
		std::cout << "estimated error: " << solver.error() << std::endl;	

	}

	{
		Eigen::SparseQR< SparseMatrix<Float>, Eigen::COLAMDOrdering<int>> solver (mat);
		x = solver.solve(rhs);

	}


	EXPECT_EQ(1, x[0]);
	EXPECT_EQ(1, x[1]);
	EXPECT_EQ(1, x[2]);
	
}

TEST(EigenSparseTest, MatrixOrder) {
	Matrix<int, 3, 3, RowMajor> A;
	A << 1,2,3, 4,5,6, 7,8, 9;
	//std::cout << A << std::endl;


	Matrix<int, 3, 3, ColMajor> B;
	B << 1,2,3, 4,5,6, 7,8, 9;
	//std::cout << B << std::endl;

	Matrix<int, 3, 3, RowMajor> C;
	C << 1,2,3, 4,5,6, 7,8, 9;

	//std::cout << A * C << std::endl;


	//std::cout << A * B << std::endl;

	

}