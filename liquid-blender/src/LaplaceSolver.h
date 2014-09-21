#include "stdafx.h"

struct LaplaceNeighbour {
	virtual int getIndex() = 0;
	virtual Float getWeight() = 0;
};

template<typename T> 
class LaplaceSolver {

public:
	LaplaceSolver(int _numVars, std::function<bool(int,T&)> _isConstrained, std::vector< std::vector<LaplaceNeighbour*> > &adjacencyList) {
		numOrigVars = _numVars;
		isConstrained = _isConstrained;

		numSolverVars=0;

		for (int i = 0; i < numOrigVars; i++) {
			T tmp;
			if (!isConstrained(i,tmp)) {
				solverToOrigMap[numSolverVars] = i;
				origToSolverMap[i] = numSolverVars;
				numSolverVars++;				
			}
		}

		VLOG(1) << "#Original Variables:" << numOrigVars;
		VLOG(1) << "#Solver Variables:" << numSolverVars;

		// Construct A and RHS
		A.resize(numSolverVars,numSolverVars);
		rhs.resize(numSolverVars);
		for (int i = 0; i < numSolverVars; i++) {
			// Populate each row
			int origIndex = solverToOrigMap[i];
			
			Float sum = 0.0;
			rhs[i] = 0.0;

			for (int j = 0; j < adjacencyList[origIndex].size(); j++) {
				T constrainedValue;
				LaplaceNeighbour *nb = adjacencyList[origIndex][j];
				sum += nb->getWeight();

				if (isConstrained(nb->getIndex(),constrainedValue)) {
					rhs[i] += nb->getWeight() * constrainedValue;
				}
				else {
					A.insert(i, origToSolverMap[nb->getIndex()]) = -nb->getWeight();					 
				}
			}
			A.insert(i,i) = sum;
		}
		
	}

	VectorX solve() {
		Eigen::ConjugateGradient<Eigen::SparseMatrix<Float> > cg;
		cg.compute(A);
		VectorX x = cg.solve(rhs);
		LOG(INFO) << "LaplaceSolver: #iterations:     " << cg.iterations() << std::endl;
		LOG(INFO) << "LaplaceSolver: estimated error: " << cg.error()      << std::endl;		

		VectorX sol(numOrigVars);

		for (int i = 0; i < numOrigVars; i++) {
			T constrainedValue;
			if (isConstrained(i,constrainedValue)) {
				sol[i] = constrainedValue;
			}
			else {
				sol[i] = x[origToSolverMap[i]];
			}
		}
		return sol;
	}


protected:

	std::map<int,int> origToSolverMap;
	std::map<int,int> solverToOrigMap;
	VectorX rhs;
	Eigen::SparseMatrix<Float> A;
	int numOrigVars;
	int numSolverVars;
	std::function<bool(int,T&)> isConstrained;
};