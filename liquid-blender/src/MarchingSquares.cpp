#include "stdafx.h"
#include "MarchingSquares2D.h"


template<>
void SDF<2>::computeSDF(const Animation<2>& anim, int time, Float edgeLength) {
	edgeLength = edgeLength;
	numPointsX = int(1.0 / edgeLength) + 1;

	int totalNumPoints = 1;

	for (int i = 0; i < 2; i++) {
		totalNumPoints *= numPointsX;
	}

	field.clear();
	field.reserve(totalNumPoints);

	for (int i = 0; i < (totalNumPoints); i++) {
		field.push_back(0.0);
	}

	for (int i = 0; i < totalNumPoints; i++) {
		int y = i / (numPointsX);
		int x = i % numPointsX;

		Vector3 point = Vector3 ( edgeLength * x, edgeLength * y, anim.getDepth(time));
		field[i] = anim.getSignedDistance(point);
	}
}


template<>
void MarchingSquares<2>::extract(string filename) {

	struct Edge2D {
		Edge2D(Vector2 start, Vector2 end): start(start), end(end) {}
		Vector2 start;
		Vector2 end;
	};


	vector < vector<Edge2D> > edges;
	int numCellsX = sdf.numPointsX - 1;

	for (int i = 0; i < numCellsX * numCellsX; i++)  {
		edges.push_back( vector<Edge2D>() );
	}
		
	for (int i = 0 ; i < numCellsX * numCellsX; i++) {
		int cellType = 0;
		int x = i % numCellsX;
		int y = i / numCellsX;

		cellType |= sdf.isInside(x,y+1);
		cellType = cellType << 1;
		cellType |= sdf.isInside(x+1,y+1);
		cellType = cellType << 1;
		cellType |= sdf.isInside(x+1,y);
		cellType = cellType << 1;
		cellType |= sdf.isInside(x,y);

		if (cellType == 0 || cellType == 15) {
			// Skip (empty cell)
		}
		else if (cellType == 1 || cellType == 14) {
			edges[i].push_back( Edge2D(computeIntersection(x,y, x+1,y), computeIntersection(x,y, x,y+1)) );
		}
		else if (cellType == 2 || cellType == 13) {
			edges[i].push_back( Edge2D(computeIntersection(x,y,x+1,y), computeIntersection(x+1,y,x+1,y+1)) );
		}
		else if (cellType == 3 || cellType == 12) {
			edges[i].push_back( Edge2D(computeIntersection(x+1,y,x+1,y+1), computeIntersection(x,y+1,x,y)) );
		}
		else if (cellType == 4 || cellType == 11) {
			edges[i].push_back( Edge2D(computeIntersection(x+1,y,x+1,y+1), computeIntersection(x+1,y+1,x,y+1)) );
		}
		else if (cellType == 5) {
			edges[i].push_back( Edge2D(computeIntersection(x+1,y+1,x,y+1), computeIntersection(x,y+1,x,y)) );
			edges[i].push_back( Edge2D(computeIntersection(x,y,x+1,y), computeIntersection(x+1,y,x+1,y+1)) );
		}
		else if (cellType == 6 || cellType == 9) {
			edges[i].push_back( Edge2D(computeIntersection(x,y,x+1,y), computeIntersection(x+1,y+1,x,y+1)) );
		}
		else if (cellType == 7 || cellType == 8) {
			edges[i].push_back( Edge2D(computeIntersection(x+1,y+1,x,y+1), computeIntersection(x,y+1,x,y)) );
		}
	}

	int numEdges = 0;
	for (int i = 0; i < numCellsX * numCellsX; i++) {		
		numEdges+=edges[i].size();		
	}

	fstream f(filename, ios::out);
	try {
		f << numEdges << endl;
		for (int i = 0; i < numCellsX * numCellsX; i++) {		
			for (int j = 0; j < edges[i].size(); j++) {
				f << edges[i][j].start[0] << " " << edges[i][j].start[1] << " " << edges[i][j].end[0] << " " << edges[i][j].end[1] << endl;
			}
		}		
		f.close();
	}
	catch (exception e) {
		throw;
	}
	
}