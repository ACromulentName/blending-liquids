#pragma once

#include "stdafx.h"
#include "Mesh2D.h"
#include "Mesh3D.h"
#include "Neighbour.h"

#include <queue>
#include <set>
#include <fstream>
#include "boost/filesystem.hpp"
#include "tbb/atomic.h"
#include <functional>
#include "json.h"


using namespace std;

struct ICPWeights;
template<int DIM> class Animation;

bool readFileIntoString(std::string& filename, std::string& data);
std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems, bool removeEmptyStrings = false);
std::vector<std::string> split(const std::string &s, char delim, bool removeEmptyStrings = false);

Float determinant(const Matrix3& mat);

class VisitorMesh2D {
public:
	VisitorMesh2D(const Mesh2D& _mesh): mesh(_mesh) {
		curIndex = 0;
		for (int i = 0; i < mesh.vertices.size(); i++) {
			visited.push_back(false);
		}
	}

	int next()
	{
		if (!visited[curIndex]) {
			int tmpIndex = curIndex;
			visited[curIndex] = true;
			curIndex = mesh.nextVertexIndex(curIndex);
			return tmpIndex;

		}
		else {
			// Find new component
			curIndex = 0;
			while (visited[curIndex]) {
				curIndex++;
				if (curIndex == mesh.vertices.size()) {
					return -1;
				}
			};

			int tmpIndex = curIndex;
			visited[curIndex] = true;
			curIndex = mesh.nextVertexIndex(curIndex);
			return tmpIndex;
		}

		return -1;
	}

	~VisitorMesh2D() {
		visited.clear();
	}

private:
	int curIndex;
	vector<bool> visited;
	const Mesh2D& mesh;

};



// Performs Laplace smoothing on data. 
template<typename T, typename S> 
void smoothLaplace(vector<T>& data, int numSmoothingPasses, function<vector<S>(int)> getNeighbourIndices) {
	
	for (int i = 0; i < numSmoothingPasses; i++) {
		vector<T> tmpData = data;

		for (int j = 0; j < data.size(); j++) {
			vector<S> neighbourIndices = getNeighbourIndices(j);
			T sum;
			Float sumWeights = 0;
			for (int k = 0; k < neighbourIndices.size(); k++) {
				if (k==0)
					sum = neighbourIndices[k].getWeight() * data[neighbourIndices[k].getIndex()];
				else
					sum += neighbourIndices[k].getWeight() * data[neighbourIndices[k].getIndex()];

				sumWeights += neighbourIndices[k].getWeight();
			}			
			tmpData[j] = sum / sumWeights;
		}

		data = tmpData;
	}
}


// Does a breadth first search starting at point given by startIndex and returns all a map of vertices and their distances (within maxRadius). 
template<int DIM> 
std::vector<Neighbour> doBFS(const unique_ptr<Animation<DIM>>& anim, int startIndex, Float maxRadius) {

	map<int, Float> curNeighbourMap;
	curNeighbourMap.clear();

	vector<Neighbour> neighbours;
	priority_queue<Float, vector<Neighbour>, std::greater<Neighbour> > q;

	unordered_set<int> visited;

	curNeighbourMap[startIndex] = 0.0;

	q.push(Neighbour(startIndex, 0.0));

	int total = 0;

	while (!q.empty()) {

		// Pop top element
		Neighbour n = q.top();
		q.pop();

		total++;

		// Termination condition
		if (n.d > maxRadius) 
			break;
		
		// Visit neighbours
		if (visited.find(n.index) == visited.end()) {
			// Mark as visited
			visited.insert(n.index);
			for (int i = 0; i < anim->neighbours[n.index].size(); i++) {
				if (anim->neighbours[n.index][i].index == Globals::InvalidNeighbour) 
					continue;

				int nbIndex = anim->neighbours[n.index][i].index;
				Float nbDist = anim->neighbours[n.index][i].d;

				auto nb = curNeighbourMap.find(nbIndex);

				if (curNeighbourMap.find(nbIndex) != curNeighbourMap.end()) {
					if (n.d + nbDist < curNeighbourMap[nbIndex]) {
						curNeighbourMap[nbIndex] = n.d + nbDist;
						
						if (visited.find(nbIndex) == visited.end())
						{				
							q.push(Neighbour(nbIndex, n.d + nbDist));
							//LOG(INFO) << "Added:" << nbIndex << " " << n.d + nbDist;
						}
					}
				}
				else {
					curNeighbourMap[nbIndex] = n.d + nbDist;
					if (visited.find(nbIndex) == visited.end())
					{				
						q.push(Neighbour(nbIndex, n.d + nbDist));
						//LOG(INFO) << "Added:" << nbIndex << " " << n.d + nbDist;
					}
				}	

			}
		}
		
	
	};
	
	
	for (auto it = curNeighbourMap.begin(); it != curNeighbourMap.end(); it++) {		
		neighbours.push_back(Neighbour(it->first, it->second));
	}

	return neighbours;
}



// Does a component-wise breadth first search starting at point given by startIndex and returns all a map of vertices and their distances. 
template<int DIM> 
std::vector<AnisotropicNeighbour> doSplitBFS(const unique_ptr<Animation<DIM>>& anim, vector<int> startIndices, Float rTime, Float rSpace) {

	Float maxRadius = max(rTime, rSpace);
	map<int, AnisotropicNeighbour > curNeighbourMap;
	curNeighbourMap.clear();

	vector<AnisotropicNeighbour> neighbours;
	priority_queue<Float, vector<AnisotropicNeighbour>, std::greater<AnisotropicNeighbour> > q;

	unordered_set<int> visited;



	for (auto index: startIndices) {
		auto nb = AnisotropicNeighbour(index, 0.0, 0.0, 0.0);
		q.push(nb);
		curNeighbourMap[index] = nb;
	}

	// Now trace forwards and backwards in time and add these to the starting set 
	for (auto index: startIndices) {
		int curIndex = index;

		int numTimesteps = 1;
		curIndex = anim->neighbours[index][PrevTime].index;
		Float dTime = abs(numTimesteps*anim->timeScaling);
		while (curIndex != Globals::InvalidNeighbour &&  dTime < rTime ) {			
			AnisotropicNeighbour nb = AnisotropicNeighbour(curIndex, dTime, dTime, 0.0);
			q.push(nb);

			curNeighbourMap[index] = nb;
			curIndex = anim->neighbours[curIndex][PrevTime].index;
			numTimesteps++;			
			dTime = abs(numTimesteps*anim->timeScaling);
		}

		numTimesteps = 1;

		curIndex = anim->neighbours[index][NextTime].index;
		dTime = abs(numTimesteps*anim->timeScaling);
		while (curIndex != Globals::InvalidNeighbour &&  dTime < rTime ) {			
			AnisotropicNeighbour nb = AnisotropicNeighbour(curIndex, dTime, dTime, 0.0);
			q.push(nb);

			curNeighbourMap[index] = nb;
			curIndex = anim->neighbours[curIndex][NextTime].index;
			numTimesteps++;			
			dTime = abs(numTimesteps*anim->timeScaling);
		}

	}

	int total = 0;

	while (!q.empty()) {

		// Pop top element
		AnisotropicNeighbour n = q.top();
		q.pop();

		total++;

		// Termination condition
		if (n.d > maxRadius) 
			break;

		// Visit neighbours
		if (visited.find(n.index) == visited.end()) {
			// Mark as visited
			visited.insert(n.index);


			vector<AnisotropicNeighbour> curNbs;
			// Existing neighbours
			for (int i = 0; i < anim->neighbours[n.index].size(); i++) {
				int nbIndex = anim->neighbours[n.index][i].index;
				Float dSq = (anim->vertices[nbIndex] - anim->vertices[n.index]).squaredNorm();
				Float dTime = abs(anim->vertices[nbIndex][DIM] - anim->vertices[n.index][DIM]);
				Float dSpace = sqrt(dSq - dTime*dTime);
				curNbs.push_back(AnisotropicNeighbour(nbIndex, sqrt(dSq), dTime, dSpace));
			}

			// We need to consider diagonal edges
			for (int timeIndex = 0; timeIndex < 2; timeIndex++) {

				if (anim->neighbours[n.index][timeIndex].index != Globals::InvalidNeighbour) {
					int timeNbIndex = anim->neighbours[n.index][timeIndex].index;
					// Get all spatial neighbours
					for (int j = 2; j < anim->neighbours[timeNbIndex].size(); j++) {
						Float dSq = (anim->vertices[ anim->neighbours[timeNbIndex][j].index ] - anim->vertices[ n.index ]).squaredNorm();
						Float dTime = abs(anim->vertices[ anim->neighbours[timeNbIndex][j].index ][DIM] - anim->vertices[ n.index ][DIM]);
						Float dSpace = sqrt(dSq - dTime*dTime); 
											
						curNbs.push_back(AnisotropicNeighbour(anim->neighbours[timeNbIndex][j].index, sqrt(dSq), dTime, dSpace));
					}
				
				}
			}


			for (int i = 0; i < curNbs.size(); i++) {
				if (curNbs[i].index == Globals::InvalidNeighbour) 
					continue;

				int nbIndex = curNbs[i].index;
				Float nbDist = curNbs[i].d;
				Float nbDistTime = curNbs[i].dTime;
				Float nbDistSpace = curNbs[i].dSpace;

				auto nb = curNeighbourMap.find(nbIndex);

				if (curNeighbourMap.find(nbIndex) != curNeighbourMap.end()) {
					if (n.d + nbDist < curNeighbourMap[nbIndex].d) {
						AnisotropicNeighbour newNb = AnisotropicNeighbour(nbIndex, n.d + nbDist, n.dTime + nbDistTime, n.dSpace + nbDistSpace);
						curNeighbourMap[nbIndex] = newNb;

						if (visited.find(nbIndex) == visited.end())
						{				
							q.push(newNb);
							//LOG(INFO) << "Added:" << nbIndex << " " << n.d + nbDist;
						}
					}
				}
				else {
					AnisotropicNeighbour newNb = AnisotropicNeighbour(nbIndex, n.d + nbDist, n.dTime + nbDistTime, n.dSpace + nbDistSpace);
					curNeighbourMap[nbIndex] = newNb;
					if (visited.find(nbIndex) == visited.end())
					{				
						q.push(newNb);
						//LOG(INFO) << "Added:" << nbIndex << " " << n.d + nbDist;
					}
				}	

			}
		}


	};

	for (auto it = curNeighbourMap.begin(); it != curNeighbourMap.end(); it++) {		
		neighbours.push_back(it->second);
	}
	return neighbours;
}