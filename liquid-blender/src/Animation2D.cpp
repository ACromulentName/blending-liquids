#include "stdafx.h"
#include "Animation.h"
#include "Geometry.h"
#include <algorithm> 
#include "EdgeLabels.h"
using namespace Geometry;


/////////////////////////////////////////////////////////////////////////////
//  Specialization for Animation2D
/////////////////////////////////////////////////////////////////////////////

bool Animation2D::loadMeshes(const string path, const int startFrame, const int endFrame) {
	bool failed = false;

	// Create meshes and velocities for each frame
	for (int i = startFrame; i <= endFrame; i++) {
		
		frames.emplace_back(new Mesh2D());
	}

	tbb::parallel_for(startFrame, endFrame+1, 1, [&](int i) 
		//for (int i = startFrame; i < endFrame+1; i++)
	{
		char meshFilename[1024];
		sprintf_s(meshFilename, 1024, "%s_%04d.mesh", path.c_str(), i);

		if (!dynamic_pointer_cast<Mesh2D>(frames[i-startFrame])->load(meshFilename)) {				
			failed = true;
		}			

		char velFilename[1024];
		sprintf_s(velFilename, 1024, "%s_%04d.vvel", path.c_str(), i);

		if (!frames[i-startFrame]->loadVelocities(velFilename)) {				
			failed = true;
		}	
	}
	);

	return !failed;
}


// Compute time scaling i.e. average vertex spacing using the first frame 
void Animation2D::computeTimeScaling()
{

	timeScaling = 0.0;

	const Mesh2D* frame0 = dynamic_pointer_cast<const Mesh2D>(frames[0]).get();
	for (int i = 0; i < frame0->vertices.size(); i++) {
		timeScaling += (frame0->vertices[i].pos - frame0->vertices[frame0->nextVertexIndex(i)].pos).norm();
	}
	timeScaling /= frame0->vertices.size();	
	Globals::AvgEdgeLength = timeScaling;
	VLOG(1) << "Average edge length: " << timeScaling ;
}




// Finds the four immediate neighbours for every vertex
void Animation2D::findNeighbours()
{	
	const Float maxSqrDistance = 9 * timeScaling * timeScaling;
	// Create entries
	neighbours.clear();
	for (int i = 0; i < vertices.size(); i++) {		
		vector<Neighbour> nbs;	
		for (int i = 0; i < 4; i++){
			nbs.push_back(Neighbour());
		}
		neighbours.push_back(nbs);
	}

	
	

#ifdef ENABLE_PARALLEL
	tbb::parallel_for(size_t(0), vertices.size(), size_t(1), [&](int i)
#else
	for (int i = 0; i < vertices.size(); i++) 
#endif
	{
		FrameIndex frameIndex = globalToFrameMap[i];
		int timeIndex = frameIndex.first;
		int spaceIndex = frameIndex.second;
		int curId = frames[timeIndex]->vertices[spaceIndex].id;

		// Neighbours prevSpace and nextSpace should always exist since the meshes are closed
		neighbours[i][PrevSpace].index = frameToGlobalMap[ FrameIndex(timeIndex, dynamic_pointer_cast<const Mesh2D>(frames[timeIndex])->prevVertexIndex(spaceIndex)) ];
		neighbours[i][NextSpace].index = frameToGlobalMap[ FrameIndex(timeIndex, dynamic_pointer_cast<const Mesh2D>(frames[timeIndex])->nextVertexIndex(spaceIndex)) ];

		// We have to search for prevTime and nextTime
		Vector2 curVel = frames[timeIndex]->velocities[spaceIndex];
		Vector3 nextPos = Vector3( vertices[i][0] + dt * curVel[0], vertices[i][1] + dt * curVel[1], vertices[i][2] + getDepth(1));
		Vector3 prevPos = Vector3( vertices[i][0] - dt * curVel[0], vertices[i][1] - dt * curVel[1], vertices[i][2] - getDepth(1));

		// Note: If the same vertex id exists in the previous/next frame, automatically assign that vertex as the neighbours. Otherwise use the kdtree to search for the closest point

		// The first frame does not have any prevTime neighbours
		if (timeIndex > 0) {
			auto iter = frames[timeIndex-1]->idToIndexMap.find(curId);
			if (iter != frames[timeIndex-1]->idToIndexMap.end()) {
				neighbours[i][PrevTime].index = frameToGlobalMap[ FrameIndex(timeIndex - 1, iter->second) ];
			}
			else {
				int closestVertexIndex = findClosestVertex(prevPos);
				Float dSqr = (vertices[closestVertexIndex] - vertices[i]).squaredNorm();

				if (dSqr < maxSqrDistance && globalToFrameMap[closestVertexIndex].first == timeIndex - 1)
				{
					neighbours[i][PrevTime].index = closestVertexIndex;
				}
			}
		}

		// The last frame does not have nextTime neighbours
		if (timeIndex < frames.size() - 1) {

			auto iter = frames[timeIndex + 1]->idToIndexMap.find(curId);
			if (iter != frames[timeIndex + 1]->idToIndexMap.end()) {
				neighbours[i][NextTime].index = frameToGlobalMap[ FrameIndex(timeIndex + 1, iter->second) ];
			}
			else {
				int closestVertexIndex = findClosestVertex(nextPos);
				Float dSqr = (vertices[closestVertexIndex] - vertices[i]).squaredNorm();

				if (dSqr < maxSqrDistance && globalToFrameMap[closestVertexIndex].first == timeIndex + 1)
				{
					neighbours[i][NextTime].index = closestVertexIndex;
				}
			}
		}

		for (int j = 0; j < 4; j++) {
			if (neighbours[i][j].index != Globals::InvalidNeighbour) {
				neighbours[i][j].d = (vertices[i] - vertices[neighbours[i][j].index]).norm();
			}
		}
	}
#ifdef ENABLE_PARALLEL
	);
#endif
	

	// Debug
	//for (int i = 0; i < vertices.size(); i++) {
	//	FrameIndex frameIndex = globalToFrameMap[i];
	//	VLOG(1) << " Vertex : " << frameIndex << " has neighbours: ";
	//	for (int j = 0; j < 4; j++) {			
	//		VLOG(1) << neighbours[i][j];			
	//	}
	//}
}


void Animation2D::computeNormals()
{
	normals.clear();
	for (int i = 0; i < vertices.size(); i++)
	{
		normals.push_back(Vector3(0,0,0));
	}

#ifdef ENABLE_PARALLEL
	tbb::parallel_for(size_t(0), vertices.size(), size_t(1), [&](int i)
#else
	for (int i = 0; i < vertices.size(); i++)
#endif

	{
		FrameIndex fIndex = globalToFrameMap[i];

		normals[i] = computeNormalAtVertex(i);
	}

#ifdef ENABLE_PARALLEL
	);
#endif

	//for (int i = 0; i < vertices.size(); i++)
	//{
	//	DLOG(INFO)<<i<<" "<<normals[i][0]<<" "<<normals[i][1]<<" "<<normals[i][2]; 
	//}
}

/*
Returns averages normal at vertex given a global index.
*/
Vector3 Animation2D::computeNormalAtVertex(int index)
{
	int adjVertices[4][2] = { 
		{NextSpace, NextTime},
		{NextTime, PrevSpace},
		{PrevSpace, PrevTime},
		{PrevTime, NextSpace}
	};

	const Vector3& v = vertices[index];
	Vector3 normal(0.0,0.0,0.0);

	Float sum = 0.0;
	for (int i = 0; i < 4; i++) {
		if (neighbours[index][ adjVertices[i][0] ].index != Globals::InvalidNeighbour && 
			neighbours[index][ adjVertices[i][1] ].index != Globals::InvalidNeighbour) {
			const Vector3& v1 = vertices[ neighbours[index][adjVertices[i][0]].index ];
			const Vector3& v2 = vertices[ neighbours[index][adjVertices[i][1]].index ];

			Vector3 n = -(v1-v).cross(v2-v);
			n.normalize();

			sum += 1;
			normal += n;
		}
	}

	normal /= sum;
	normal.normalize();
	return normal;
}

/*
Finds the closest point and normal to a given point (not necessarily vertex)
*/
Correspondence<2> Animation2D::findClosestPoint(const VectorN& queryPoint) const{
	int closestVertexIndex = findClosestVertex(queryPoint);
	return findClosestPointAroundVertex(queryPoint, closestVertexIndex);
}

/*
Finds the closest point and normal to a given query point given that the closest vertex is known 
*/
Correspondence<2> Animation2D::findClosestPointAroundVertex(const Vector3& queryPoint, int closestVertexIndex) const {
	/*
	We need to run a set of tests here to avoid weird cases: 
	1) point to point
	2) point to edge (project onto each edge)
	3) point to face (project onto each face)
	*/

	Correspondence<2> corr;
	corr.valid = true;

	// Point to point
	Float closestDist = (vertices[closestVertexIndex] - queryPoint).norm();
	Vector3 closestPoint = vertices[closestVertexIndex];
	Vector3 closestNormal = normals[closestVertexIndex];
	const Vector3& v0 = vertices[closestVertexIndex];

	int adjVertices[4] = {PrevSpace, NextSpace, PrevTime, NextTime};

	int triangles[4][2] = { 
		{NextTime, NextSpace},
		{PrevSpace, NextTime},
		{PrevTime, PrevSpace},
		{NextSpace, PrevTime}
	};

	// Project onto edges
	for (int i = 0; i < 4; i++) {
		if (neighbours[closestVertexIndex][adjVertices[i]].index != Globals::InvalidNeighbour) {
			int nbIndex = neighbours[closestVertexIndex][adjVertices[i]].index;
			Vector3 projPt;
			Float interpolant;

			if (projectOntoEdge(queryPoint, vertices[closestVertexIndex], vertices[nbIndex], projPt, interpolant)) {
				Float d = (projPt - queryPoint).norm();
				if (d < closestDist) {
					closestDist = d;
					closestPoint = projPt;					
					closestNormal = Geometry::nlerp(normals[closestVertexIndex], normals[nbIndex], interpolant);
				}
			}
		}		
	}
	
	// Project onto faces
	for (int i = 0; i < 4; i++) {
		if (neighbours[closestVertexIndex][ triangles[i][0] ].index != Globals::InvalidNeighbour && 
			neighbours[closestVertexIndex][ triangles[i][1] ].index != Globals::InvalidNeighbour) {
				const Vector3& v1 = vertices[ neighbours[closestVertexIndex][triangles[i][0]].index ];
				const Vector3& v2 = vertices[ neighbours[closestVertexIndex][triangles[i][1]].index ];

				Vector3 projPt;
				Vector3 projNormal;
				if (projectOntoFace(queryPoint, v0, v1, v2, projPt, projNormal)) {
					Float d = (projPt - queryPoint).norm();
					if (d < closestDist) {
						closestDist = d;
						closestPoint = projPt;					
						closestNormal = projNormal;
				}					
			}
		}
	}

	corr.p = closestPoint;
	corr.n = closestNormal;
	return corr;
}

// Nothing required for 2D animations
void Animation2D::setupExtraction() {

}

void Animation2D::extractSlice(string outFilename, Float time, bool perturbVertices) {
	struct Edge2D {
		Edge2D(Vector3 start, Vector3 end): start(start), end(end) {}
		Vector3 start;
		Vector3 end;
	};
	//this->normals.clear();
	//for (int i = 0; i < this->vertices.size(); i++)
	//{
	//	this->normals.push_back(Vector3(0,0,0));
	//}

	vector<Edge2D>  edges;
	
	fstream errorFile(outFilename+".error", ios::out);

	for (int i = 0; i < vertices.size(); i++) {
		int numIntersections = 0;

		Vector3 intersections[4];
		bool valid[4];

		for (int j = 0; j < 4; j++)
			valid[j] = false;

		int nextScurT = neighbours[i][NextSpace].index;
		int nextSnextT = neighbours[nextScurT][Neighbours::NextTime].index;
		int curSnextT = neighbours[i][Neighbours::NextTime].index;

		if (curSnextT == Globals::InvalidNeighbour) {
			if (nextSnextT != Globals::InvalidNeighbour)
				curSnextT = neighbours[nextSnextT][PrevSpace].index;
		}

		int nextTnextS = -1;

		// Try 4 edges
		// (t, s) -> (t + 1, s)
		if (curSnextT != Globals::InvalidNeighbour) {
			Vector3 a = vertices[i];
			Vector3 b = vertices[curSnextT];

			valid[0] = Geometry::doesPlaneIntersectEdge(time, a, b, intersections[0]);
			if (valid[0])
				numIntersections++;

		}
		if (nextScurT != Globals::InvalidNeighbour) {
			if (nextSnextT != Globals::InvalidNeighbour) {
				// (t, s+1) -> (t + 1, s + 1)
				Vector3 a = vertices[nextScurT];
				Vector3 b = vertices[nextSnextT];

				valid[1] = Geometry::doesPlaneIntersectEdge(time, a, b, intersections[numIntersections]);
				if (valid[1]) 
					numIntersections++;
			}
			else
			{
				if (curSnextT != Globals::InvalidNeighbour) {
					nextTnextS = neighbours[curSnextT][NextSpace].index;
					if (nextTnextS != Globals::InvalidNeighbour) {
						// (t, s+1) -> (t + 1, s + 1)
						Vector3 a = vertices[nextScurT];
						Vector3 b = vertices[nextTnextS];
						valid[1] = Geometry::doesPlaneIntersectEdge(time, a, b, intersections[numIntersections]);
						if (valid[1])
							numIntersections++;
					}
				}
			}

			// (t, s) -> (t, s+1)
			Vector3 a = vertices[i];
			Vector3 b = vertices[nextScurT];

			valid[2] = Geometry::doesPlaneIntersectEdge(time, a, b, intersections[numIntersections]);
			if (valid[2]) 
				numIntersections++;
		}

		

		// Try : (t, s+1) -> ( s+1, t+1 )
		if (curSnextT != Globals::InvalidNeighbour && nextSnextT != Globals::InvalidNeighbour)
		{
			Vector3 a = vertices[curSnextT];
			Vector3 b = vertices[nextSnextT];
				valid[3] = Geometry::doesPlaneIntersectEdge(time, a, b, intersections[numIntersections]);
			if (valid[3])
				numIntersections++;

			
		}
		// Try : (t+1, s+1) -> (t, s+1) 
		else if (curSnextT != Globals::InvalidNeighbour) {
			nextTnextS = neighbours[curSnextT][NextSpace].index;


			if (nextTnextS != Globals::InvalidNeighbour) {
				// (t + 1, s) -> (t + 1, s + 1)
				Vector3 a = vertices[curSnextT];
				Vector3 b = vertices[nextTnextS];
				valid[3] = Geometry::doesPlaneIntersectEdge(time, a, b, intersections[numIntersections]);
				if (valid[3])
					numIntersections++;
			}
		}

		if (numIntersections >= 2) 
		{
			Vector3 tmpIntersections[2];
			int count = 0;
			for (int j = 0; j < 4; j++) {
				if (count == 2) break;
				if (valid[j]) {
					tmpIntersections[count] = intersections[count];
					count++;
				}
			}
			edges.push_back(Edge2D(tmpIntersections[0], tmpIntersections[1]));
			
		}
		if (numIntersections > 0 && numIntersections != 2) {
			errorFile << globalToFrameMap[i].first << " " <<  globalToFrameMap[i].second << " " << vertices[i] << " " <<numIntersections << endl;
			errorFile << curSnextT << " " << globalToFrameMap[curSnextT].first << " " <<  globalToFrameMap[curSnextT].second << " " << vertices[curSnextT] << endl;
			errorFile << nextScurT << " " << globalToFrameMap[nextScurT].first << " " <<  globalToFrameMap[nextScurT].second << " " << vertices[nextScurT] <<endl;
			errorFile << nextSnextT <<" " << globalToFrameMap[nextSnextT].first << " " <<  globalToFrameMap[nextSnextT].second <<  " " << vertices[nextSnextT] << endl;			
			errorFile << nextTnextS << " " << globalToFrameMap[nextTnextS].first << " " <<  globalToFrameMap[nextTnextS].second <<  " " << vertices[nextTnextS] <<endl;			

			for (int j = 0; j < 4; j++) {
				errorFile << intersections[j] << " " << valid[j] << endl;
			}
			errorFile << endl;
		
		}
	}

	errorFile.close();
	
	fstream f(outFilename+".ms", ios::out);
	try {
		f << edges.size() << endl;
		for (int i = 0; i < edges.size(); i++) {
				f << edges[i].start[0] << " " << edges[i].start[1] << " " << edges[i].end[0] << " " << edges[i].end[1] << endl;
			}		
		f.close();
	}
	catch (exception e) {
		throw;
	}
}

void Animation2D::identifyTopologyChanges() {	

	//// Set up for detecting topology changes
	//topoChangeVertices.clear();	
	//for (int i = 0; i < vertices.size(); i++) {
	//	topoChangeVertices.push_back(false);		
	//}

	for (int i = 0; i < vertices.size(); i++) {			
		if  ( (!isVertexOnStartFrame(i) && neighbours[i][PrevTime].index == Globals::InvalidNeighbour) ||
			(!isVertexOnEndFrame(i) && neighbours[i][NextTime].index == Globals::InvalidNeighbour)) {
			tbb::concurrent_hash_map<int, bool>::accessor a;
			topoChangeVertices.insert(a, pair<int, bool>(i, true));					
		}
	}

	for (int i = 0; i < vertices.size(); i++) {
		tbb::concurrent_hash_map<int, bool>::const_accessor tmpBool;
		if (topoChangeVertices.find(tmpBool,i)){
			continue;
		}
		
			

		vector<size_t> indices = findNClosestVertices(vertices[i], 2);
		int closestNbIndex = indices[1];
		Float d = (vertices[closestNbIndex] - vertices[i]).norm();
		bool found = false;
		// Test neighbours to see if closest vertex is 
		for (int j = 0; j < 2; j++) {
			if (neighbours[i][j].index != Globals::InvalidNeighbour) {
				if (neighbours[i][j].index == closestNbIndex) {
					found = true;
					break;
				}
			}
		}

		if (!found && d < 0.75 * timeScaling) {
			tbb::concurrent_hash_map<int, bool>::accessor a;
			topoChangeVertices.insert(a, pair<int, bool>(i, true));					
		}
	}

	vector<int> tmp;
	for (int i = 0; i < vertices.size(); i++) {
		for (int j = 0; j < 2; j++) {
			if (neighbours[i][j].index != Globals::InvalidNeighbour) {
				tbb::concurrent_hash_map<int, bool>::const_accessor tmpBool;
				if (topoChangeVertices.find(tmpBool,neighbours[i][j].index)){
					tmp.push_back(i);
				}
			}
		}
	}

	for (int i = 0; i < tmp.size(); i++) {
		tbb::concurrent_hash_map<int, bool>::accessor a;
		topoChangeVertices.insert(a, pair<int, bool>(i, true));		
	}

	//fstream f("topo.txt", ios::app);
	//f << "Topology change vertices:" << endl;
	//for (int i = 0; i < vertices.size(); i++) {
	//	if (topoChangeVertices[i])
	//		f << i << endl;
	//}
	//f.close();
}


int Animation2D::getNextSpaceIndex(GlobalIndex index) const {
	FrameIndex fIndex = globalToFrameMap.find(index)->second;
	int time = fIndex.first;
	int space = fIndex.second;
	return frameToGlobalMap.find(FrameIndex(time, dynamic_pointer_cast<const Mesh2D>(frames[time])->nextVertexIndex(space)))->second;
}


int Animation2D::getPrevSpaceIndex(GlobalIndex index) const {
	FrameIndex fIndex = globalToFrameMap.find(index)->second;
	int time = fIndex.first;
	int space = fIndex.second;
	return frameToGlobalMap.find(FrameIndex(time, dynamic_pointer_cast<const Mesh2D>(frames[time])->prevVertexIndex(space)))->second;
}

