#include "stdafx.h"
#include "Animation.h"
#include "Geometry.h"
#include <algorithm> 
#include "EdgeLabels.h"
#include <Bench/BenchTimer.h>
#include "Mesh3D.h"
using namespace Geometry;

/////////////////////////////////////////////////////////////////////////////
//  Specialization for Animation3D
/////////////////////////////////////////////////////////////////////////////



// Rotate about X axis
void Animation3D::rotateX(Float theta) {
	Float rad = theta * 3.14159 / 180.0;
	MatrixN rot;

	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			rot(i,j) = 0.0;
		}
	}
	rot(0,0) = 1;
	rot(1,1) = cos(rad);
	rot(1,2) = -sin(rad);
	rot(2,1) = sin(rad);
	rot(2,2) = cos(rad);
	rot(3,3) = 1;

	tbb::parallel_for(size_t(0), size_t(vertices.size()), size_t(1), [&](int i)	{
		vertices[i] = rot * vertices[i];
	}
	);

}


// Load meshes and associated velocity fields

bool Animation3D::loadMeshes(const string path,const int startFrame, const int endFrame) {
	bool failed = false;

	
	// Create meshes and velocities for each frame
	for (int i = startFrame; i <= endFrame; i++) {
		frames.emplace_back(new Mesh3D());
	}


	LOG(INFO) << "About to load meshes ... ";
	tbb::parallel_for(startFrame, endFrame+1, 1, [&](int i) 	
	//for (int i = startFrame; i < endFrame+1; i++)
	{
		char meshFilename[1024];
		sprintf_s(meshFilename, 1024, "%s_%04d.bobj.gz", path.c_str(), i);

		VLOG(1) << "About to dynamic cast ... ";

		try {
			if (!dynamic_pointer_cast<Mesh3D>(frames[i-startFrame])->loadBinaryObj(meshFilename)) {				
				LOG(ERROR) << "Failed while loading";
				failed = true;			
			}			
		}
		catch (exception e) {
			LOG(ERROR) << e.what();
			throw e;
		}
		char velFilename[1024];
		sprintf_s(velFilename, 1024, "%s_%04d.vvel", path.c_str(), i);

		if (!dynamic_pointer_cast<Mesh3D>(frames[i-startFrame])->loadVelocities(velFilename)) {				
			failed = true;			
		}	
	}
	);

	return !failed;
}

// Compute time scaling i.e. average vertex spacing using the first frame 

void Animation3D::computeTimeScaling()
{	
	timeScaling = 0.0;
	const Mesh3D* mesh = dynamic_pointer_cast<Mesh3D>(frames[0]).get();
	Float sum = 0.0;
	int count = 0;
	for (int i = 0; i < mesh->vertices.size(); i++) {
		vector<int> oneRing = mesh->getOneRing(i);

		FrameIndex fIndex = globalToFrameMap[i];
		if (oneRing.size() > 10)
			VLOG(1) << "Vertex: " << i << " : " << oneRing.size();

		for (int j = 0; j < oneRing.size(); j++) {
			sum += (mesh->vertices[oneRing[j]].pos - mesh->vertices[i].pos).norm();
			count++;
		}
	}
	timeScaling = sum / count;
	Globals::AvgEdgeLength = timeScaling;
	VLOG(1) << "Average edge length: " << timeScaling ;
}

// Computes the closest point and 
Correspondence<3> Animation3D::findClosestPointOnTetrahedron(const VectorN& p, const VectorN& a, const VectorN& b, const VectorN& c, const VectorN& d) const {
	Correspondence<3> corr;

	Vector4 n = Geometry::computeOrthogonalVector(a-d,b-d,c-d);

	

	Float dist = (p-d).dot(n);
	corr.n = n;	
	corr.p = p - dist*n;
	corr.d = abs(dist);

	// Find max component
	Float maxComponent = 0;
	int maxIndex = 0;
	for (int i = 0; i < 4; i++) {
		if (abs(n[i])> maxComponent) {
			maxComponent = abs(n[i]);
			maxIndex = i;
		}
	}

	// Project down 
	Vector3 newA,newB,newC,newD,newP;
	int count = 0;
	for (int i = 0; i < 4; i++) {
		if (i != maxIndex) {
			newA[count] = a[i];
			newB[count] = b[i];
			newC[count] = c[i];
			newD[count] = d[i];
			newP[count] = p[i];
			count++;
		}		
	}

	Matrix3 T;
	Vector3 rhs;
	for (int i = 0; i < 3; i++) {
		T(i,0) = newA[i] - newD[i];
		T(i,1) = newB[i] - newD[i];
		T(i,2) = newC[i] - newD[i];
		rhs[i] = newP[i] - newD[i];
 	}

	Vector3 bary = T.inverse() * rhs;

	corr.valid = true;
	for (int i = 0; i < 3; i++) {
		if (bary[i] < -1e-6 || bary[i] > 1+1e-6) {
			corr.valid = false;
			break;
		}
	}
	

	VLOG(1) << "Tet corr:" << corr.n;

	return corr;
};



Correspondence<3> Animation3D::findClosestPointOnTriangleFan(const VectorN& p, int curIndex) const{	
	FrameIndex fIndex = globalToFrameMap.at(curIndex);
	int time = fIndex.first;
	int spaceIndex = fIndex.second;

	GlobalIndex curTriIndex[3];
	GlobalIndex nextTriIndex[3];
	Correspondence<3> corr;
	corr.valid = false;
	corr.d = 1e6;
	
	curTriIndex[0] = curIndex;
	nextTriIndex[0] = neighbours[curIndex][NextTime].index;

	if (nextTriIndex[0] != Globals::InvalidNeighbour) {
		vector<int> oneRing = dynamic_pointer_cast<const Mesh3D>(frames[time])->getOneRing(spaceIndex);

		for (int i = 0; i < oneRing.size()-1; i++) {

			FrameIndex fIndex (time,oneRing[i]);
			curTriIndex[1] = frameToGlobalMap.at(fIndex);

			 fIndex = FrameIndex(time,oneRing[i+1]);
			curTriIndex[2] = frameToGlobalMap.at(fIndex);

			bool validPrism = true;
			for (int j = 1; j < 3; j++) {
				nextTriIndex[j] = neighbours[curTriIndex[j]][NextTime].index;
				if (nextTriIndex[j] == Globals::InvalidNeighbour) {
					validPrism = false;
					break;
				}
			}

			if (nextTriIndex[0]==nextTriIndex[1] || nextTriIndex[1] == nextTriIndex[2] || nextTriIndex[2] == nextTriIndex[0]) {
				validPrism = false;
				break;
			}

			if (validPrism) {
				// Three tets per prism

				//doesPlaneIntersectTetrahedron(time, curTriangle[0], curTriangle[1], curTriangle[2], futureTriangle[0], tmpTriangles);
				//doesPlaneIntersectTetrahedron(time, futureTriangle[0], curTriangle[1], curTriangle[2], futureTriangle[1], tmpTriangles);
				//doesPlaneIntersectTetrahedron(time, futureTriangle[0], futureTriangle[1], futureTriangle[2], curTriangle[2], tmpTriangles);
				Correspondence<3> tetCorrs[3];
				tetCorrs[0] = findClosestPointOnTetrahedron(p, vertices[ curTriIndex[0] ], vertices[ curTriIndex[1] ], vertices[ curTriIndex[2] ], vertices[ nextTriIndex[0] ]);
				tetCorrs[1] = findClosestPointOnTetrahedron(p, vertices[ curTriIndex[1] ], vertices[ nextTriIndex[1] ], vertices[ nextTriIndex[2] ], vertices[ nextTriIndex[0] ] );				
				tetCorrs[2] = findClosestPointOnTetrahedron(p, vertices[ curTriIndex[1] ], vertices[ nextTriIndex[2] ], vertices[ curTriIndex[2] ], vertices[ nextTriIndex[0] ]);				

				for (int k = 0; k < 3; k++) {
					if (tetCorrs[k].valid && tetCorrs[k].d < corr.d) {
						corr = tetCorrs[k];
					}
				}
			}
		}
	}
	return corr;
}


/*
Finds the closest point and normal to a given point (not necessarily vertex)
*/

Correspondence<3> Animation3D::findClosestPoint(const VectorN& queryPoint) const{
	const GlobalIndex closestVertexIndex = findClosestVertex(queryPoint);
	Correspondence<3> prevFanCorr, curFanCorr, corr;

	//corr.p = vertices[closestVertexIndex];
	//corr.n = normals[closestVertexIndex];
	//return corr;
	
	corr.p = vertices[closestVertexIndex];
	corr.d = (queryPoint - vertices[closestVertexIndex]).norm();
	corr.n = normals[closestVertexIndex];
	corr.valid = true;

	//return corr;

	if (neighbours[closestVertexIndex][PrevTime].index != Globals::InvalidNeighbour) {
		prevFanCorr = findClosestPointOnTriangleFan(queryPoint,neighbours[closestVertexIndex][PrevTime].index);
		if (prevFanCorr.valid) {
			if (curFanCorr.d < corr.d) 
				corr = prevFanCorr;
		}
	}
	if (neighbours[closestVertexIndex][NextTime].index != Globals::InvalidNeighbour) {
		curFanCorr = findClosestPointOnTriangleFan(queryPoint,closestVertexIndex);

		if (curFanCorr.valid) {
			if (curFanCorr.d < corr.d) {
				corr = curFanCorr;
			}
		}
	}

	VLOG(1) << "Closest vertex normal:" << closestVertexIndex << " " << normals[closestVertexIndex].format(CommaFmt);
	VLOG(1) << "Corr normal:" << corr.n;
	
	return corr;
}


void Animation3D::identifyTopologyChanges() {
	LOG(INFO) << "Identifying topology changes and thin sheets";

	tbb::parallel_for(size_t(0), vertices.size(), size_t(1), [&](int i) 
	{			
		if  (!isVertexOnEndFrame(i) && neighbours[i][NextTime].index == Globals::InvalidNeighbour) {
			tbb::concurrent_hash_map<int, bool>::accessor a;
			topoChangeVertices.insert(a, pair<int, bool>(i, true));			
			
		}
		else {
			const int maxNbs = 50;
			vector<size_t> indices = findNClosestVertices(vertices[i], maxNbs);
			FrameIndex fIndex = globalToFrameMap.at(i);
			//VLOG(1) << fIndex.first << "," << fIndex.second << " normal:" << normals[i].format(CommaFmt) << endl;
			int numOppNbs = 0; // number of nbs with opposite normals
			for (int j = 0; j < indices.size(); j++) {
				int nbIndex = indices[j];
				if (globalToFrameMap.at(nbIndex).first == fIndex.first) {
					if (normals[nbIndex].dot(normals[i]) < 0.0) {
						numOppNbs++;
						break;
					}
				}
			}

			if (numOppNbs >= 1) {
				tbb::concurrent_hash_map<int, bool>::accessor a;
				topoChangeVertices.insert(a, pair<int, bool>(i, true));			

				VLOG(1) << "Thin sheet at:" << fIndex.first << "," << fIndex.second;
			}
		}
	}
	);

	tbb::concurrent_vector<int> nbVertices;
	tbb::parallel_for(size_t(0), vertices.size(), size_t(1), [&](int i) 
	{			
		for (int j = 0; j < neighbours[i].size(); j++) {
			int nbIndex = neighbours[i][j].index;
			tbb::concurrent_hash_map<int, bool>::const_accessor tmpBool;
			if (topoChangeVertices.find(tmpBool, nbIndex))
				nbVertices.push_back(i);
		}
	});

	for (int i = 0; i < nbVertices.size(); i++) {
		tbb::concurrent_hash_map<int, bool>::accessor a;
		topoChangeVertices.insert(a, pair<int, bool>(nbVertices[i], true));	
	}

	
}


void Animation3D::findNeighbours() {
	LOG(INFO)<<"Finding neighbours on the spacetime mesh";
	const Float maxSqrDistance = 4 * timeScaling * timeScaling;
	tbb::mutex   countMutex;

	neighbours.clear();
	neighbours.reserve(vertices.size());

	for (int i = 0; i < vertices.size(); i++) 	
	{
		neighbours.emplace_back(vector<Neighbour>());
	}

	int numResampled = 0;

	#ifdef ENABLE_PARALLEL
	tbb::parallel_for(size_t(0), vertices.size(), size_t(1), [&](int i)	
	#else
	for (int i = 0; i < vertices.size(); i++) 
	#endif
	{
		//VLOG(1) << "Vertex:" << i << " / " << vertices.size();
		FrameIndex fIndex = globalToFrameMap[i];
		int time = fIndex.first;
		int spaceIndex = fIndex.second;


		const Mesh3D* mesh = dynamic_pointer_cast<Mesh3D>(frames[time]).get();
		int curId = mesh->vertices[spaceIndex].id;
		vector<Neighbour> nbs;
		nbs.reserve(30);

		// Create entries for neighbours in time		
		Neighbour n;
		n.index = Globals::InvalidNeighbour;
		Vector3 curVel(0,0,0);
		nbs.push_back(n);
		nbs.push_back(n);

		// Predict forward and backward positions
		curVel = frames[time]->velocities[spaceIndex];
		Vector4 nextPos = Vector4( vertices[i][0] + dt * curVel[0], vertices[i][1] + dt * curVel[1], vertices[i][2] + dt * curVel[2], vertices[i][3] + getDepth(1));

		Vector4 prevPos = Vector4( vertices[i][0] - dt * curVel[0], vertices[i][1] - dt * curVel[1], vertices[i][2] - dt * curVel[2], vertices[i][3] - getDepth(1));


		// Skip first frame for prevTime
		if (time > 0) {
			auto iter = frames[time - 1]->idToIndexMap.find(curId);
			if (iter != frames[time - 1]->idToIndexMap.end()) {
				nbs[PrevTime].index = frameToGlobalMap[FrameIndex(time-1, iter->second)];
				Float dSqr = (vertices[nbs[PrevTime].index] - vertices[i]).squaredNorm();
				nbs[PrevTime].d = sqrt(dSqr);
			}
			else  {
				//Eigen::BenchTimer timer;
				//timer.start();

				int closestVertexIndex = findClosestVertex(prevPos);

				//timer.stop();
				//VLOG(2)<<"Closest index: " << closestVertexIndex;
				//VLOG(2)<<"a:" << vertices[i].format(CommaFmt);
				//VLOG(2)<<"b:" << vertices[closestVertexIndex].format(CommaFmt);
				//VLOG(2)<<"Time:"<<timer.total();


				Float dSqr = (vertices[closestVertexIndex] - prevPos).squaredNorm();

				if (dSqr < maxSqrDistance && globalToFrameMap[closestVertexIndex].first == time - 1) {


					nbs[PrevTime].index = closestVertexIndex;
					nbs[PrevTime].d =  (vertices[ closestVertexIndex ] -  vertices[i]).norm();
				}
				else {
					nbs[PrevTime].index = Globals::InvalidNeighbour;
					nbs[PrevTime].d = 1e6;
				}
			}
		}

		// Skip last frame for nextTime
		if (time <  frames.size() - 1) {
			auto iter = frames[time + 1]->idToIndexMap.find(curId);
			if (iter != frames[time + 1]->idToIndexMap.end()) {
				nbs[NextTime].index = frameToGlobalMap[ FrameIndex(time + 1, iter->second) ];
				Float dSqr = (vertices[nbs[NextTime].index] - vertices[i]).squaredNorm();
				nbs[NextTime].d = sqrt(dSqr);

			}
			else {

				int closestVertexIndex = findClosestVertex(nextPos);				
				Float dSqr = (vertices[closestVertexIndex] - nextPos).squaredNorm();


				if (dSqr < maxSqrDistance && globalToFrameMap[closestVertexIndex].first == time + 1) {
					//VLOG(1) <<"Time : " << time << " Vertex:" << i << " was resampled " << vertices[i] << " " << vertices[closestVertexIndex].format(CommaFmt) << " " << curVel;
					nbs[NextTime].index = closestVertexIndex;
					nbs[NextTime].d = (vertices[ closestVertexIndex ] -  vertices[i]).norm();
				}
				else {					


					vector<size_t> closestVertices = findNClosestVertices(nextPos, 10);					
					int count = 0;

					for (int k = 0; k < closestVertices.size(); k++) {

						FrameIndex fIndex = globalToFrameMap.at(closestVertices[k]);						
						if (fIndex.first == time+1) {							

							nbs[NextTime].index = closestVertices[k];
							nbs[NextTime].d = (vertices[ closestVertices[k] ] -  vertices[i]).norm();
							count++;
							break;
						}
					}

					if (count != 1)	{
						//LOG(ERROR) << "Could not find nearby vertex for: (" << time << "," << spaceIndex << ") " << vertices[i] << " | " << nextPos << "|" << curVel << "|" << dt << endl;
						nbs[NextTime].index = Globals::InvalidNeighbour;
						nbs[NextTime].d = 1e6;

					}
				}
				//numResampled++;				
			}
		}

		// Create entries for spatial neighbours		
		vector<int> oneRing = mesh->getOneRing(spaceIndex);
		for (int j = 0; j < oneRing.size(); j++) {
			Neighbour n;
			GlobalIndex globalIndex = frameToGlobalMap[FrameIndex(time, oneRing[j])];
			Float d = (vertices[i] - vertices[globalIndex]).norm();
			n.index = globalIndex;
			n.d = d;			
			nbs.emplace_back(n);
		}

		neighbours[i] = nbs;

		
	}
	#ifdef ENABLE_PARALLEL
	);
	#endif

	
}

// Compute spacetime normals in parallel

void Animation3D::computeNormals() {
	normals.clear();
	for (int i = 0; i < vertices.size(); i++)
	{
		normals.push_back(Vector4(0,0,0,0));
	}

#ifdef ENABLE_PARALLEL
	tbb::parallel_for(size_t(0), vertices.size(), size_t(1), [&](int i)	
#else
	for (int i = 0; i < vertices.size(); i++)
#endif

	{		
		
		normals[i] = computeNormalAtVertex(i);
	}

#ifdef ENABLE_PARALLEL
	);
#endif

	//writeNormalsToDisk(path);
}


// Computes normal at vertex using global index

Vector4 Animation3D:: computeNormalAtVertex(int index) {
	FrameIndex frameIndex = globalToFrameMap[index];
	int timeIndex = frameIndex.first;
	int spaceIndex = frameIndex.second;


	vector<int> oneRing = dynamic_pointer_cast<Mesh3D>(frames[timeIndex])->getOneRing(spaceIndex);

	Vector4 finalNormal(0,0,0,0);

	if (oneRing.size() == 0) 
		return finalNormal;

	for (int i = 0; i < oneRing.size()-1; i++) {

		int gIndex0 = frameToGlobalMap[FrameIndex(timeIndex, oneRing[i])];

		if (frameToGlobalMap.find(FrameIndex(timeIndex, oneRing[(i+1)%oneRing.size()])) == frameToGlobalMap.end())
			continue;

		int gIndex1 = frameToGlobalMap[FrameIndex(timeIndex, oneRing[(i+1)%oneRing.size()])];

		Vector4 st0 = vertices[gIndex0] - vertices[index];

		Vector4 st1 = vertices[gIndex1] - vertices[index];

		Vector4 st2;

		if (neighbours[index][NextTime].index != Globals::InvalidNeighbour) {
			st2 = vertices[neighbours[index][NextTime].index] - vertices[index];
			finalNormal += Geometry::computeOrthogonalVector(st0,st2,st1);
		}
		else if (neighbours[index][PrevTime].index != Globals::InvalidNeighbour) {
			st2 = vertices[index] - vertices[neighbours[index][PrevTime].index];
			finalNormal += Geometry::computeOrthogonalVector(st0,st2,st1);
		}	

	}


	finalNormal /= oneRing.size();
	Float norm = finalNormal.norm();

	if (norm < 1e-15) {
		//LOG(ERROR) << "Vertex : " << timeIndex << "," << spaceIndex << " has normal :" << finalNormal << " " << norm;
	}
	finalNormal.normalize();


	return finalNormal;

}



//////////////////////////////////////
 
bool Animation3D::doesFutureTriangleExist(int time, int faceIndex, vector<Geometry::Vertex4>& curTriangle, vector<Geometry::Vertex4>& futureTriangle) {
	futureTriangle.clear();
	curTriangle.clear();

	int nextSpaceIndices[3];

	for (int k = 0; k < 3; k++) {
		int sIndex = dynamic_pointer_cast<Mesh3D>(frames[time])->faces[faceIndex].v[k];
		FrameIndex fIndex(time, sIndex);
		GlobalIndex gIndex = frameToGlobalMap[fIndex];

		curTriangle.push_back(Geometry::Vertex4(vertices[gIndex], gIndex));

		if (neighbours[gIndex][NextTime].index == Globals::InvalidNeighbour) {			
			//VLOG(1) << time << "," << sIndex << " Invalid future neighbour";
				return false;			
		}				
		else {
			
			nextSpaceIndices[k] = globalToFrameMap[neighbours[gIndex][NextTime].index].second;
			futureTriangle.push_back(Geometry::Vertex4(vertices[ neighbours[gIndex][NextTime].index ], neighbours[gIndex][NextTime].index));		
		}
		

		//VLOG(1) << "Future index:" << faceIndex << " " << time << " " << sIndex << " " << gIndex << " " << neighbours[gIndex][NextTime].index;
	}

	return true;
	auto generateFaceHash = [] (int a, int b, int c) -> size_t {
		size_t seed = 0;
		boost::hash_combine(seed, a);
		boost::hash_combine(seed, b);
		boost::hash_combine(seed, c);
		return seed;
	};


	Mesh3D* nextFrame = dynamic_pointer_cast<Mesh3D>(frames[time+1]).get();	
	
	return !(nextFrame->faceHashMap.find(generateFaceHash(nextSpaceIndices[0],nextSpaceIndices[1],nextSpaceIndices[2])) == nextFrame->faceHashMap.end());
			
}


void Animation3D::setupExtraction() {
	
	frameEdgeLabels.clear();
	for (int i = 0; i < frames.size(); i++) {
		frameEdgeLabels.emplace_back(EdgeLabels());
	}

	tbb::parallel_for(size_t(0), size_t(frames.size()), size_t(1), [&](int i)  
	{		
		labelEdges(*dynamic_pointer_cast<Mesh3D>(frames[i]).get(), frameEdgeLabels[i]);
		
	});
}

void fillHoles(Mesh3D& mesh) {
	unordered_set < pair<int,int>, boost::hash< std::pair<int, int> >  > boundaryEdges;
	unordered_set<int> boundaryVerts;

	LOG(INFO) << "Building corner table";
	mesh.buildCornerTable();
	//
	//LOG(INFO) << "Finding holes";

	int numCorners = 3 * mesh.faces.size();
	for (int i = 0; i < numCorners; i++) {
		int c = i;

		if (mesh.b(c)) {
			vector<int> bv = mesh.boundaryVertices(c);


			bool valid = true;

			for (int j = 0; j < bv.size(); j++) {
				int cur = bv[j];
				int next = bv[(j+1)%bv.size()];
				auto edge = pair<int,int> ( min(cur,next), max(cur,next));
				if (boundaryEdges.find(edge) != boundaryEdges.end()) {
					valid = false;
					break;					
				}				
			}

			if (valid) {
				//LOG(INFO) << "Boundary at corner:" << c << endl;
				//LOG(INFO) << "Vertex:" << mesh.v(c) << endl;
				Vertex3D vert;
				vert.pos = Vector3(0,0,0);

				for (int j = 0; j < bv.size() ; j++) {
					int cur = bv[j];
					int next = bv[(j+1)%bv.size()];
					auto edge = pair<int,int> ( min(cur,next), max(cur,next));
					boundaryEdges.insert(edge);
					//LOG(INFO) << bv[j] << " ";

					boundaryVerts.insert(bv[j]);
					vert.pos += mesh.vertices[bv[j]].pos;
				}

				vert.pos /= bv.size();								
				mesh.vertices.push_back(vert);

				for (int j = 0; j < bv.size() ; j++) {
					int cur = bv[j];
					int next = bv[(j+1)%bv.size()];

					Face3D face;					
					face.v[0] = next;
					face.v[1] = cur;
					face.v[2] = mesh.vertices.size()-1;			


					boundaryVerts.insert(face.v[2]);
					mesh.faces.push_back(face);
				}

				//LOG(INFO) << endl;
			}			
		}
	}

	LOG(INFO) << "Building 2nd corner table";
	mesh.buildCornerTable();


	LOG(INFO) << "Smoothing ... ";
	// Performs Laplace smoothing on data. 

	
	int numSmoothingPasses = 6;

	Float t_lambda=0.6307;
	Float t_mu=-0.6732;

	// Try taubin smoothing

	for (int i = 0; i < numSmoothingPasses; i++) {

		auto tmpData = mesh.vertices;

		tbb::parallel_for(size_t(0), size_t(mesh.vertices.size()), size_t(1), [&](int j)  
		{			
			Vertex3D sum;
			Float sumWeights = 0;
			vector<int> adjacentVertices = mesh.getOneRing(j);

			{
				bool nearBoundaryVert = false;
				for (int k = 0; k < adjacentVertices.size(); k++) {

					//if (boundaryVerts.find(adjacentVertices[k]) != boundaryVerts.end())
					//{
					//	nearBoundaryVert = true;
					//	break;
					//}
					if (k==0)
						sum.pos = mesh.vertices[adjacentVertices[k]].pos - mesh.vertices[j].pos;
					else
						sum.pos += mesh.vertices[adjacentVertices[k]].pos - mesh.vertices[j].pos;

				}			

				//&& boundaryVerts.find(j) == boundaryVerts.end() && !nearBoundaryVert
				if (adjacentVertices.size() > 0 ) {
					if (i%2 == 0)
						tmpData[j].pos =  mesh.vertices[j].pos + t_lambda * sum.pos / adjacentVertices.size();
					else 
						tmpData[j].pos =  mesh.vertices[j].pos + t_mu * sum.pos / adjacentVertices.size();
				}

			}			

		});

		mesh.vertices = tmpData;
	}




	LOG(INFO) << "Done!";
}


void Animation3D::extractSlice(string outFilename, Float time, bool perturb) {

	//LOG(INFO) << "Extracting frames ... ";

	if (perturb) {
		perturbVertices();
	}

	fstream f;
	//f.open(outFilename+".invalid",ios::out);

	// Construct tetrahedra for the space-time volume
	vector<Geometry::Triangle> triangles;

	tbb::mutex triangleMutex;

	//LOG(INFO)<<"Extracting frame at time: "<<time;

	int faceCount = 0;
	for (int frame = 0; frame < frames.size()-1; frame++) {		
		if (min(minTime[frame],minTime[frame+1]) < time && max(maxTime[frame],maxTime[frame+1]) > time) {			
			// Traverse mesh

			//LOG(INFO) << "Using triangles from frame:" << frame;

#ifdef ENABLE_PARALLEL
			tbb::parallel_for(size_t(0), size_t(dynamic_pointer_cast<Mesh3D>(frames[frame])->faces.size()), size_t(1), [&](int faceIndex)  			
#else
			for (int faceIndex = 0; faceIndex < frames[frame]->faces.size(); faceIndex++) 
#endif
			{

				// Should we try to slice this into tetrahedra?
				vector<Vertex4> curTriangle,futureTriangle;				

				if (doesFutureTriangleExist(frame, faceIndex, curTriangle, futureTriangle)) {					
					if ( 
						(curTriangle[0].pos[3] < time && futureTriangle[0].pos[3] > time) ||
						(curTriangle[1].pos[3] < time && futureTriangle[1].pos[3] > time) ||
						(curTriangle[2].pos[3] < time && futureTriangle[2].pos[3] > time)) {

							// Create tetrahedra
							// Run the doesPlaneIntersectTetrahedron test and add triangles

							vector<Geometry::Triangle> tmpTriangles;
							tmpTriangles.clear();

							// Look for a sequence 001 or 110
							const EdgeLabels& edgeLabels = frameEdgeLabels[frame];
							const Mesh3D* mesh = dynamic_pointer_cast<Mesh3D>(frames[frame]).get();
							bool labels[3];

							//f << "Face:" << faceIndex << endl;
							for (int k = 0; k < 3; k++) {
								
								findEdge(edgeLabels, mesh->faces[faceIndex].v[k], mesh->faces[faceIndex].v[(k+1)%3], labels[k]);
								//f << mesh.faces[faceIndex].v[k] << " : " << labels[k] << endl;								
							}

							int indexA = 0;
							for (int k = 0; k < 3; k++) {
								if (labels[k] == labels[(k+1)%3]) {
									indexA = k;
									break;
								}
							}
							int indexB = (indexA+1)%3;
							int indexC = (indexA+2)%3;

							if (labels[indexA] == 0) {
								doesPlaneIntersectTetrahedron(time, futureTriangle[indexA], futureTriangle[indexB], futureTriangle[indexC], curTriangle[indexA], tmpTriangles, faceIndex, 0);
								doesPlaneIntersectTetrahedron(time, curTriangle[indexB], futureTriangle[indexB], futureTriangle[indexC], curTriangle[indexA], tmpTriangles, faceIndex, 1);
								doesPlaneIntersectTetrahedron(time, curTriangle[indexB], futureTriangle[indexC], curTriangle[indexC], curTriangle[indexA], tmpTriangles, faceIndex, 2);
							}
							else {
								doesPlaneIntersectTetrahedron(time, curTriangle[indexA], curTriangle[indexB], curTriangle[indexC], futureTriangle[indexA], tmpTriangles, faceIndex, 3);
								doesPlaneIntersectTetrahedron(time, curTriangle[indexB], curTriangle[indexC], futureTriangle[indexA], futureTriangle[indexB], tmpTriangles, faceIndex, 4);
								doesPlaneIntersectTetrahedron(time, curTriangle[indexC], futureTriangle[indexB], futureTriangle[indexA], futureTriangle[indexC], tmpTriangles, faceIndex, 5);
							}

							// Orient the triangle s
							for (int k = 0; k < tmpTriangles.size(); k++) {
								Vector4 origA = curTriangle[1].pos - curTriangle[0].pos;
								Vector4 origB = curTriangle[2].pos - curTriangle[0].pos;								

								Vector3 origA3 (origA[0],origA[1],origA[2]);
								Vector3 origB3 (origB[0],origB[1],origB[2]);

								Vector3 origN = origA3.cross(origB3);

								Vector4 tmpA = tmpTriangles[k].b.pos - tmpTriangles[k].a.pos;
								Vector4 tmpB = tmpTriangles[k].c.pos - tmpTriangles[k].a.pos;

								Vector3 tmpA3 (tmpA[0],tmpA[1],tmpA[2]);
								Vector3 tmpB3 (tmpB[0],tmpB[1],tmpB[2]);
								Vector3 tmpN = tmpA3.cross(tmpB3);

								if (tmpN.dot(origN) < 0.0) {
									swap(tmpTriangles[k].a, tmpTriangles[k].b);
								}

								 
							}


							{
								tbb::mutex::scoped_lock lock(triangleMutex); 
								for (int k = 0; k < tmpTriangles.size(); k++) {									
									triangles.emplace_back(tmpTriangles[k]);									
								}
							}
					}
				}
				else {
					// Do nothing with unused triangles
				
				}
			}		
#ifdef ENABLE_PARALLEL
			);
#endif
		}		
	}


	
	
	// Combine triangles into a mesh
	Mesh3D outMesh;	
	// Which space-time edge did this vertex come from?
	typedef pair<GlobalIndex,GlobalIndex> IntersectionEdge;
	typedef pair<int,int> MeshEdge;

	map < IntersectionEdge , int > vertexMap;
	map < MeshEdge, bool > edgeMap; // True implies good (non-boundary edge)
	
	vector< Face3D > faces;
	int vertexCount = 0;

	unordered_set<int> boundaryVertices;

	vector< vector<int> > adjacentVertices;


	//f.open(outFilename+".debug", ios::out);

	for (int i = 0; i < triangles.size(); i++) {
		Vertex3D a,b,c;		
		for (int j = 0; j < 3; j++) {
			a.pos[j] = triangles[i].a.pos[j];
			b.pos[j] = triangles[i].b.pos[j];
			c.pos[j] = triangles[i].c.pos[j];
		}

		
		IntersectionEdge edgeA = IntersectionEdge(triangles[i].a.index0, triangles[i].a.index1);
		IntersectionEdge edgeB = IntersectionEdge(triangles[i].b.index0, triangles[i].b.index1);
		IntersectionEdge edgeC = IntersectionEdge(triangles[i].c.index0, triangles[i].c.index1);
		



		if (vertexMap.find(edgeA) == vertexMap.end()) {
			vertexMap[edgeA] = vertexCount;
			outMesh.vertices.push_back(a);
			adjacentVertices.push_back(vector<int>());
			vertexCount++;

		}
		if (vertexMap.find(edgeB) == vertexMap.end()) {
			vertexMap[edgeB] = vertexCount;
			outMesh.vertices.push_back(b);
			adjacentVertices.push_back(vector<int>());
			vertexCount++;

		}
		if (vertexMap.find(edgeC) == vertexMap.end()) {
			vertexMap[edgeC] = vertexCount;
			outMesh.vertices.push_back(c);
			adjacentVertices.push_back(vector<int>());
			vertexCount++;
		}


		Face3D face;		
		face.v[0] = vertexMap.at(edgeA);
		face.v[1] = vertexMap.at(edgeB);
		face.v[2] = vertexMap.at(edgeC);

		// Identify boundaries
		for (int k = 0; k < 3; k++) {
			MeshEdge meshEdge = MeshEdge(min(face.v[k],face.v[(k+1)%3]),max(face.v[k],face.v[(k+1)%3]));
			if (edgeMap.find(meshEdge) != edgeMap.end()) {
				edgeMap[meshEdge] = true;
			}
			else {
				edgeMap[meshEdge] = false;
			}
		}


		// Build up adjacency for vertices
		if (find(adjacentVertices[face.v[0]].begin(),adjacentVertices[face.v[0]].end(), face.v[1]) == adjacentVertices[face.v[0]].end()) {
			adjacentVertices[face.v[0]].push_back(face.v[1]);
			adjacentVertices[face.v[1]].push_back(face.v[0]);
		}
		if (find(adjacentVertices[face.v[1]].begin(),adjacentVertices[face.v[1]].end(), face.v[2]) == adjacentVertices[face.v[1]].end()) {
			adjacentVertices[face.v[1]].push_back(face.v[2]);
			adjacentVertices[face.v[2]].push_back(face.v[1]);
		}
		if (find(adjacentVertices[face.v[2]].begin(),adjacentVertices[face.v[2]].end(), face.v[0]) == adjacentVertices[face.v[2]].end()) {
			adjacentVertices[face.v[2]].push_back(face.v[0]);
			adjacentVertices[face.v[0]].push_back(face.v[2]);
		}

		// Skip degenerate triangles
		if ( face.v[0] == face.v[1] || face.v[1] == face.v[2] || face.v[0] == face.v[2] || triangles[i].origTriIndex == -1 ) {
 		} else {			
			int frame =  globalToFrameMap[triangles[i].a.index0].first;
			const EdgeLabels& edgeLabels = frameEdgeLabels[frame];
			const Mesh3D* mesh = dynamic_pointer_cast<Mesh3D>(frames[frame]).get();
			bool labels[3];
			for (int k = 0; k < 3; k++) {
				if (triangles[i].origTriIndex >= mesh->faces.size()) continue;
				findEdge(edgeLabels, mesh->faces[triangles[i].origTriIndex].v[k], mesh->faces[triangles[i].origTriIndex].v[(k+1)%3], labels[k]);
			}

			//f << outMesh.faces.size() << " " << triangles[i].origTriIndex << " " << triangles[i].type << endl;
			//for (int k = 0; k < 3; k++)
			//	f << mesh.faces[triangles[i].origTriIndex].v[k] << " " << labels[k] << endl;

			//f << "a:" << triangles[i].a.index0 << " " << triangles[i].a.index1 << endl;
			//f << "b:" << triangles[i].b.index0 << " " << triangles[i].b.index1 << endl;
			//f << "c:" << triangles[i].c.index0 << " " << triangles[i].c.index1 << endl;
			//f << "Cur:" << " (" << globalToFrameMap[triangles[i].a.index0].first << "," << globalToFrameMap[triangles[i].a.index0].second << ") "
			//	<< " (" << globalToFrameMap[triangles[i].b.index0].first << "," << globalToFrameMap[triangles[i].b.index0].second << ") "
			//	<< " (" << globalToFrameMap[triangles[i].c.index0].first << "," << globalToFrameMap[triangles[i].c.index0].second << ") " << endl;

			//f << "Fut:" << " (" << globalToFrameMap[triangles[i].a.index1].first << "," << globalToFrameMap[triangles[i].a.index1].second << ") "
			//	<< " (" << globalToFrameMap[triangles[i].b.index1].first << "," << globalToFrameMap[triangles[i].b.index1].second << ") "
			//	<< " (" << globalToFrameMap[triangles[i].c.index1].first << "," << globalToFrameMap[triangles[i].c.index1].second << ") " << endl;


			outMesh.faces.push_back(face);
		}
	}


	fillHoles(outMesh);	
	outMesh.writeBinaryObj(outFilename + ".bobj.gz");	
	//outMesh.writeObj(outFilename + ".obj");
}


void Animation3D::debugResamplingEvents() {
	vector<Geometry::Vertex4> curTriangle,futureTriangle;				

	for (int i = 0; i < frames.size(); i++) {
		char filename[1024];
		sprintf_s(filename, 1024, "frame_%04d.invalid", i);
		fstream f(filename,ios::out);
		const Mesh3D* frame = dynamic_pointer_cast<Mesh3D>(frames[i]).get();
		f << "maya.cmds.select([";
		for (int faceIndex = 0; faceIndex < frame->faces.size(); faceIndex++) {
			if (!doesFutureTriangleExist(i, faceIndex, curTriangle, futureTriangle)) {	
				f << "\"" << "AnimationA_Mesh.f[" << faceIndex << "]" << "\"" << ",";
			}
		}
		f << "])" << endl;
		f.close();
	}

	return;

}

