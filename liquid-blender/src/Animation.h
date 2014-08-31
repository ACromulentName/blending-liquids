#pragma once

#include "stdafx.h"
#include "KDTree.hpp"
#include "Mesh.h"
#include "Mesh2D.h"
#include "Mesh3D.h"
#include "Geometry.h"
#include <fstream>
#include <memory>
#include <set>
#include "Utility.hpp"
#include <unordered_map>


using namespace std;

typedef pair<int,int> FrameIndex;
typedef int GlobalIndex;
typedef vector< GlobalIndex > LineFeature;
typedef vector< LineFeature > LineFeatures;

enum Neighbours {PrevTime, NextTime, PrevSpace, NextSpace};

template<int DIM>
struct Correspondence {
	typedef Eigen::Matrix<Float, DIM+1, 1> VectorN;
	

	Correspondence() {
		valid = false;
		d = 0.0;
	}

	VectorN p; // Corresponding point
	VectorN n; // Corresponding normal	
	Float d;

	bool valid;
};

struct Bary {
	GlobalIndex gIndex[3];
	Float bary[3];
};

typedef pair<int,int> Edge;


struct FrameIndexHasher {
public:
	size_t operator()(FrameIndex const& v) const throw() {
		size_t seed = 0;
		boost::hash_combine(seed, v.first);
		boost::hash_combine(seed, v.second);
		return seed;
	}
};

typedef unordered_map <Edge, bool, FrameIndexHasher> EdgeLabels;

template<int DIM>
class Animation
{
public:
	typedef typename Eigen::Matrix<Float, DIM+1, 1> VectorN;
	typedef KDTreeVectorOfVectorsAdaptor< vector < VectorN >, Float, DIM+1>  KDTreeND;
	typedef typename Eigen::Matrix<Float, DIM+1, DIM+1> MatrixN;


	Animation() {
		rotX = 0;
		LOG(INFO) << "Creating animation";
	}

	virtual ~Animation() {
		cleanup();
	}

	virtual bool load(const string path, const int startFrame, const int endFrame, bool isBlendMode=false) { 
		cleanup();

		// Load simulation data
		char simFilename[1024];
		sprintf_s(simFilename, 1024, "%s.sim", path.c_str());

		// Read in timestep
		ifstream file;
		try {
			file.exceptions ( ifstream::eofbit | ifstream::failbit | ifstream::badbit );
			file.open(simFilename, ios::in);

			file >> dt;
			file.close();
		}
		catch (exception e) {
			LOG(ERROR) << "Exception occurred while loading " << simFilename << " : " << e.what();
			return false;
		}

		// Load in parallel 
		// C++11 lambdas are neat!
		bool isLoadSuccessful = loadMeshes(path, startFrame, endFrame);

		if (!isLoadSuccessful)
			return false;

		numFrames = frames.size();

		this->path = path;
		computeTimeScaling();
		createGlobalIndexing();

		createKDTree();
		findNeighbours();

		computeNormals();	
		if (!isBlendMode)
			identifyTopologyChanges();
		
		LOG(INFO) << "Animation loaded. # Vertices:" << vertices.size() << endl;

		return true;
	}


	// This is the minimal set of methods that need to be implemented by any derived Animation class.
	// This returns the closest point (vector) and normal as a pair

	virtual bool loadMeshes(const string path,const int startFrame, const int endFrame) = 0; 
	virtual Correspondence<DIM> findClosestPoint(const VectorN&) const = 0;
	virtual VectorN computeNormalAtVertex(int i) = 0;	

	virtual void computeTimeScaling() = 0;		
	virtual void findNeighbours() = 0;
	virtual void computeNormals() = 0;
	virtual void identifyTopologyChanges() = 0;

	virtual void setupExtraction() = 0;
	virtual void extractSlice(string outFilename, Float time, bool perturbVertices=false) = 0;

	virtual void cleanupFrames() {
		for (int i = 0; i < frames.size(); i++) {
			frames[i]->cleanup();
		}
	}

	virtual void cleanup() {
		LOG(INFO) << "Cleaning up animation";

		for (int i = 0; i < frames.size(); i++) {
			if (frames[i] != nullptr) {
				frames[i]->cleanup();				
			}
		}

		frames.clear();
		vertices.clear();
		frameToGlobalMap.clear();
		globalToFrameMap.clear();
		normals.clear();
	}

	virtual void createKDTree() {
		//assert(vertices.size() > 0);
		//kdTree.reset();
		if (vertices.size() > 0) {
			LOG(INFO)<<"Building kd-tree";
		
			kdTree = unique_ptr<KDTreeND>(new KDTreeND(DIM+1, vertices, 10));
			kdTree->index->buildIndex();
		}
	}

	// These return closest vertices (indices)
	virtual size_t findClosestVertex(const VectorN& queryPoint) const { 
		return findNClosestVertices(queryPoint, 1) [0] ;	
	}

	virtual vector<size_t> findNClosestVertices(const VectorN& queryPoint, const int n) const { 
		const size_t numResults = n;
		std::vector<size_t>   retIndices(numResults);
		std::vector<Float> sqrDistances(numResults);

		nanoflann::KNNResultSet<Float> resultSet(numResults);

		resultSet.init(&retIndices[0], &sqrDistances[0] );
		kdTree->index->findNeighbors(resultSet, &queryPoint[0], nanoflann::SearchParams(10));

		return retIndices;
	}


	virtual Animation<DIM>* makeCopy()=0;

	// Explicit copy method
	virtual void copy(Animation<DIM>* a) {
		a->frames = frames;
		a->neighbours = neighbours;
		a->normals = normals;
		a->vertices = vertices;
		a->timeScaling = timeScaling;
		a->dt = dt;
		a->frameOffsets = frameOffsets;		
		a->frameToGlobalMap = frameToGlobalMap;
		a->globalToFrameMap = globalToFrameMap;
		a->resampledVertices = resampledVertices;
	}


	virtual bool isVertexOnStartFrame(GlobalIndex index) const {
		auto it = globalToFrameMap.find(index);
		return (it->second.first == 0);
	}

	virtual bool isVertexOnEndFrame(GlobalIndex index) const {
		auto it = globalToFrameMap.find(index);
		return (it->second.first == frames.size() - 1);
	}

	virtual Float getDepth(Float frameIndex) const {
		return timeScaling * frameIndex; 
	}

	// Returns the index of the line feature that has vertex with index
	virtual pair<int,int> getLineFeature(int index) const {
		for (int i = 0; i < lineFeatures.size(); i++) {
			for (int j = 0; j < lineFeatures[i].size(); j++){
				if (lineFeatures[i][j] == index)
					return pair<int,int> (i, j);
			}				
		}
		return pair<int,int> (-1,-1);
	}
	
	// Is the vertex with GlobalIndex index part of a line feature?
	virtual bool isLineFeature(int index) const {		
		for (int i = 0; i < lineFeatures.size(); i++) {
			for (int j = 0; j < lineFeatures[i].size(); j++){
				if (lineFeatures[i][j] == index)
					return true;
			}
		}
		return false;
	}

	// Returns signed distance of a point using normal 
	// @deprecated
	virtual Float getSignedDistance(VectorN& point) const {
		Correspondence<DIM> corr = findClosestPoint(point);
		Float d = (corr.p - point).norm();
		VectorN v = (corr.p - point) / d;
		if (v.dot(corr.n) > 1e-5)
			return d;
		return -d;
	}


	// Create vertices in space-time (global indexing)	
	virtual void createGlobalIndexing() 
	{	
		vertices.clear();

		//frameToGlobalMap.reserve(frames.size() * frames[0]->vertices.size());
		//globalToFrameMap.reserve(frames.size() * frames[0]->vertices.size());
		GlobalIndex globalIndex = 0;
		for (int i = 0; i < frames.size(); i++) {
			Float z = getDepth(i);
			const Mesh<DIM>* frame = frames[i].get();
			for (int j = 0; j < frame->vertices.size(); j++) {
				//Vector3 v (frame.vertices[j].pos[0], frame.vertices[j].pos[1], z);
				VectorN v;
				for (int k = 0; k < DIM; k++)
					v[k] = frame->vertices[j].pos[k];

				v[DIM] = z;
				vertices.emplace_back(v);
				FrameIndex frameIndex(i, j);
				frameToGlobalMap[frameIndex] = globalIndex;
				globalToFrameMap[globalIndex] = frameIndex;
				globalIndex++;
			}
		}

		frameOffsets.clear();
		int sum = 0;
		for (int i = 0; i < frames.size(); i++) {
			frameOffsets.push_back(sum);
			sum+= frames[i]->vertices.size();
		}
	}

	// Compute time bounds for each frame and store in minTime/maxTime
	virtual Float computeTimeBounds() {		
		minTime.clear();
		maxTime.clear();
		Float absMaxTime = 0.0;

		for (int i = 0; i < frames.size(); i++) {
			minTime.push_back(1e6);
			maxTime.push_back(-1e-6);
		}

#ifdef ENABLE_PARALLEL
		tbb::parallel_for(size_t(0), size_t(frames.size()), size_t(1), [&](int i)		
		//for (int i = 0; i < frames.size(); i++)
#else
		
#endif
		{
			for (int v = 0; v < frames[i]->vertices.size(); v++) {
				GlobalIndex gIndex = frameToGlobalMap[FrameIndex(i,v)];
				Float curTime = vertices[gIndex][DIM];
				if (curTime < minTime[i]) {
					minTime[i] = curTime;
				}
				if (curTime > maxTime[i]) {
					maxTime[i] = curTime;
				}

			}			
		}
#ifdef ENABLE_PARALLEL
		);
#endif
	
		for (int i = 0; i < frames.size(); i++) {
			VLOG(1)<<"Frame: "<< i << " Time: " << getDepth(i) << " Min: " << minTime[i] << " Max: " << maxTime[i];


			if (maxTime[i] > absMaxTime)
				absMaxTime = maxTime[i];
		}

		return absMaxTime;
	}

	// Perturbs any vertices that lie on the exact time slice
	virtual void perturbVertices() {
		Float eps = 1e-10;
		for (int i = 0; i < frames.size(); i++) {
			Float time = getDepth(i);
			for (int v = 0; v < vertices.size(); v++) {
				if (std::abs(vertices[v][DIM] - time) < eps) {
					Float perturbation = -1e-5;
					if (i == frames.size()-1) {
						perturbation *= -1;
					}					
					vertices[v][DIM] += perturbation;
				}
			}
		}
	}

	// Translate the spacetime mesh 
	virtual void translate(VectorN vec) {
		tbb::parallel_for(size_t(0), size_t(vertices.size()), size_t(1), [&](int i)	{
			vertices[i] += vec;
		}
		);
	}

	virtual bool writeNormalsToDisk(string filename) {

		for (int i = 0; i < frames.size(); i++) {
			std::fstream f;
			char fullPath[2048];
			sprintf_s(fullPath, 2048, "%s_%04d.normals", filename.c_str(), i);
			
			f.open(fullPath,ios::out);
			Json::StyledWriter writer;
			Json::Value root;   // will contains the root value after parsing.
			Json::Value jsonNodes;

			for (int j = 0; j < frames[i]->vertices.size(); j++) {
				int gIndex = frameToGlobalMap.at(FrameIndex(i,j));

				VectorN normal = normals[gIndex];

				Json::Value jsonNormal;
							
				jsonNormal["spaceIndex"] = j;			

				Json::Value arrayNode;
				for (int k=0; k < DIM+1; k++) {
					arrayNode.append(normal[k]);
				}
				jsonNormal["normal"] = arrayNode;
				jsonNodes.append(jsonNormal);
			}

			root["normals"] = jsonNodes;
			std::string outputString = writer.write(root);
			f << outputString << endl;
			f.close();
		}

		return true;
	}

	virtual void setRotationX(Float rotX) {
		rotX = rotX;
	}

	virtual void setTranslation(const VectorN& t) {
		translateVec = t;
	}


	virtual void smooth(int numPasses) {
		vector<VectorN> tmpVertices;
		unsmoothedVertices.reserve(vertices.size());		
		tmpVertices.reserve(vertices.size());

		unsmoothedVertices = vertices;
		tmpVertices = vertices;

		for (int k = 0; k < numPasses; k++) {
			tbb::parallel_for(size_t(0), vertices.size(), size_t(1), [&](int i)
			{
				VectorN sumVector = vertices[i];

				int numNeighbours = 0;				
				for (int j = 0; j < neighbours[i].size(); j++) {
					GlobalIndex nbIndex = neighbours[i][j].index;
					if (nbIndex != Globals::InvalidNeighbour) {
							numNeighbours++;
							sumVector += vertices[nbIndex];												
					}
				}

				if (numNeighbours > 0)
					tmpVertices[i] = (sumVector / numNeighbours);
			});			
			vertices = tmpVertices;
		}
	}

	virtual void unsmoothVertices() {
		vertices = unsmoothedVertices;
		createKDTree();
	}

	vector< shared_ptr <Mesh<DIM>> > frames;
	vector < VectorN > vertices;
	vector < VectorN > unsmoothedVertices;
	tbb::concurrent_hash_map<int, bool> topoChangeVertices;
	
	unique_ptr<KDTreeND> kdTree;

	vector < VectorN > normals;
	vector< vector< Neighbour> > neighbours;
	map<GlobalIndex, FrameIndex> globalToFrameMap;
	map<FrameIndex, GlobalIndex> frameToGlobalMap;

	int numFrames;

	// Time bounds for each frame
	vector<Float> minTime;
	vector<Float> maxTime;

	Float timeScaling; 
	Float dt;
	
	string path;
	LineFeatures lineFeatures;	
	
	Float rotX;
	VectorN translateVec;

	vector<EdgeLabels> frameEdgeLabels;  // Edge labels for each frame of data
	unordered_map< GlobalIndex, Bary > resampledVertices;	
	vector<int> frameOffsets;

	//unordered_map<FrameIndex, GlobalIndex, FrameIndexHasher> frameToGlobalMap;
	//unordered_map<GlobalIndex, FrameIndex> globalToFrameMap;

	private:
	// Do not allow accidental copies because copying is expensive
	Animation(const Animation<DIM>& a);

};

// 2D Animations
class Animation2D : public Animation<2> {

public:
	Animation2D() : Animation<2>() {

	}

	virtual ~Animation2D() {

	}

	virtual Animation2D* makeCopy() {
		Animation2D* anim = new Animation2D();
		Animation<2>::copy(anim);
		return anim;
	}

	// 2D version of base class methods

	virtual bool loadMeshes(const string path,const int startFrame, const int endFrame); 
	virtual Correspondence<2> findClosestPoint(const VectorN&) const ;
	virtual VectorN computeNormalAtVertex(int i) ;
	virtual Correspondence<2> findClosestPointAroundVertex(const Vector3& queryPoint, int index) const ;

	virtual void computeTimeScaling();		
	virtual void findNeighbours();
	virtual void computeNormals();
	virtual void identifyTopologyChanges();

	virtual void setupExtraction();
	virtual void extractSlice(string outFilename, Float time, bool perturbVertices=false);
	

	// Methods specific to Animation2D (helpers)
	int getNextSpaceIndex(GlobalIndex index) const;
	int getPrevSpaceIndex(GlobalIndex index) const;
};

// 3D Animations
class Animation3D : public Animation<3> {

public:	
	Animation3D() : Animation<3>() {

	}

	virtual ~Animation3D() {

	}

	virtual Animation3D* makeCopy() {
		Animation3D* anim = new Animation3D();
		Animation<3>::copy(anim);
		return anim;
	}

	// 3D versions of base class methods
	virtual bool loadMeshes(const string path,const int startFrame, const int endFrame); 
	virtual Correspondence<3> findClosestPoint(const VectorN&) const ;	
	virtual VectorN computeNormalAtVertex(int i) ;	
	
	virtual void computeTimeScaling();		
	virtual void findNeighbours();
	virtual void computeNormals();
	virtual void identifyTopologyChanges();
	virtual void setupExtraction();
	virtual void extractSlice(string outFilename, Float time, bool perturbVertices=false);

	
	virtual void debugResamplingEvents();	
	virtual bool doesFutureTriangleExist(int time, int faceIndex, vector<Geometry::Vertex4>& curTriangle, vector<Geometry::Vertex4>& futureTriangle);
	
	Correspondence<3> findClosestPointOnTriangleFan(const VectorN& p, int curIndex) const;
	Correspondence<3> findClosestPointOnTetrahedron(const VectorN& p, const VectorN& a, const VectorN& b, const VectorN& c, const VectorN& d) const;
	void rotateX(Float);

};