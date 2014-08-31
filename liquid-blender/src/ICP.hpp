#pragma once 

#include "stdafx.h"
#include "Animation.h"
#include "Utility.hpp"
#include "Configuration.hpp"
#include "Output.hpp"
#include "zlib.h"
#include <Bench/BenchTimer.h>

#include <tbb/mutex.h>
#include <cmath>
#include "gtest/gtest.h"
#include <Eigen/SparseCholesky>
#include <Eigen/PardisoSupport>
#include <Eigen/IterativeLinearSolvers>
//#include <Eigen/SparseQR>

#include <Eigen/OrderingMethods>
#include <Eigen/SparseQR>

#include "..\LaplaceSolver\LaplaceSolver.h"

#define USE_NORMAL_LSQR

struct Error {
	Float avgSmooth;
	Float avgPlane;
	Float avgPoint;
	Float maxSmooth;
	Float maxPlane;
	Float maxPoint;
	int numUnmatchedVertices;
	
	Error() {
		avgSmooth = 0;
		avgPlane = 0;
		avgPoint = 0;
		maxSmooth = 0;
		maxPlane = 0;
		maxPoint = 0;
		numUnmatchedVertices = 0;
	}
};

struct VertexError {

	VertexError(int _index, Float _error) : index(_index), error(_error) {
	}

	bool operator >(const VertexError& ve) {
		return error > ve.error;
	}

	bool operator <(const VertexError& ve) {
		return error < ve.error;
	}

	int index;
	Float error;
};

struct NeighbourNode : public LaplaceNeighbour {
	NeighbourNode(int nodeIndex, Float weight) : nodeIndex(nodeIndex), weight(weight) {

	}

	NeighbourNode() {}
	int getIndex() { return nodeIndex;}	
	Float getWeight() { return 1.0;}

	int nodeIndex;
	Float weight;
};

template <int DIM>
class ICP {
public:
	typedef Eigen::Matrix<Float, DIM+1, 1> VectorN;
	typedef Eigen::Matrix<Float, DIM+1, DIM+1> MatrixN;

	typedef Eigen::Triplet<Float> TripletF;

	ICP(const Configuration<DIM> & cfg): cfgName(cfg.configName), animA(cfg.animA), animB(cfg.animB), 
		samplingRate(cfg.samplingRate), 
		subsamplingRate(cfg.subsamplingRate), 
		pointCorrs(cfg.pointCorrs), 
		lineCorrs(cfg.lineCorrs), 
		spaceCorrs(cfg.spaceCorrs),
		maxIterations(cfg.maxIterations), 
		saveEvery(cfg.saveEvery),
		initialWeights(cfg.weights),
		configPath(cfg.configPath),
		enforceEqualLength(cfg.enforceEqualLength),
		usePreviousSolution(cfg.usePreviousSolution),
		prevSolPath(cfg.prevSolPath),
		enforceHardConstraints(cfg.enforceHardConstraints),
		generateProjectedOutput(cfg.generateProjectedOutput)
	{
		
		reset();
	}


	ICP(const unique_ptr<Animation<DIM>>& animA, const unique_ptr<Animation<DIM>>& animB, const int samplingRate, const int subsamplingRate, 
		const int maxIterations,
		const int saveEvery,
		const PointCorrs& pointCorrs, 
		const LineCorrs& lineCorrs, 	

		const ICPWeights weights=ICPWeights(),
		bool enforceEqualLength = true, 		
		bool enforceHardConstraints = false): 
	animA(animA), animB(animB), samplingRate(samplingRate), subsamplingRate(subsamplingRate), pointCorrs(pointCorrs), lineCorrs(lineCorrs),
		initialWeights(weights),
		enforceEqualLength(enforceEqualLength),
		maxIterations(maxIterations),
		saveEvery(saveEvery),
		enforceHardConstraints(enforceHardConstraints),
		generateProjectedOutput(generateProjectedOutput)
	{
		LOG(INFO)<< "Anim A: #Vertices: " << animA->vertices.size() << endl;
		LOG(INFO)<< "Anim B: #Vertices: " << animB->vertices.size() << endl;
	}


	void reset() {
		// Reset weights
		curWeights = initialWeights;
		vertexToNodeIndexMap.clear();
	}

	void initialize() {

		energyFile.open("energy.txt", ios::out);
		enableESmooth = true;
		enableEPoint = true;
		enableEPlane = true;
	
		dumpDeformationNodesToDisk = false;
		iteration = 1;
		curAnim = unique_ptr<Animation<DIM>>(animA->makeCopy());
		maxCorrDist = animA->timeScaling * 25 ;
		minCorrAngle = 0.1;
		corrs.clear();
		corrs.reserve(animA->vertices.size());
		tmpCorrs.reserve(animA->vertices.size());

		allVertexList.reserve(animA->vertices.size());


		for (int i = 0; i < animA->vertices.size(); i++) {
			corrs.push_back(Correspondence<DIM>());
			allVertexList.push_back(i);
		}

		setupRadii();

		// Make every pointCorrespondence a node
		// Sample line features with nodes
		createNodesOnUserCorrs();

		// Sample the mesh / create deformation nodes.
		if (usePreviousSolution) {
			if (!readSolution(prevSolPath)) {
				LOG(ERROR) << "Unable to read previous solution!";
			}
		}
		else {
			createNodes();
		}
		
		
		// Compute node smoothness weights and vertex weights		
		precomputeWeights();

		// Detect vertices that are a part of a small component and do not have user correspondences
		flagTinyComponentVertices();

		if (!usePreviousSolution) {
			// Find initial correspondences for lines
			findInitialCorrespondences();

		}
		
		// Diffuse correspondences (TODO: only for new correspondences?)
		diffuseCorrespondences();

		// Determine free variables and create a mapping 
		computeNodeOffsets();

		if (!enforceHardConstraints) 
			lineCorrs.clear();

	}
	
	void flagTinyComponentVertices() {
		tbb::parallel_for(size_t(0), animA->vertices.size(), size_t(1), [&](int i)
		//for (int i = 0; i < animA->vertices.size(); i++) 
		{
			FrameIndex fIndex = animA->globalToFrameMap.at(i);
			const Mesh<DIM>* mesh = animA->frames[fIndex.first].get();
			if (mesh->vertices[fIndex.second].isTinyComponent) {
				bool isNearUserCorr = false;
				for (int j = 0; j < vertexWeights[i].size(); j++) {
					if (isNodeUserCorrespondence(vertexWeights[i][j].nodeIndex)) {
						isNearUserCorr = true;
						break;
					}
				}
				if (!isNearUserCorr) {
					tbb::concurrent_hash_map<int, bool>::accessor a;
					tinyComponentVertices.insert(a, pair<int, bool>(i, true));	
					
					//LOG(INFO) << "Found tiny component vertex:" << fIndex.first << "," << fIndex.second;
				}
			}
		});
	}

	// Computes the offsets for each node in the list of free variables and also computes numFreeVariables (i.e. numColumns)
	void computeNodeOffsets() {
		nodeOffsets.clear();
		numFreeVariables = 0;
		for (int i = 0; i < nodes.size(); i++) {
			nodeOffsets.push_back(numFreeVariables);
			GlobalIndex index = nodes[i];
			if (enforceHardConstraints && isPointCorrespondence(index)) 
				numFreeVariables += (DIM+1)*(DIM+1);
			else if (isSpaceCorrespondence(index))
				numFreeVariables += (DIM+1)*(DIM+1) + 1;
			else if (enforceEqualLength && (isNodeOnStartFrame(i) || isNodeOnEndFrame(i))) 
				numFreeVariables += ((DIM+1)*(DIM+1) + DIM);
			else
				numFreeVariables += ((DIM+1)*(DIM+1) + DIM+1);

			
		}		
		numCols = numFreeVariables;
	}

	// Attempt at parallel sampling - doesn't seem to work correctly at the moment (DO NOT USE!)
	void samplePoissonDiscParallel(map<int,bool>& indexMap, int numNodes, Float maxRadius) {
		boost::random::mt19937 gen;
		int numLeft = indexMap.size();
		boost::random::uniform_int_distribution<> dist(0, indexMap.size());
		const int maxTrials = 2;

		unordered_set<int> curSamples;

		typedef tbb::spin_mutex MutexType;
		MutexType curSamplesMutex;

		tbb::parallel_for(tbb::blocked_range<int>(0,numNodes), [&] (tbb::blocked_range<int>& r) {
			for (int i = r.begin(); i < r.end(); i++) {


				// Generate a random sample			
				bool validDart = false;

				int numTrials = 0;
				int dartIndex = 0;

				while (!validDart && numTrials < maxTrials && numLeft > 0) {
					dartIndex = dist(gen);

					LOG(INFO)<<"Generating dart:" << dartIndex;

					// i.e. if indexMap[dartIndex] is true, we have a valid position
					if (indexMap[dartIndex]) {
						validDart = true;

						{														
							for (auto it = curSamples.begin(); it != curSamples.end(); it++) {
								if ((animA->vertices[dartIndex] - animA->vertices[*it]).squaredNorm() < maxRadius * maxRadius) {
									validDart = false;
									break;
								}
							}
						}
						
						if (validDart)
						{
							{
								MutexType::scoped_lock lock(curSamplesMutex);
								curSamples.insert(dartIndex);
							}
							break;					
						}
					}

					numTrials++;
				}

				if (validDart) {
						
					createNode(dartIndex);

					// Now run BFS and remove all nodes within maxRadius
					auto neighbours = doBFS(animA, dartIndex, maxRadius);

					for (int j = 0; j < neighbours.size(); j++) {
						indexMap[neighbours[j].index] = false;
						numLeft--;
					}				
					
					// Mutex
					{
						MutexType::scoped_lock lock(curSamplesMutex);
						curSamples.erase(dartIndex);
					}
				}
			
			}
		});
	}

	// Samples the spacetime surface animA using dart throwing.
	// TODO: This is not fast and should be parallelized.
	void samplePoissonDisc(map<int,bool>& indexMap, int numNodes, Float maxRadius) {

		boost::random::mt19937 gen;

		for (int i = 0; i < numNodes; i++) {
			// Generate new dart
			if (indexMap.size() <= 0) return;


			boost::random::uniform_int_distribution<> dist(0, indexMap.size()-1);
			int dartIndex = dist(gen);

			int count = 0;
			auto it = indexMap.begin();
			for (; count < dartIndex; it++, count++);			
			int vertexIndex = it->first;

			VLOG(1) << "Dart index:" << dartIndex;
			VLOG(1) << "Sample:" << i <<  " creating nodes at index:" << vertexIndex ; 
			Eigen::BenchTimer timer;



			if (!isNode(vertexIndex)) {
				createNode(vertexIndex);

				// Now run BFS and remove all nodes within maxRadius
				auto neighbours = doBFS(animA, vertexIndex, maxRadius);
				
				//timer.reset();
				//timer.start();

				int numRemoved = 0;
				for (int j = 0; j < neighbours.size(); j++) {
					auto it = indexMap.find(neighbours[j].index);
					if (it != indexMap.end()) {
						indexMap.erase(it);
						numRemoved++;
					}					
				}

				//timer.stop();
				//VLOG(2) << "Delete Time:" << timer.total();
				//VLOG(2) <<" # neighbours in radius:" << neighbours.size() << " # vertices removed:" << numRemoved << " # left:" << indexMap.size();
			}
		}
		
	}

	void createNodes() {
		LOG(INFO) << "Creating nodes";
		map<int, bool> allVerticesMap;
		map<int, bool> endVerticesMap;



		int numVertices = animA->vertices.size();
		int numNodes = numVertices / samplingRate;

		VLOG(1) << "Number of desired nodes:" << numNodes;

		Float maxRadius = (2.414) * animA->timeScaling * pow(samplingRate,0.3333);

		VLOG(1) << "BFS Radius: " << maxRadius;

		// Create vertex lists 
		for (int i = 0; i < animA->vertices.size(); i++) {
			FrameIndex frameIndex = animA->globalToFrameMap.find(i)->second;
			if (frameIndex.first == 0 
				|| frameIndex.first == animA->frames.size()-1
				) {
				endVerticesMap[i] = true;
			}		
			allVerticesMap[i] = true;
		}

		// Remove all vertices near user nodes
		for (int i = 0; i < nodes.size(); i++) {		
			vector<Neighbour> neighbours = doBFS(animA, nodes[i], maxRadius);

			for (int j = 0; j < neighbours.size(); j++) {
				if (allVerticesMap.find(neighbours[j].index) != allVerticesMap.end())
					allVerticesMap.erase(allVerticesMap.find(neighbours[j].index));
			}

			for (int j = 0; j < neighbours.size(); j++) {
				if (endVerticesMap.find(neighbours[j].index) != endVerticesMap.end())
					endVerticesMap.erase(endVerticesMap.find(neighbours[j].index));
			}

		}

		// Sample first and last frames
		samplePoissonDisc(endVerticesMap, numNodes, maxRadius);

		VLOG(1)<<"# nodes on start and end frames:" << nodes.size();

		// Remove nearby vertices
		for (int i = 0; i < nodes.size(); i++) {
			if (allVerticesMap.find(nodes[i]) != allVerticesMap.end()) {
				vector<Neighbour> neighbours = doBFS(animA, nodes[i], maxRadius);

				for (int j = 0; j < neighbours.size(); j++) {
					if (allVerticesMap.find(neighbours[j].index) != allVerticesMap.end())
						allVerticesMap.erase(allVerticesMap.find(neighbours[j].index));
				}
			} 
		}

		numNodes -= nodes.size();
		

		//numNodes *= 2;
		// Sample the rest of the domain
		VLOG(1) << "Number of desired nodes:" << numNodes;
		samplePoissonDisc(allVerticesMap, numNodes, maxRadius);

		LOG(INFO) << "# created nodes:" << nodes.size();

	}	

	void writeDeformationNodes(string path) {
		fstream f;
		f.open(path + ".nodes",ios::out);
		
		Json::StyledWriter writer;
		Json::Value root;   // will contains the root value after parsing.
		Json::Reader reader;

		Json::Value jsonNodes;


		for (int i = 0; i < nodes.size(); i++) {
			FrameIndex fIndex = animA->globalToFrameMap.at(nodes[i]);
			Json::Value deformationNode;
			deformationNode["index"] = i;
			deformationNode["timeIndex"] = fIndex.first;
			deformationNode["spaceIndex"] = fIndex.second;
			
			Json::Value vertices;

			for (int j = 0; j < nodeToVertexWeights[i].size(); j++) {
				Json::Value vertex;
				FrameIndex tmpIndex = animA->globalToFrameMap.at(nodeToVertexWeights[i][j].nodeIndex);
				vertex["timeIndex"] = tmpIndex.first;
				vertex["spaceIndex"] = tmpIndex.second;
				vertex["weight"] = nodeToVertexWeights[i][j].weight;
				vertices.append(vertex);				
			}
			nodeToVertexWeights[i].clear();
			deformationNode["vertices"] = vertices;
			jsonNodes.append(deformationNode);
		}

		nodeToVertexWeights.clear();
		
		root["nodes"] = jsonNodes;
		string outputString = writer.write(root);
		f << outputString << endl;
		f.close();
	}

	void sampleUniformTime();	
	void createNodesOnUserCorrs();

	// Computes weight using a quintic kernel
	Float computeWeight(Float d, Float r) {
		if (d > r) 
			return 0.0;
		return max(Float(0), Float(pow((1 - (d*d)/(r*r)),3)));
	}

	void setupRadii();
	
	// Compute contribution of nodes to each vertex, as well as to other nodes
	void precomputeWeights() {
		tbb::mutex simpleMutex;

		// Create vectors for weights
		nodeWeights.reserve(nodes.size());
		for (int i = 0; i < nodes.size(); i++) {
			vector<LaplaceNeighbour*> nbs;
			nbs.reserve(30);
			nodeWeights.push_back(nbs);
		}

		vertexWeights.reserve(animA->vertices.size());
		for (int i = 0; i < animA->vertices.size(); i++) {
			vector<NeighbourNode> nbs;
			nbs.reserve(30);
			vertexWeights.push_back(nbs);
		}
		
		maxNumNeighbours = 0;
		tbb::parallel_for(size_t(0), nodes.size(), size_t(1), [&](int i)
		//for (int i = 0; i < nodes.size(); i++)
		{
			if (dumpDeformationNodesToDisk) 
				nodeToVertexWeights.push_back( vector<NeighbourNode>() );

			vector< Neighbour > nbVertices = doBFS(animA, nodes[i], nodeRadius);
			Float smoothnessSum = 0.0;
			
			for (auto nb : nbVertices) {
				Float vertexWeight = computeWeight(nb.d, vertexRadius);

				if (nb.d <= vertexRadius && vertexWeight > 0.01) {
					tbb::mutex::scoped_lock lock(simpleMutex);					
					vertexWeights[nb.index].push_back(NeighbourNode(i,vertexWeight));

					if (dumpDeformationNodesToDisk) 
						nodeToVertexWeights[i].push_back(NeighbourNode(nb.index,vertexWeight));
				}

				// Don't include self for node smoothness
				if (nb.index != nodes[i] && isNode(nb.index)) {
					//Float nodeWeight = 1.0;
					Float nodeWeight = computeWeight(nb.d, nodeRadius);
					nodeWeights[i].push_back(new NeighbourNode(vertexToNodeIndexMap[nb.index],nodeWeight));
					smoothnessSum += nodeWeight;
				}
			}
			
			FrameIndex f = animA->globalToFrameMap.at(nodes[i]);
			//LOG(INFO)<<"Node: " << f.first << "," << f.second << " # nbs: " << nodeWeights[i].size();

			// Normalize node smoothness weights			
			for (int j = 0; j < nodeWeights[i].size(); j++) {
				NeighbourNode* nb = dynamic_cast<NeighbourNode*>(nodeWeights[i][j]);
				nb->weight = sqrt(nb->weight/smoothnessSum);
				FrameIndex fIndex = animA->globalToFrameMap.at(nodes[nb->nodeIndex]);
				//LOG(INFO)<<"Neighbour: " << fIndex.first << " " << fIndex.second << " " << nb.weight << endl;
			}

			
			
			
		}
		);

		LOG(INFO) << "Done precomputing node weights" ;
		maxNumNeighbours = 1;


		int avgNumNeighbours = 0;
		
		// Normalize vertex weights
		tbb::parallel_for(size_t(0),  animA->vertices.size(), size_t(1), [&](int i) 
		//for (int i = 0; i < animA->vertices.size(); i++) 
		{

			auto it = animA->globalToFrameMap.find(i);
			
			Float sum = 0.0;
			avgNumNeighbours += vertexWeights[i].size();
			
			for (int j = 0; j < vertexWeights[i].size(); j++) {
				sum += vertexWeights[i][j].weight;
				int index = nodes[ vertexWeights[i][j].nodeIndex ];
				auto it2 = animA->globalToFrameMap.find(index);

			}

			if (vertexWeights[i].size() == 0)
				LOG(ERROR) << "Vertex:" << i << " # nodes:" << vertexWeights[i].size() << endl;
			for (int j = 0; j < vertexWeights[i].size(); j++) {
				vertexWeights[i][j].weight /= sum;
				
			}
									
		}		
		);
		
		for (int i = 0; i < animA->vertices.size(); i++) {
			maxNumNeighbours = max(int(vertexWeights[i].size()), int(maxNumNeighbours));
		}
		
		//f.close();
		//maxNumNeighbours = avgNumNeighbours / animA->vertices.size();
		

	}

	void widenLines();

	void findInitialCorrespondences() {

		// Widen lines 
		widenLines();

		for (int i = 0; i < nodes.size(); i++) {			

			if (isNodeUserCorrespondence(i)) 
			{
				if (isSpaceCorrespondence(nodes[i])) 
				{
					b[i] = spaceCorrs[nodes[i]];
					VLOG(1) << "User Space Corr Node:" << i << " " << nodes[i] << " " << b[i] ;		
				}
				else if (isPointCorrespondence(nodes[i]))
				{
					b[i] = animB->vertices[ pointCorrs[nodes[i]] ] - animA->vertices[ nodes[i] ];					
					VLOG(1) << "User Point Corr Node:" << i << " " << nodes[i] << " " << b[i] ;		
				}

				

				Float d = b[i].norm() ;
				if (d > maxCorrDist) {
					maxCorrDist = d;
				}
			}
			else if (isNodeOnStartFrame(i)) {
				b[i][DIM] = 0;
			}
			else if (isNodeOnEndFrame(i)) {
				b[i][DIM] = animB->getDepth(animB->frames.size()) - animA->vertices[nodes[i]][DIM]; 
			}
			else {
				//Globals.getDepth(time * len(self.anims[1].frames) / len(self.anims[0].frames)) - Globals.getDepth(time)  
				auto it = animA->globalToFrameMap.find(nodes[i]);
				int timeA = it->second.first;
				int spaceA = it->second.second;
				for (int j = 0; j < DIM; j++)
					b[i][j] = 0;

				Float timeB = 1.0 * timeA * ( (1.0 * animB->frames.size()) /(1.0 * animA->frames.size()) );
				b[i][DIM] = animB->getDepth(timeB) - animA->vertices[nodes[i]][DIM];
				
			}
			
		}		


	}

	void diffuseCorrespondences() {
		runSteadyStateDiffusion();
		//runExplicitDiffusion(5000);

		VLOG(1) << "Diffused:";
		for (int i = 0; i < nodes.size(); i++) {
			if (isNodeUserCorrespondence(i)) {
				for (int j = 0; j < nodeWeights[i].size(); j++)
					VLOG(1) << b[nodeWeights[i][j]->getIndex()].format(CommaFmt);
			}
			
		}
	}

	bool isNodeConstrained(int nodeIndex, Float& value, int dim=0) {
		GlobalIndex gIndex = nodes[nodeIndex];
		if (isPointCorrespondence(gIndex) ) {
			value = b[nodeIndex][dim];
			return true;
		}
		else if ((isNodeOnStartFrame(nodeIndex) || isNodeOnEndFrame(nodeIndex))  && dim == DIM) {
			value = b[nodeIndex][dim];
			return true;
		}

		return false;
	}

	// Solves the diffusion Laplace equation to steady state
	void runSteadyStateDiffusion() {
		for (int dim = 0; dim < DIM+1; dim++) {
			using namespace std::placeholders;
			LaplaceSolver<Float> solver(nodes.size(), std::bind(&ICP<DIM>::isNodeConstrained, this, _1,_2, dim), nodeWeights);

			VectorX sol = solver.solve();

			// Copy solution 
			for (int i = 0; i < nodes.size(); i++) {
				b[i][dim] = sol[i];
			}
		}
		
	}

	void runExplicitDiffusion(const int numSteps) {



		// Make a copy
		auto tmpB = b;

		Float diffusionTimestep = 0.01;		

		for (int step = 0; step < numSteps; step++) {

			//for (int i = 0; i < nodes.size(); i++) {
			tbb::parallel_for(size_t(0), nodes.size(), size_t(1), [&](int i) 
			{				
				FrameIndex fIndexI = animA->globalToFrameMap.at(nodes[i]);
				int timeI = fIndexI.first;

				VectorN sum = VectorN::Zero();
				Float sumWeights = 0.0;
				for (int j = 0; j < nodeWeights[i].size(); j++) {				
					int gIndexJ = nodes[nodeWeights[i][j]->nodeIndex];
					FrameIndex fIndexJ = animA->globalToFrameMap.at(gIndexJ);
					Float scaling = 1.0;

					if (fIndexJ.first == timeI)
						scaling = 1.0;

					Float weight = scaling * nodeWeights[i][j]->weight;

					sum += weight * (b[nodeWeights[i][j]->nodeIndex] - b[i]);										
				}

				tmpB[i] = sum / sumWeights ;
			});
			
			//for (int i = 0; i < nodes.size(); i++) {
			tbb::parallel_for(size_t(0), nodes.size(), size_t(1), [&](int i) 
			{
				if (isNodeUserCorrespondence(i));
				else if (enforceEqualLength && (isNodeOnStartFrame(i) || isNodeOnEndFrame(i)) ) {
						for (int k = 0; k < DIM; k++)
							b[i][k] += tmpB[i][k];					
				}
				else
					b[i] += tmpB[i];				
			});
		}


	}


	// Creates a node given GlobalIndex index on animA
	bool createNode(int index) {
		if (isNode(index)) 
			return false;

		vertexToNodeIndexMap[index] = nodes.size();
		nodes.push_back(index);

		MatrixN mat;
		mat.setIdentity();
		VectorN v = VectorN::Zero();

		A.push_back(mat);
		b.push_back(v);
		

		//LOG(INFO)<<"Created node at:"<<animA->vertices[index].format(CommaFmt);		
		return true;
	}

	ICPWeights& getWeights() {
		return curWeights;
	}

	VectorN computeEFitPointVector(GlobalIndex globalIndex) {
		VectorN e = curWeights.weightPoint * (curAnim->vertices[globalIndex] - corrs[globalIndex].p);

		return e;
	}

	Float computeEFitPlane(GlobalIndex globalIndex) {
		//Float e = curWeights.weightPlane * abs(corrs[globalIndex].n.dot(curAnim->vertices[globalIndex] - corrs[globalIndex].p));

		Float e = curWeights.weightPlane * (corrs[globalIndex].n.dot(curAnim->vertices[globalIndex] - corrs[globalIndex].p));

		return e;
	}

	// Note: Computes energy for node with nodeIndexI relative to its neighbour index j (i.e. jth element in its neighbour list)
	VectorN computeESmoothVector(int nodeIndexI, int j) {
		NeighbourNode* nbNodeJ = dynamic_cast<NeighbourNode*>(nodeWeights[nodeIndexI][j]);
		int nodeIndexJ = nbNodeJ->nodeIndex;
		Float nodeWeight = nbNodeJ->weight;

		VectorN xI = animA->vertices[nodes[nodeIndexI]];
		VectorN xJ = animA->vertices[nodes[nodeIndexJ]];

		VectorN e =  curWeights.weightSmooth * nodeWeight * (A[nodeIndexI]*(xJ-xI) + xI + b[nodeIndexI] - xJ - b[nodeIndexJ]);

		return e;
	}
	
	void computeJFitPoint(GlobalIndex globalIndex) {
		VectorN term = curAnim->vertices[globalIndex] - corrs[globalIndex].p;						VectorN v = animA->vertices[globalIndex];

		for (int i = 0; i < vertexWeights[globalIndex].size(); i++) {
			int nodeIndex = vertexWeights[globalIndex][i].nodeIndex;
			Float nodeWeight = vertexWeights[globalIndex][i].weight;
			int nodeOffset = nodeOffsets[nodeIndex];
			
			for (int row = 0; row < DIM+1; row++) {

				for (int col = 0; col < DIM+1; col++) {				
					addToJacobian(jFitPointMap[globalIndex] + row, nodeOffset + (DIM+1)*col + row, 
						curWeights.weightPoint  * nodeWeight * (animA->vertices[globalIndex][col] - animA->vertices[nodes[nodeIndex]][col]));				
				}

				int numRows = DIM+1;
				int gIndexNode = nodes[nodeIndex];
				if (isPointCorrespondence(gIndexNode))  {					
					if (enforceHardConstraints) {
						numRows = 0; // Completely constrained
					} 
				}
				else if (isSpaceCorrespondence(gIndexNode)) {
					numRows = 1; // Only allow time to change
				}
				else {					
					if (enforceEqualLength && (isNodeOnStartFrame(nodeIndex) || isNodeOnEndFrame(nodeIndex))) 
						numRows = DIM;
				}

				if (row < numRows) {
					addToJacobian(jFitPointMap[globalIndex] + row, nodeOffset + (DIM+1)*(DIM+1) + row, curWeights.weightPoint * nodeWeight);
				}

			}			
		}
	}

	void computeJFitPlane(GlobalIndex globalIndex) {
				
		VectorN normal = corrs[globalIndex].n;
		Float sign = 1;

		//if (normal.dot(curAnim->vertices[globalIndex] - corrs[globalIndex].p) < 0)
		//	sign  = -1;


		for (int i = 0; i < vertexWeights[globalIndex].size(); i++) {
			int nodeIndex = vertexWeights[globalIndex][i].nodeIndex;
			Float nodeWeight = vertexWeights[globalIndex][i].weight;
			int nodeOffset = nodeOffsets[nodeIndex];


			for (int col = 0; col < DIM+1; col++) {
				for (int row = 0; row < DIM+1; row++) {
					addToJacobian(jFitPlaneMap[globalIndex], nodeOffset + (DIM+1)*col + row, 
						sign * curWeights.weightPlane * nodeWeight * normal[row] * 
						(animA->vertices[globalIndex][col] - animA->vertices[nodes[nodeIndex]][col]));
				}
			}

			int numRows = DIM+1;
			int gIndexNode = nodes[nodeIndex];
			if (isPointCorrespondence(gIndexNode))  {					
				if (enforceHardConstraints) {
					numRows = 0; // Completely constrained
				} 
			}
			else if (isSpaceCorrespondence(gIndexNode)) {
				numRows = 1; // Only allow time to change
			}
			else {					
				if (enforceEqualLength && (isNodeOnStartFrame(nodeIndex) || isNodeOnEndFrame(nodeIndex))) 
					numRows = DIM;
			}

			for (int row = 0; row < numRows; row++) {
				addToJacobian(jFitPlaneMap[globalIndex], nodeOffset + (DIM+1)*(DIM+1) + row, sign * curWeights.weightPlane * nodeWeight * normal[row]);
			}

		}

	}

	void computeJSmooth(int nodeIndexI, int j) {

		NeighbourNode* nbNode = dynamic_cast<NeighbourNode*>(nodeWeights[nodeIndexI][j]);
		int nodeIndexJ = nbNode->nodeIndex;
		Float nodeWeight = nbNode->weight;

		VectorN xI = animA->vertices[nodes[nodeIndexI]];
		VectorN xJ = animA->vertices[nodes[nodeIndexJ]];

		VectorN term = A[nodeIndexI] * (xJ-xI) + xI + b[nodeIndexI] - (xJ + b[nodeIndexJ]);
		
		Float netWeight = curWeights.weightSmooth * nodeWeight;
		int nodeOffsetI = nodeOffsets[nodeIndexI];
		int nodeOffsetJ = nodeOffsets[nodeIndexJ];

		// Deal with each row
		for (int k = 0; k < DIM+1; k++) {
			// Handle the A_i
			for (int col = 0; col < DIM+1; col++) {
				addToJacobian(curRow + k, nodeOffsetI + (DIM+1)*col + k, netWeight * (xJ[col] - xI[col]));
			}

			// Handle b_i
			int numRows = DIM+1;
			int gIndexNode = nodes[nodeIndexI];
			if (isPointCorrespondence(gIndexNode))  {					
				if (enforceHardConstraints) {
					numRows = 0; // Completely constrained
				} 
			}
			else if (isSpaceCorrespondence(gIndexNode)) {
				numRows = 1; // Only allow time to change
			}
			else {					
				if (enforceEqualLength && (isNodeOnStartFrame(nodeIndexI) || isNodeOnEndFrame(nodeIndexI))) 
					numRows = DIM;
			}

			if (k < numRows) {
				addToJacobian(curRow + k, nodeOffsetI + (DIM+1)*(DIM+1) + k, netWeight);
			}

			// Handle b_j
			numRows = DIM+1;
			gIndexNode = nodes[nodeIndexJ];
			if (isPointCorrespondence(gIndexNode))  {					
				if (enforceHardConstraints) {
					numRows = 0; // Completely constrained
				} 
			}
			else if (isSpaceCorrespondence(gIndexNode)) {
				numRows = 1; // Only allow time to change
			}
			else {					
				if (enforceEqualLength && (isNodeOnStartFrame(nodeIndexJ) || isNodeOnEndFrame(nodeIndexJ))) 
					numRows = DIM;
			}			

			if (k < numRows) {
				addToJacobian(curRow + k, nodeOffsetJ + (DIM+1)*(DIM+1) + k, -netWeight);
			}
			
		}

	}

	void subsample(vector<int>& subsamples) {
		LOG(INFO) << "Creating samples on the mesh";
		LOG(INFO) << "Attempting to subsample " << animA->vertices.size()/ subsamplingRate << " out of " << animA->vertices.size();
		map<int, bool> indexMap;		

		int numVertices = animA->vertices.size();		
		int numSubsamples = animA->vertices.size() / subsamplingRate;
		
		// Create vertex lists 
		for (int i = 0; i < animA->vertices.size(); i++) {
			indexMap[i] = true;
		}

		map<int,bool> subsampleMap;
		boost::random::mt19937 gen;

		const Float percentageRetained = 0.1;
		int numRetained = percentageRetained * numSubsamples;
		
		if (perVertexErrorList.size() > 0) {
			make_heap(perVertexErrorList.begin(), perVertexErrorList.end());

			for (int i = 0; i < numRetained; i++) {

				VertexError ve = perVertexErrorList.front();				
				subsamples.push_back(ve.index);
				subsampleMap[ve.index] = true;
				pop_heap(perVertexErrorList.begin(), perVertexErrorList.end());
				perVertexErrorList.pop_back();
			}
			
		}
		else {
			numRetained = 0;
		}

		

		// Subsample on a per-frame basis

		int frame = 0;
		for (int i = 0; i < numSubsamples - numRetained; i++) {
			boost::random::uniform_int_distribution<> dist(0, animA->frames[frame]->numVertices - 1);

			bool unique = true;
			int dartIndex;

			do {
				int vertexIndex = dist(gen);
				dartIndex = animA->frameToGlobalMap.at(FrameIndex(frame,vertexIndex));
			} while (subsampleMap.find(dartIndex)!=subsampleMap.end());

			subsamples.push_back(dartIndex);
			subsampleMap[dartIndex] = true;

			frame++;			
			if (frame == animA->numFrames)
				frame = 0;
		}
		LOG(INFO) << "Number of subsampled vertices:" << subsamples.size();
		
		return;

		// // Subsample randomly (global)
		//map<int,bool> subsampleMap;
		//for (int i = 0; i < numSubsamples; i++) {
		//	boost::random::uniform_int_distribution<> dist(0, indexMap.size()-1);

		//	bool unique = true;
		//	int dartIndex;
		//	
		//	do {
		//		dartIndex = dist(gen);
		//	} while (subsampleMap.find(dartIndex)!=subsampleMap.end());

		//	subsamples.push_back(dartIndex);
		//	subsampleMap[dartIndex] = true;

		//}
		//LOG(INFO) << "Number of subsampled vertices:" << subsamples.size();

		//// Geodesic subsampling 
		//Float maxRadius = (2.414) * animA->timeScaling * pow(subsamplingRate,0.3333);
		//LOG(INFO) << "Subsampling Radius: " << maxRadius;
		//for (int i = 0; i < numSubsamples; i++) {
		//	// Generate new dart
		//	if (indexMap.size() <= 0) break;

		//	//LOG(INFO) << "Size of map:" << indexMap.size();
		//	//for (auto it = indexMap.begin(); it != indexMap.end(); it++) {
		//	//	LOG(INFO) << it->first;
		//	//}
		//	boost::random::uniform_int_distribution<> dist(0, indexMap.size()-1);
		//	int dartIndex = dist(gen);

		//	int count = 0;
		//	auto it = indexMap.begin();
		//	for (; count < dartIndex; it++, count++);			
		//	int vertexIndex = it->first;

		//	//LOG(INFO) << "Dart index:" << dartIndex;
		//	//LOG(INFO) << "Sample:" << i <<  " creating nodes at index:" << vertexIndex ; 
		//	subsamples.push_back(vertexIndex);

		//	// Now run BFS and remove all nodes within maxRadius
		//	auto neighbours = doBFS(animA, vertexIndex, maxRadius);

		//	int numRemoved = 0;
		//	for (int j = 0; j < neighbours.size(); j++) {

		//		if (indexMap.find(neighbours[j].index) != indexMap.end()) {
		//			indexMap.erase(indexMap.find(neighbours[j].index));
		//			numRemoved++;
		//		}					
		//	}
		//}
		//LOG(INFO) << "Number of subsampled vertices:" << subsamples.size();

		

	}

	void computeEnergy() {
		curRow = 0;		
		Float normESmooth, normEPlane, normEPoint;
		normESmooth = 0;
		normEPlane = 0;
		normEPoint = 0;
		

		if (enableESmooth) {
			for (int i = 0; i < nodes.size(); i++) {
				for (int j = 0; j < nodeWeights[i].size(); j++) {
					VectorN eSmoothVector = computeESmoothVector(i,j);
					for (int k = 0; k < DIM+1; k++) {
						energy[curRow] = eSmoothVector[k];		
						normESmooth += energy[curRow]*energy[curRow];
						curRow++;

						
					}
				}
			}
		}


		maxEPointNorm = 0;
		maxEPlaneNorm = 0;
		perVertexErrorList.clear();

		for (int c = 0; c < subsamples.size(); c++) {

			int i = subsamples[c];
			if (!corrs[i].valid) 
				continue;

			if (enableEPoint) {
				VectorN eFitPointVector = computeEFitPointVector(i);	
				Float curEPoint = eFitPointVector.norm();
				maxEPointNorm = max(maxEPointNorm,curEPoint);
				
				perVertexErrorList.push_back(VertexError(i, curEPoint));

				for (int k = 0; k < DIM+1; k++) {
					energy[curRow] = eFitPointVector[k];
					normEPoint += energy[curRow]*energy[curRow];
					curRow++;					
				}
				
			}
			if (enableEPlane) {				
				energy[curRow] = computeEFitPlane(i);		
				
				maxEPlaneNorm = max(maxEPlaneNorm, energy[curRow]);

				normEPlane += energy[curRow]*energy[curRow];
				curRow++;

			}
		}


		energyFile << normESmooth << " ,  " << normEPoint << " , " << normEPlane << " , " << normESmooth+normEPoint+normEPlane<<endl;

		prevNormEnergy = curNormEnergy;
		curNormEnergy = normESmooth+normEPoint+normEPlane;
		//LOG(INFO)<<"Energy: ";
		//for (int i = 0; i < numRows; i++) {
		//	LOG(INFO) << energy[i] << " ";
		//}
		//LOG(INFO) << endl;

		//LOG(INFO) << energy << endl;
	}

	void computeJacobian() {				
		curRow = 0;

		if (enableESmooth) {
			for (int i = 0; i < nodes.size(); i++) {
				for (int j = 0; j < nodeWeights[i].size(); j++) {
					computeJSmooth(i,j);			
					curRow += (DIM+1);
				}
			}
		}



#ifdef ENABLE_PARALLEL
		tbb::parallel_for(size_t(0), subsamples.size(), size_t(1), [&](int c) 
#else
		for (int i = 0; i < animA->vertices.size(); i++) 
#endif
		{
			int i = subsamples[c];
			if (corrs[i].valid) 
			{
				if (enableEPoint) {
					computeJFitPoint(i);	
				}
				//curRow++;
				if (enableEPlane)
					computeJFitPlane(i);			
				//curRow++;
			}
		}
#ifdef ENABLE_PARALLEL
	);
#endif
		

		finalizeJacobian();
	}

	//void multiplySparseMatrices(const Eigen::SparseMatrix<Float>& lhs,const Eigen::SparseMatrix<Float>& rhs, Eigen::SparseMatrix<Float>& res, Float tolerance) {
	void multiplySparseMatrices(const Eigen::SparseMatrix<Float>& lhs,const Eigen::SparseMatrix<Float>& rhs, MatrixX& res, Float tolerance) {		

		int rows = lhs.innerSize();
		int cols = rhs.outerSize();		
		//res.resize(lhs.innerSize(), rhs.outerSize());
		int expectedNNZ = lhs.nonZeros() + rhs.nonZeros();
		int totalNNZ = 0;
		//double ratioColRes = double(numElements)/double(lhs.rows()*rhs.cols());
		double ratioColRes = 0.5;
		tbb::parallel_for(size_t(0), size_t(cols), size_t(1), [&](int j) 
		//for (int j = 0; j < cols; j++)
		{
			Eigen::internal::AmbiVector<Float,int> tempVector(rows);
			tempVector.init(ratioColRes);
			tempVector.setZero();
			for (Eigen::SparseMatrix<Float>::InnerIterator rhsIt(rhs, j); rhsIt; ++rhsIt)
			{
				// FIXME should be written like this: tmp += rhsIt.value() * lhs.col(rhsIt.index())
				tempVector.restart();
				Float x = rhsIt.value();
				for (Eigen::SparseMatrix<Float>::InnerIterator lhsIt(lhs, rhsIt.index()); lhsIt; ++lhsIt)
				{
					tempVector.coeffRef(lhsIt.index()) += lhsIt.value() * x;
				}
			}
			int count = 0;
			for (Eigen::internal::AmbiVector<Float,int>::Iterator it(tempVector,tolerance); it; ++it) {
				//res.insert(it.index(),j) = it.value();

				res(it.index(),j) = it.value();
				count++;
			}
			totalNNZ += count;			
		}
		);

	}

	void doGaussNewton() {				
		
		//Eigen::SparseMatrix<Float, Eigen::ColMajor, int> dampingMatrix(numFreeVariables, numFreeVariables);
		VectorX x0(numFreeVariables);		
		VectorX x(numFreeVariables);
		VectorX rhs(numFreeVariables);		


		for (int i = 0; i < numFreeVariables; i++) {
			x0[i] = 0.0;
		}

		//Eigen::SparseMatrix<Float, Eigen::ColMajor, int> tmpJ = jacobian;

		//lhs = (jacobian.transpose() * jacobian).pruned(1e-6) + dampingMatrix;


#ifdef USE_NORMAL_LSQR
		LOG(INFO) << "Multiplying matrices ... " << endl;
		multiplySparseMatrices(jacobianTranspose, jacobian, normalMatrix, 1e-16);

		LOG(INFO) << "Multiplying matrices ... DONE " << endl;

		normalMatrix += dampingMatrix;
#endif

		rhs = jacobianTranspose * energy;
		x0 = flatten();


		

		LOG(INFO) << "Solving matrices ... " << endl;
		//Eigen::PardisoLDLT< Eigen::SparseMatrix<Float>>  solver;		
		//solver.compute(normalMatrix);		
		//x = solver.solve(rhs);
				
		Eigen::ConjugateGradient<MatrixX> cg;
		cg.compute(normalMatrix);
		cg.setTolerance(1e-2);
		x = cg.solveWithGuess(rhs, x0);
		VLOG(1) << "#iterations:     " << cg.iterations() << std::endl;
		VLOG(1) << "estimated error: " << cg.error()      << std::endl;

		unflatten(x0 - x);

		deformAnimation(curAnim,A,b);
		//energyFile << "Energy after solve:" << endl;
		computeEnergy();				
		//energyFile << "---------------" << endl;

		LOG(INFO) << "Solving matrices... DONE" << endl;		
	}
	
	void clearJacobian() {
		tripletList.clear();
		//jacobian.setFromTriplets(tripletList);
		
	}

	void addToJacobian(int i, int j, Float val) {
		//tripletList.push_back(TripletF(i,j,val));				
		jacobian.insert(i,j) = val;
		jacobianTranspose.insert(j,i) = val;
	}

	void finalizeJacobian() {
		//jacobian.setFromTriplets(tripletList.begin(), tripletList.end());
		//jacobian.makeCompressed();
		//jacobianTranspose.makeCompressed();
		
	}

	void computeErrors() {
		Error e;
		curRow = 0;

		for (int i = 0; i < nodes.size(); i++) {
			//LOG(INFO) << "Node:" << i << endl;
			//LOG(INFO) << A[i] << endl;
			//LOG(INFO) << b[i] << endl;

			for (int j = 0; j < nodeWeights[i].size(); j++) {
				e.avgSmooth += energy[curRow] * energy[curRow];
				e.maxSmooth = max(e.maxSmooth, energy[curRow]);
				curRow++;
			}
		}
		

		Float invalidCorrError = 1.0;
		for (int i = 0; i < animA->vertices.size(); i++) {
			if (!corrs[i].valid) 
			{
				e.avgPoint += invalidCorrError * invalidCorrError;
				e.avgPlane += invalidCorrError * invalidCorrError;
			}
			else {
				e.avgPoint += (energy[curRow] * energy[curRow]);
				e.maxPoint = max(e.maxPoint, energy[curRow]);
				curRow++;
				e.avgPlane += (energy[curRow] * energy[curRow]);
				e.maxPlane = max(e.maxPlane, energy[curRow]);
				curRow++;
			}
		}

		e.avgSmooth = sqrt(e.avgSmooth) / nodes.size();
		e.avgPoint  = sqrt(e.avgPoint) / animA->vertices.size();
		e.avgPlane  = sqrt(e.avgPlane) / animA->vertices.size();

		errors.push_back(e);
	}

	void doIteration() {
		LOG(INFO) << iteration << " / " << maxIterations;
		LOG(INFO) << "Deforming animation";
		deformAnimation(curAnim, A, b);

		
		subsamples.clear();		
		subsample(subsamples);
		

		LOG(INFO) << "Computing correspondences";

		computeCorrespondences();		
		
		LOG(INFO) << "Setting up matrices";
		setupMatrices();
		

		LOG(INFO) << "Computing energy";
		computeEnergy();

		LOG(INFO) << "Computing Jacobian";
		computeJacobian();

		LOG(INFO) << "Computing Gauss Newton";
		doGaussNewton();



		//computeErrors();


		if (iteration % saveEvery == 0 || iteration == maxIterations) {
			deformAnimation(curAnim, A, b);
			computeCorrespondences(true);

			//writeVertexDataToDisk();
			writeSolution();

			saveRegistration();	// Unprojected

			//Smooth correspondences
			LOG(INFO) << "Smoothing ...";		
			smoothCorrespondences(3);
			LOG(INFO) << "Smoothing ... " << "DONE" ;

			for (int i = 0; i < curAnim->vertices.size(); i++) {
				if (corrs[i].valid)
					curAnim->vertices[i] = corrs[i].p;
			}
			saveRegistration(true);	// Projected

		}
		if (iteration > 0) {
			if ( abs(curNormEnergy - prevNormEnergy) < 0.01 * prevNormEnergy ) {			
				LOG(INFO) << "Lowering weights" << endl;
				curWeights.weightSmooth /= sqrt(2);
			}
		}

		iteration++;
	}

	~ICP() {
		cleanup();
	}

	void cleanup() {
		LOG(INFO) << "Cleaning up ICP";

		nodes.clear();
		A.clear();
		b.clear();		
		corrs.clear();
		jacobian.data().clear();
		jacobianTranspose.data().clear();
		
		for (int i = 0; i < nodeWeights.size(); i++) {
			for (int j = 0; j < nodeWeights[i].size(); j++) 
				delete nodeWeights[i][j];

			nodeWeights[i].clear();
		}
		nodeWeights.clear();

		for (int i = 0; i < vertexWeights.size(); i++) {
			vertexWeights[i].clear();
		}
		vertexWeights.clear();
	}

	// Find correspondences for each vertex and creates a map: corr[indexA] -> Correspondence<DIM> 
	void computeCorrespondences(bool useFullMesh=false) {
		
		vector<int> *vertexList;
		if (useFullMesh) {
			vertexList = &allVertexList;
		}
		else
		{
			vertexList = &subsamples;
		}

		tbb::parallel_for(size_t(0), vertexList->size(), size_t(1), [&](int i) 	 
		{
			corrs[i].valid = false;
		});

#ifdef ENABLE_PARALLEL
		tbb::parallel_for(size_t(0), vertexList->size(), size_t(1), [&](int c) 		
#else
		for (int i = 0; i < animA->vertices.size(); i++) 
#endif
		{
			int i = vertexList->at(c);
			Correspondence<DIM> corr;
				
			handleUserCorr(i, corr);

		
			
			if (!corr.valid) {

				tbb::concurrent_hash_map<int, bool>::const_accessor tmpBool;
				if (animA->topoChangeVertices.find(tmpBool,i)  || 
					(tinyComponentVertices.find(tmpBool,i) )) {
					corr.valid = false;
					corrs[i] = corr;
				}
				
				else {					
					corr = animB->findClosestPoint(curAnim->vertices[i]);
					VectorN curNormal = curAnim->normals[i];
					Float d = (corr.p - curAnim->vertices[i]).norm();
					
					
					FrameIndex fIndex = curAnim->globalToFrameMap.at(i);



					// Try basic test
					if (d < maxCorrDist 
						&& 
						curNormal.dot(corr.n) > minCorrAngle
						) {
						corrs[i] = corr;
					}
					else {			
						auto fIndex = animA->globalToFrameMap.find(i);
						Float stepSize = abs(min(animA->timeScaling, animB->timeScaling));
						Float curStep = stepSize;
						bool found = false;

						
						while (curStep < maxCorrDist) {
							VectorN plusPoint = curAnim->vertices[i] + curStep * curNormal;
							VectorN minusPoint = curAnim->vertices[i] - curStep * curNormal;

							corr = animB->findClosestPoint(plusPoint);

							d = (corr.p - curAnim->vertices[i]).norm();							
							if (d < maxCorrDist && curNormal.dot(corr.n) > minCorrAngle) 
							{
								corrs[i] = corr;
								found = true;								
								break;
							}
							
							corr = animB->findClosestPoint(minusPoint);
						
							d = (corr.p - curAnim->vertices[i]).norm();							
							if (d < maxCorrDist && curNormal.dot(corr.n) > minCorrAngle) {
								corrs[i] = corr;
								found = true;																
								break;
							}
						

							curStep += stepSize;
						}

						if (!found)
						{							
							corrs[i] = Correspondence<DIM> ();
						}
					}
				}

				if (animA->isVertexOnStartFrame(i)) {
					corrs[i].p[DIM] = animB->getDepth(0);
				}
				if (animA->isVertexOnEndFrame(i)) {
					corrs[i].p[DIM] = animB->getDepth(animB->frames.size()-1);
				}
			}

		}
#ifdef ENABLE_PARALLEL
		);
#endif


	}


	void smoothCorrespondences(int numPasses) {

		tmpCorrs = corrs;


		for (int k = 0; k < numPasses; k++) {
			tbb::parallel_for(size_t(0), animA->vertices.size(), size_t(1), [&](int i)
			{
				//
				VectorN sumVector;

				for (int j = 0; j < DIM+1; j++)
					sumVector[j] = 0.0;

				int numNeighbours = 0;

				if (corrs[i].valid) {
					sumVector = corrs[i].p - animA->vertices[i];
					numNeighbours++;
				}
				

				
				bool firstElement = true;
				for (int j = 0; j < animA->neighbours[i].size(); j++) {
					GlobalIndex nbIndex = animA->neighbours[i][j].index;
					if (nbIndex != Globals::InvalidNeighbour) {
						if (corrs[nbIndex].valid) {
							numNeighbours++;
							sumVector += corrs[nbIndex].p - animA->vertices[nbIndex];						
						}
					}
				}

				if (numNeighbours > 0)
				{
					tmpCorrs[i].p = (sumVector / numNeighbours) + animA->vertices[i];
					tmpCorrs[i].valid = corrs[i].valid;

					//Correspondence<DIM> corr = animB->findClosestPoint(tmpCorrs[i].p);
					//if (corr.valid)
					//	tmpCorrs[i].p = corr.p;

				}

			});			
			corrs = tmpCorrs;
		}
	}

	// Handles correspondence for vertex i 
	void handleUserCorr(int i, Correspondence<DIM>& corr) {
		if (isPointCorrespondence(i)) {				
			corr.valid = true;
			corr.p = animB->vertices[ pointCorrs[i] ];
			corr.n = animB->normals[ pointCorrs[i] ];
			corrs[i] = corr;

			if (enforceHardConstraints && isNode(i)) {
				int nodeIndex = vertexToNodeIndexMap[i];
				//b[nodeIndex] = corr.p - animA->vertices[i];
			}
		}
		else if (isSpaceCorrespondence(i)) {						
				corr.valid = true;
				corr.p = animA->vertices[i] + spaceCorrs[i];				
				corr.n = animB->normals[ animB->findClosestVertex(corr.p) ];
				corrs[i] = corr;			
		}
		else if (animA->isLineFeature(i) && isLineCorrespondence(animA->getLineFeature(i).first)) {

			Float closestDist = 1e3;
			int closestVertexIndex = -1;



			pair<int,int> lineData = animA->getLineFeature(i);
			int lineIndexA = lineData.first;
			int indexAlongLineA = lineData.second;

			int vIndexA =  animA->lineFeatures[lineIndexA][0];
			auto it = animA->globalToFrameMap.find(vIndexA);
			int startA = it->second.first;
			int endA = startA + animA->lineFeatures[lineIndexA].size() - 1;

			for (int j = 0; j < lineCorrs[lineIndexA].size(); j++) {
				LineCorr lc = lineCorrs[lineIndexA][j];
				int lineIndexB = lc.first;
				int offset = lc.second;

				int vIndexB = animB->lineFeatures[lineIndexB][0];
				auto it = animB->globalToFrameMap.find(vIndexB);
				int startB = it->second.first;
				int endB = startB + animB->lineFeatures[lineIndexB].size() - 1;

				int time = startA + indexAlongLineA + offset;
				if (time >= startB && time <= endB) {
					corr.valid = true;			
					closestVertexIndex = animB->lineFeatures[lineIndexB][time - startB];

				}

				//VLOG(1) <<"Vertex: " << i << " line : " << lineIndexA << " " << " offset : " << offset << " time: " << time <<  " startB: " << startB << " endB:" << endB << endl;

			}

			// Uncomment this if you want to resort to closest point matching for lines
			//vector<int> pointsB = linePointCorrs[ animA->getLineFeature(i).first ];
			//for (int j = 0; j < pointsB.size(); j++) {						
			//	int indexB = pointsB[j];
			//	Float d = (animB->vertices[indexB] - curAnim->vertices[i]).norm();
			//	if (d < closestDist) {
			//		closestVertexIndex = indexB;
			//		closestDist = d;
			//	}
			//}

			if (closestVertexIndex != -1)
			{
				corr.p = animB->vertices[closestVertexIndex];
				corr.n = animB->normals[closestVertexIndex];
				corrs[i] = corr;
			}
		}
	}

	// Do this for every iteration
	void setupMatrices() {
		jSmoothMap.clear();
		jFitPlaneMap.clear();
		jFitPointMap.clear();
		numRows = 0;		

		int termCount = 0;

		if (enableESmooth) {
			// JSmooth
			for (int i = 0; i < nodes.size(); i++) {			
				numRows += (DIM+1)*nodeWeights[i].size();
				for (int j = 0; j < nodeWeights[i].size(); j++) {
					jSmoothMap[ pair<int, int> (i, j) ] = termCount;
					termCount+=(DIM+1);
				}
			}
		}

		// JPoint, JPlane
		for (int c = 0; c < subsamples.size(); c++) {

			int i = subsamples[c];

			if (!corrs[i].valid) 
				continue;
			else {			
				if (enableEPoint) {
					jFitPointMap[i] = termCount;
					termCount+=(DIM+1);
					numRows+=(DIM+1);
				}
				if (enableEPlane) {
					jFitPlaneMap[i] = termCount;
					termCount++;
					numRows++;
				}
			}
			
		}

		//jacobian = Eigen::SparseMatrix<Float, Eigen::RowMajor, int>(numRows, numCols);		
		//jacobianTranspose = Eigen::SparseMatrix<Float, Eigen::ColMajor, int>(numCols, numRows);		

		jacobian = Eigen::SparseMatrix<Float, Eigen::RowMajor, int>(numRows, numCols);		
		jacobianTranspose = Eigen::SparseMatrix<Float, Eigen::ColMajor, int>(numCols, numRows);		

		VLOG(1) << "# non-zeros:" << jacobian.nonZeros() << endl;
		//jacobian.makeCompressed();
		//jacobianTranspose.makeCompressed();

		VLOG(1) << "Max neighbours:" << maxNumNeighbours << endl;
		VLOG(1) << "# rows:" << numRows << endl;
		VLOG(1) << "Expected:" << numRows * maxNumNeighbours * ((DIM+1)*(DIM+1)+(DIM+1)) << endl;
		jacobian.reserve(Eigen::VectorXi::Constant(numRows, maxNumNeighbours * ((DIM+1)*(DIM+1)+(DIM+1))) );
		jacobianTranspose.reserve(Eigen::VectorXi::Constant(numRows,maxNumNeighbours * ((DIM+1)*(DIM+1)+(DIM+1))));

		
		
		int numElements = maxNumNeighbours * ((DIM+1)*(DIM+1)+(DIM+1)) ;

		static bool firstRun = true;
		
		if (firstRun) {
			normalMatrix.resize(numFreeVariables, numFreeVariables);
			dampingMatrix.resize(numFreeVariables, numFreeVariables);		
			firstRun = false;

			Float damping = 1e-6;
			for (int i = 0; i < numFreeVariables; i++) {			
				dampingMatrix(i,i) = damping;			
				//dampingMatrix.insert(i,i) = damping;
			}

		}

		//normalMatrix.reserve(Eigen::VectorXi::Constant(numFreeVariables, numFreeVariables));

		energy = VectorX(numRows);
		curRow = 0;

	}

	// Deforms animA using A,b and stores the mesh in dst.
	void deformAnimation(const unique_ptr<Animation<DIM>>& dst, const vector<MatrixN>& A, const vector<VectorN>& b ) 
	{		

		// Now apply deformation to vertices of src and store them in dst
		//for (int i = 0; i < animA->vertices.size(); i++) 
		tbb::parallel_for(size_t(0), size_t(animA->vertices.size()), size_t(1), [&](int i) 	
		{
			dst->vertices[i] = VectorN::Zero();
			for (int j = 0; j < vertexWeights[i].size(); j++) {
				int nodeIndex = vertexWeights[i][j].nodeIndex;
				Float nodeWeight = vertexWeights[i][j].weight;
				VectorN nodePos = animA->vertices[nodes[nodeIndex]];
				
				FrameIndex fIndex = animA->globalToFrameMap.at(i);

				
				if (animA->frames[fIndex.first]->vertices[fIndex.second].isTinyComponent) {
					
				}
				else  
					dst->vertices[i] += nodeWeight * (A[nodeIndex] * (animA->vertices[i] - nodePos) + nodePos + b[nodeIndex]);
			}
		});

		computeDropletCentroids(dst);
	}

	void computeDropletCentroids(const unique_ptr<Animation<DIM>>& dst) {
		for (int i = 0; i < animA->frames.size(); i++) {
			for (int j = 0; j < animA->frames[i]->droplets.size(); j++) {
				Droplet<DIM> droplet = animA->frames[i]->droplets[j];

				// Compute the centroid in space-time coordinates
				VectorN centroid;
				for (int k = 0; k < DIM; k++)
					centroid[k] = droplet.centroid[k];
				centroid[DIM] = animA->getDepth(i);
				
				int closestVertexIndex = animA->findClosestVertex(centroid);

				// Use the weights from this point to deform the centroid

				VectorN newCentroid = VectorN::Zero();

				for (int k = 0; k < vertexWeights[closestVertexIndex].size(); k++) {
					int nodeIndex = vertexWeights[closestVertexIndex][k].nodeIndex;
					Float nodeWeight = vertexWeights[closestVertexIndex][k].weight;
					VectorN nodePos = animA->vertices[nodes[nodeIndex]];

					newCentroid += nodeWeight * (A[nodeIndex] * (centroid - nodePos) + nodePos + b[nodeIndex]);
				}

				// Now move all the points using the centroid

				for (int k = 0; k < droplet.vertices.size(); k++) {
					GlobalIndex gIndex = animA->frameToGlobalMap.at(FrameIndex(i, droplet.vertices[k]));
					dst->vertices[gIndex] = newCentroid + animA->vertices[gIndex] - centroid;
				}

			}
		}
	}

	// Returns transform for deformation node i.
	void getNode(int i, MatrixN& a, VectorN& b) {
		a = A[i];
		b = b[i];
	}

	// Set transform for a specific deformation node
	void setNode(int i, const MatrixN& a, const VectorN& b) {
		A[i] = a;
		b[i] = b;
	}

	// Get transforms for all deformation nodes
	void getNodes(vector<MatrixN>& a, vector<VectorN>& b) {
		a = A;
		b = b;
	}

	// Set transforms for all deformation nodes
	void setNodes(const vector<MatrixN>& a,const vector<VectorN>& b) {
		A = a;
		b = b;
	}

	size_t getNumNodes() {
		return nodes.size();
	}
	
	// Unflattens into A and b
	void unflatten(const VectorX& x) {
		int count = 0;
		for (int i = 0; i < nodes.size(); i++) {			

			for (int col = 0; col < DIM+1; col++) {
				for (int row = 0; row < DIM+1; row++) {
					A[i](row,col) = x[count];
					count++;
				}
			}

			int numRows = DIM+1;
			int gIndexNode = nodes[i];
			if (isPointCorrespondence(gIndexNode))  {					
				if (enforceHardConstraints) {					
				} 
			}
			else if (isSpaceCorrespondence(gIndexNode)) {				
				b[i][DIM] = x[count];
				count++;
			}
			else {					
				if (enforceEqualLength && (isNodeOnStartFrame(i) || isNodeOnEndFrame(i))) 
					numRows = DIM;

				for (int row = 0; row < numRows; row++) {
					b[i][row] = x[count];
					count++;
				}	
			}
		
		}		
	}

	// Flattens A, b into vector of free variables
	VectorX flatten() {
		VectorX x(numFreeVariables);
		int count = 0;
		for (int i = 0; i < nodes.size(); i++) {			
				
			for (int col = 0; col < DIM+1; col++) {
				for (int row = 0; row < DIM+1; row++) {
					x[count] = A[i](row,col);
					count++;
				}
			}

			int numRows = DIM+1;
			int gIndexNode = nodes[i];
			if (isPointCorrespondence(gIndexNode))  {					
				if (enforceHardConstraints) {					
				} 
			}
			else if (isSpaceCorrespondence(gIndexNode)) {				
				x[count] = b[i][DIM]; 
				count++;
			}
			else {					
				if (enforceEqualLength && (isNodeOnStartFrame(i) || isNodeOnEndFrame(i))) 
					numRows = DIM;

				for (int row = 0; row < numRows; row++) {
					x[count] = b[i][row]; 
					count++;
				}		
			}
	
		}				
		return x;
	}

	// Finite difference version of JSmooth
	VectorX computeJSmoothFD(int nodeIndexI, int j) {
		saveState();
		Float eps = 1e-6;
		VectorX jacobian(numFreeVariables);
		for (int i = 0 ; i < numFreeVariables; i++) 
			jacobian[i] = 0.0;

		VectorX initialGuess = flatten();
		Float initialEnergy = computeESmooth(nodeIndexI, j);

		for (int i = 0; i < numFreeVariables; i++) {
			VectorX newGuess = initialGuess;
			newGuess[i] += eps;

			unflatten(newGuess);
			Float newEnergy = computeESmooth(nodeIndexI,j);			
			jacobian[i] = (newEnergy - initialEnergy) / eps;
		}
		restoreState();

		return jacobian;
	}

	// Finite difference version of JFitPlane
	VectorX computeJFitPlaneFD(GlobalIndex globalIndex) {
		saveState();
		Float eps = 1e-9;
		VectorX jacobian(numFreeVariables);
		for (int i = 0 ; i < numFreeVariables; i++) 
			jacobian[i] = 0.0; 

		VectorX initialGuess = flatten();
		deformAnimation(curAnim, A, b);
		Float initialEnergy = computeEFitPlane(globalIndex);

		for (int i = 0; i < numFreeVariables; i++) {
			VectorX newGuess = initialGuess;
			newGuess[i] += eps;

			unflatten(newGuess);			
			deformAnimation(curAnim, A, b);

			Float newEnergy = computeEFitPlane(globalIndex);		
			jacobian[i] = (newEnergy - initialEnergy) / eps;
		}

		restoreState();
		return jacobian;
	}

	// Finite difference version of JFitPoint
	VectorX computeJFitPointFD(GlobalIndex globalIndex) {
		saveState();
		Float eps = 1e-9;
		VectorX jacobian(numFreeVariables);
		for (int i = 0 ; i < numFreeVariables; i++) 
			jacobian[i] = 0.0; 

		VectorX initialGuess = flatten();
		deformAnimation(curAnim, A, b);
		Float initialEnergy = computeEFitPoint(globalIndex);

		for (int i = 0; i < numFreeVariables; i++) {
			VectorX newGuess = initialGuess;
			newGuess[i] += eps;

			unflatten(newGuess);			
			deformAnimation(curAnim, A, b);

			Float newEnergy = computeEFitPoint(globalIndex);		
			jacobian[i] = (newEnergy - initialEnergy) / eps;
		}

		restoreState();
		return jacobian;
	}

	// Saves state of free variables
	void saveState() {
		tmpA = A;
		tmpB = b;
	}

	// Restores state of free variables
	void restoreState() {
		A = tmpA;
		b = tmpB;
	}

	const vector<Error>& getErrors() {
		return errors;
	}
	
	const unique_ptr<Animation<DIM>>& getCorrAnim() {
		return curAnim;
	}

	// Writes all vertex data to disk
	bool writeVertexDataToDisk() {
		LOG(INFO) << "Writing vertex data to disk ... " << animA->vertices.size() ;
		boost::filesystem::path cfgDirPath("vdata/"+cfgName);
		boost::filesystem::create_directory(cfgDirPath);

		//for (int i = 0; i < animA->numFrames; i++)
		tbb::parallel_for(size_t(0), size_t(animA->numFrames), size_t(1), [&](int i) 		
		{	
			char iterDir[1024];
			sprintf_s(iterDir, 1024, "vdata/%s/iter%02d",cfgName.c_str(), iteration);
			boost::filesystem::path iterDirPath(iterDir);
			boost::filesystem::create_directory(iterDirPath);

			char fullPath[1024];
			sprintf_s(fullPath, 1024, "vdata/%s/iter%02d/frame_%04d.vdata", cfgName.c_str(), iteration, i);

			gzFile gzf;
			gzf = gzopen(fullPath, "wb");


			// Write # vertices
			gzwrite(gzf, &animA->frames[i].numVertices, sizeof(animA->frames[i]->numVertices));

			for (int j = 0; j < animA->frames[i]->numVertices ; j++) {
				int gIndex = animA->frameToGlobalMap.at(FrameIndex(i,j));

				Correspondence<DIM> c = corrs[gIndex];
				gzwrite(gzf, &j, sizeof(j)); // spaceIndex
				gzwrite(gzf, &c.valid, sizeof(c.valid)); // Is corr valid

				// Corr
				for (int k = 0; k < DIM+1; k++)
					gzwrite(gzf, &c.p[k], sizeof(Float));

				// EPoint
				Float normalizedEPoint =  computeEFitPointVector(gIndex).norm() / maxEPointNorm;
				gzwrite(gzf, &normalizedEPoint, sizeof(Float));

				Float normalizeEPlane = abs(computeEFitPlane(gIndex)) / maxEPlaneNorm;

				bool isVertexNode = isNode(gIndex);
				// Is this vertex a node
				gzwrite(gzf, &isVertexNode, sizeof(isVertexNode));

				if (isVertexNode) {
					int nodeIndex = vertexToNodeIndexMap[gIndex];					
					const VectorN& curB = b[nodeIndex];
					for (int k = 0; k < DIM+1; k++) {				
						gzwrite(gzf, &curB[k], sizeof(Float));
					}
				}		
			}		
			gzclose(gzf);
		}
		);

		LOG(INFO) << "Done!" << endl;
		return true;
	}

	// Writes solution to disk (deformation nodes)
	bool writeSolution() {
		LOG(INFO) << "Writing solution to disk ... ";
		char fullPath[1024];
		sprintf_s(fullPath, 1024, "sol/%s_%04d.sol", cfgName.c_str(), iteration);

		try {
			gzFile gzf;
			gzf = gzopen(fullPath, "wb");
			//fstream f(fullPath, ios::out);

			//f << nodes.size() << endl;
			int numNodes = nodes.size();
			gzwrite(gzf, &numNodes, sizeof(numNodes));

			for (int i = 0; i < nodes.size(); i++) {

				// Global index
				FrameIndex fIndex = animA->globalToFrameMap.at(nodes[i]);
				//f << fIndex.first << " " << fIndex.second << " ";
				gzwrite(gzf, &fIndex.first, sizeof(fIndex.first));
				gzwrite(gzf, &fIndex.second, sizeof(fIndex.second));

				// A
				for (int r = 0; r < DIM+1; r++) {
					for (int c = 0; c < DIM+1; c++) {
						//f << A[i](r,c) << " ";
						gzwrite(gzf, &A[i](r,c), sizeof(Float));
					}
				} 
				//f << endl;

				// b
				for (int c = 0; c < DIM+1; c++) {
					//f << b[i][c] << " ";
					gzwrite(gzf, &b[i][c], sizeof(Float));
				}
				//f << endl;
			}

			//f.close();
			gzclose(gzf);
		}
		catch (exception e) {
			throw e;
			return false;
		}

		return true;
		LOG(INFO) << "DONE";
	}
	
	// Read solution from disk and create deformation nodes
	bool readSolution(string filename) {
		try {
			int numNodes;
			//fstream f(filename, ios::in);
			gzFile gzf;
			gzf = gzopen(filename.c_str(), "rb1");
			if (!gzf)
				LOG(INFO)<<"cannot open file " << filename;

			// read vertices

			if (gzread(gzf, &numNodes, sizeof(numNodes)) != sizeof(numNodes))
				LOG(INFO) << "Error reading file" << filename;
				
			 //f >> numNodes;

			 for (int i = 0; i < numNodes; i++) {
				 int frame;
				 int spaceIndex;
				 //f >> frame >> spaceIndex;

				 if (gzread(gzf, &frame, sizeof(frame)) != sizeof(frame))
					 LOG(INFO) << "Error reading file" << filename;

				 if (gzread(gzf, &spaceIndex, sizeof(spaceIndex)) != sizeof(spaceIndex))
					 LOG(INFO) << "Error reading file" << filename;

				 GlobalIndex gIndex = animA->frameToGlobalMap.at(FrameIndex(frame,spaceIndex));
				 MatrixN tmpA;
				 VectorN tmpB;

				 for (int r = 0; r < DIM+1; r++) {
					 for (int c = 0; c < DIM+1; c++) {
						 if (gzread(gzf, &tmpA(r,c), sizeof(Float)) != sizeof(Float))
							 LOG(INFO) << "Error reading file" << filename;

						 //f >> tmpA(r,c);
					 }
				 } 
				 
				 for (int c = 0; c < DIM+1; c++) {
					 //f >> tmpB[c];
					 if (gzread(gzf, &tmpB[c], sizeof(Float)) != sizeof(Float))
						 LOG(INFO) << "Error reading file" << filename;
				 }
				 
				 // Create a node here if one doesn't already exist
				 if (!isNode(gIndex)) {
					 createNode(gIndex);
				 }

				 // Set transform
				 int nodeIndex = vertexToNodeIndexMap[gIndex];
				 A[nodeIndex] = tmpA;
				 b[nodeIndex] = tmpB;
			 }

			 gzclose(gzf);

		}
		catch (exception e) {
			throw e;
			return false;
		}

		return true;
	}

	void saveRegistration(bool saveProjected=false)  {		

		LOG(INFO) << "Writing corrs to disk ... ";
		boost::filesystem::path p = boost::filesystem::canonical(boost::filesystem::path(configPath));
		string tmpPath = p.string();
		char filename[1024];


		// Write out original corr file
		if (saveProjected) {
			sprintf_s(filename, "registrations/%s-projected_%04d.corr",cfgName.c_str(), iteration);	
		}
		else {
			sprintf_s(filename, 1024, "registrations/%s_%04d.corr",cfgName.c_str(), iteration);	
		}
		
		writeCorrFile(filename, tmpPath, saveProjected);	

		LOG(INFO) << "Writing corrs to disk ... " << "DONE!";

	}

protected:
	bool isNodeOnStartFrame(int nodeIndex) {
		return animA->isVertexOnStartFrame(nodes[nodeIndex]);
	}

	bool isNodeOnEndFrame(int nodeIndex) {
		return animA->isVertexOnEndFrame(nodes[nodeIndex]);
	}

	bool isNode(int index) {
		return (vertexToNodeIndexMap.find(index) != vertexToNodeIndexMap.end());
	}
	
	bool isNodeUserCorrespondence(int nodeIndex);

	// Is point with index globalIndex a user defined correspondence?
	bool isPointCorrespondence(int globalIndex) {
		return (pointCorrs.find(globalIndex) != pointCorrs.end());
	}

	bool isSpaceCorrespondence(int globalIndex) {
		return (spaceCorrs.find(globalIndex) != spaceCorrs.end());
	}
	// Has this line index in animA been mapped to an index in B?
	bool isLineCorrespondence(int lineIndex) {
		if (lineCorrs.find(lineIndex) != lineCorrs.end())
			return true;
		return false;
	}


	void writeCorrFile(string filename, string path, bool writeProjected=false);

	vector< GlobalIndex > nodes;
	vector <MatrixN> A;
	vector <VectorN> b;
	vector <MatrixN> tmpA;
	vector <VectorN> tmpB;

	vector< Correspondence<DIM> > corrs;
	vector< Correspondence<DIM> > tmpCorrs;
	vector < vector<NeighbourNode> >  vertexWeights; // For all vertices
	vector < vector<LaplaceNeighbour*> > nodeWeights; // For all nodes (smoothness)
	map< GlobalIndex , int > vertexToNodeIndexMap;
	tbb::concurrent_hash_map<int, bool> tinyComponentVertices;
	
	
	Float prevNormEnergy, curNormEnergy;
	VectorX energy;
	Eigen::SparseMatrix<Float, Eigen::RowMajor, int> jacobian;
	Eigen::SparseMatrix<Float, Eigen::ColMajor, int> jacobianTranspose;	
	MatrixX normalMatrix;

	int numRows, numCols; // Rows and columns in the Jacobian
	int curRow;

	bool enableESmooth, enableEPoint, enableEPlane;

	vector < Error > errors;

	const unique_ptr<Animation<DIM>>& animA;
	const unique_ptr<Animation<DIM>>& animB;

	const int samplingRate;
	ICPWeights curWeights;
	const ICPWeights initialWeights;
	bool enforceEqualLength;
	bool enforceHardConstraints;
	bool usePreviousSolution;
	string prevSolPath;
	const int maxIterations;
	const int saveEvery;
	string configPath;

	unique_ptr<Animation<DIM>>  curAnim; // Deformed version of anim A
	Float nodeRadius;
	Float vertexRadius;
	int maxNumNeighbours;
	map < int , vector<int> > linePointCorrs;
	map < int , vector<LineCorr> > lineCorrs; // 1 - many 
	map < int, VectorN > spaceCorrs;
	map < int, int > pointCorrs; // 1-1 
	Float minCorrAngle;
	Float maxCorrDist;
	vector<int> nodeOffsets;
	int numFreeVariables;
	map< pair<int, int> , int > jSmoothMap;
	map< int, int > jFitPlaneMap;
	map< int, int > jFitPointMap;
	
	int iteration;
	MatrixX dampingMatrix;

	vector < vector<NeighbourNode> > nodeToVertexWeights;
	bool dumpDeformationNodesToDisk;
	fstream energyFile;
	std::string cfgName;

	vector<int> subsamples;
	int subsamplingRate;
	Float maxEPointNorm, maxEPlaneNorm;
	
	bool generateProjectedOutput;

	vector<int> allVertexList;
	vector<VertexError> perVertexErrorList;

	FRIEND_TEST(ICPTest, Constructor);
	FRIEND_TEST(ICPTest, DISABLED_Correspondences);
	FRIEND_TEST(ICPTest, DISABLED_ESmooth);
	FRIEND_TEST(ICPTest, DISABLED_JSmooth);
	FRIEND_TEST(ICPTest, JFitPlane);
	FRIEND_TEST(ICPTest, JFitPoint);

	FRIEND_TEST(ICPTest, Drop);
};

typedef ICP<2> ICP2D;