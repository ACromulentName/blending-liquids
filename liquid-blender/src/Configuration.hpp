#pragma once

#include "stdafx.h"
#include "Mesh2D.h"
#include "Animation.h"
#include <queue>
#include <set>
#include <fstream>
#include "boost/filesystem.hpp"
#include "tbb/atomic.h"
#include <functional>
#include "json.h"
using namespace std;



typedef map < int, int > PointCorrs; // 
typedef pair<int, int> LineCorr; // <lineIndex, offset>
typedef map < int, vector<LineCorr> > LineCorrs; //1-many


struct ICPWeights
{
	Float weightRigid;
	Float weightPlane;
	Float weightPoint;
	Float weightSmooth;

	ICPWeights() {
		weightRigid = 0;
		weightSmooth = sqrt(1.0);
		weightPoint = sqrt(0.1);
		weightPlane = sqrt(1);
	}
};


template<int DIM>
struct Configuration {
	typedef Eigen::Matrix<Float, DIM+1, 1> VectorN;
	typedef map < int, VectorN > SpaceCorrs;

	Configuration() {
		enforceEqualLength = true;
		enforceHardConstraints = false;
	}

	Configuration(string configFilename, bool onlyLoadAnimA = false) : Configuration() {		
		configPath = configFilename;
		load(configFilename, onlyLoadAnimA);
	}

	bool load(string configFilename, bool blend = false) {
		enforceEqualLength = true;

		configPath = configFilename;

		int start = 0;
		int end = 0;
		for (int i = configFilename.length() - 1; i >=0; i--) {
			if (configFilename[i] == '.')
				end = i;
			if (configFilename[i] == '/' || configFilename[i] == '\\') {
				start = i+1;
				break;
			}
		}

		configName = configFilename.substr(start, end-start);

		std::string fileContents;
		Json::Value root;   // will contains the root value after parsing.
		Json::Reader reader;	

		bool readSuccess = readFileIntoString(configFilename, fileContents);
		if (!readSuccess) {
			LOG(ERROR) << "Failed to read config file into string: "<< configFilename;
			return false;
		}

		bool parsingSuccessful = reader.parse(fileContents, root );
		if ( !parsingSuccessful )
		{
			// report to the user the failure and their locations in the document.
			LOG(ERROR) << "Failed to parse configuration\n" 
				<< reader.getFormattedErrorMessages();
			return false;
		}

		if (!root["enforceEqualLength"].isNull())
			enforceEqualLength = root["enforceEqualLength"].asBool();


		if (root["animA"].isNull()) {
			LOG(ERROR) << "No entry for animation A!";
			return false;
		}
		else
		{
			std::string path = root["animA"]["path"].asString();
			int startFrame = root["animA"]["startFrame"].asInt();
			int endFrame = root["animA"]["endFrame"].asInt();

			VectorN translateVec;
			if (!root["animA"]["translate"].isNull()) {				
				for (int i = 0; i < root["animA"]["translate"].size(); i++) {
					translateVec[i] = root["animA"]["translate"][i].asFloat();					
				}

				LOG(INFO) << "Anim A translate vector:" << translateVec;
			}

			Float rotateX = 0.0;
			if (!root["animA"]["rotateX"].isNull()) {
				rotateX = root["animA"]["rotateX"].asFloat();
				LOG(INFO) << "Anim A rotation:" << rotateX;
			}

			animA->setRotationX(rotateX);
			animA->setTranslation(translateVec);

			bool success = animA->load(path, startFrame, endFrame, blend);
			if (!success)
				return false;



		}

		if (blend)
			return true;

		if (root["animB"].isNull()) {
			LOG(ERROR) << "No entry for animation B!";
			return false;
		}
		else
		{
			std::string path = root["animB"]["path"].asString();
			int startFrame = root["animB"]["startFrame"].asInt();
			int endFrame = root["animB"]["endFrame"].asInt();

			VectorN translateVec;
			if (!root["animB"]["translate"].isNull()) {				
				for (int i = 0; i < root["animB"]["translate"].size(); i++) {
					translateVec[i] = root["animB"]["translate"][i].asFloat();					
				}

				LOG(INFO) << "Anim B translate vector:" << translateVec;
			}

			Float rotateX  = 0.0;
			if (!root["animB"]["rotateX"].isNull()) {
				rotateX = root["animB"]["rotateX"].asFloat();
				LOG(INFO) << "Anim B rotation:" << rotateX;
			}

			animB->setRotationX(rotateX);
			animB->setTranslation(translateVec);
			animB->load(path, startFrame, endFrame, blend);

			
		}


		generateProjectedOutput = root.get("generateProjectedOutput", false).asBool();
		maxIterations = root.get("numIterations", 10).asInt();
		enforceHardConstraints = root.get("enforceHardConstraints", false).asBool();

		LOG(INFO) << "Enforce hard constraints: " << enforceHardConstraints;
		if (root["weights"].isNull()) {
			LOG(ERROR) << "No weights specified!";
			return false;
		}
		else {
			weights.weightRigid = root["weights"]["rigid"].asFloat();
			weights.weightSmooth = root["weights"]["smooth"].asFloat();
			weights.weightPoint = root["weights"]["point"].asFloat();
			weights.weightPlane = root["weights"]["plane"].asFloat();

			weights.weightRigid = sqrt(weights.weightRigid);
			weights.weightSmooth = sqrt(weights.weightSmooth);
			weights.weightPoint = sqrt(weights.weightPoint);
			weights.weightPlane = sqrt(weights.weightPlane);
		}

		if (root["samplingRate"].isNull()) {
			LOG(ERROR) << "Sampling rate is missing!";
			return false;
		}
		else {
			samplingRate = root["samplingRate"].asInt();
		}

		if (root["subsamplingRate"].isNull()) {
			LOG(ERROR) << "Subsampling rate is not specified!";
			subsamplingRate = 1;
		}
		else {
			subsamplingRate = root["subsamplingRate"].asInt();
		}

		if (root["saveEvery"].isNull()) {
			LOG(ERROR) << "saveEvery is missing!";
			saveEvery = maxIterations;
		}
		else {
			saveEvery = root["saveEvery"].asInt();
		}

		if (root["usePreviousSolution"].isNull()) {
			usePreviousSolution = false;
		}
		else {
			usePreviousSolution = root["usePreviousSolution"]["enabled"].asBool();
			if (usePreviousSolution) {
				prevSolPath = root["usePreviousSolution"]["path"].asString();
				LOG(INFO) << "Using previous solution : " << prevSolPath;
			}
			
		}

		if (root["correspondences"].isNull()) {
			LOG(ERROR) << "Correspondence file is missing!";
			return false;
		}
		else {
			std::string corrFile = root["correspondences"]["path"].asString();
			bool isEnabled = root["correspondences"]["enabled"].asBool();

			bool reverseCorrs = false;
			if (!root["correspondences"]["reverse"].isNull()) {
				reverseCorrs = root["correspondences"]["reverse"].asBool();
			}


			LOG(INFO) << "Ucorrs enabled:" << isEnabled;

			if (isEnabled)
				loadUserCorrespondences<DIM>(corrFile, animA, animB, pointCorrs, lineCorrs, spaceCorrs, samplingRate, reverseCorrs);
		}


		return true;
	}

	unique_ptr<Animation<DIM>> animA;
	unique_ptr<Animation<DIM>> animB;

	int maxIterations;
	int samplingRate;
	int subsamplingRate;
	ICPWeights weights;
	string configName;
	PointCorrs pointCorrs;
	LineCorrs lineCorrs;
	SpaceCorrs spaceCorrs;
	bool enforceEqualLength;
	bool enforceHardConstraints;
	string configPath;
	int saveEvery;

	bool usePreviousSolution;
	bool generateProjectedOutput;
	string prevSolPath;

	private:
		Configuration(const Configuration<DIM>&) {}
		Configuration& operator=(const Configuration<DIM>&) { return *this;}
};

template<>
Configuration<2>::Configuration() {
	animA = unique_ptr<Animation2D>(new Animation2D());
	animB = unique_ptr<Animation2D>(new Animation2D());
}

template<>
Configuration<3>::Configuration() {
	animA = unique_ptr<Animation3D>(new Animation3D());
	animB = unique_ptr<Animation3D>(new Animation3D());
}

template<int DIM>
bool loadUserCorrespondences(string filename, const unique_ptr<Animation<DIM>> &animA, const unique_ptr<Animation<DIM>> & animB, PointCorrs& pointCorrs, LineCorrs& lineCorrs, map < int, Eigen::Matrix<Float, DIM+1, 1> > & spaceCorrs, int samplingRate, bool reverseCorrs) {
	
	typedef Eigen::Matrix<Float, DIM+1, 1> VectorN;
	fstream f(filename.c_str(), ios::in);
	Json::Value root;  
	Json::Reader reader;
	string fileContents;
	bool readSuccess = readFileIntoString(filename, fileContents);


	if (!readSuccess) {
		LOG(ERROR) << "Failed to read ucorr file into string : " << filename;
		return false;
	}

	bool parsingSuccessful = reader.parse(fileContents, root );
	if ( !parsingSuccessful )
	{
		// report to the user the failure and their locations in the document.
		LOG(ERROR) << "Failed to parse ucorr file: \n "
			<< reader.getFormattedErrorMessages();
		return false;
	}

	bool useGlobalIndexing = false;

	if (!root["useGlobalIndexing"].isNull()) {
		useGlobalIndexing = root["useGlobalIndexing"].asBool();
	}

	LOG(INFO) << "Using global indexing for point corrs:" << useGlobalIndexing;

	if (root["pointCorrs"].isNull()) {
		LOG(WARNING) << "No point correspondences found";
	}
	else {
		Json::Value jsonPointCorrs = root["pointCorrs"];
		int numPointCorrs = jsonPointCorrs.size();
		for (int i = 0; i < numPointCorrs; i++) {

			if (useGlobalIndexing) {
				int indexA = jsonPointCorrs[i]["indexA"].asInt();
				int indexB = jsonPointCorrs[i]["indexB"].asInt();

				if (reverseCorrs)
					swap(indexA, indexB);

				pointCorrs[indexA] = indexB;
			}
			else {
				int timeIndexA = jsonPointCorrs[i]["timeIndexA"].asInt();
				int spaceIndexA = jsonPointCorrs[i]["spaceIndexA"].asInt();
				int timeIndexB = jsonPointCorrs[i]["timeIndexB"].asInt();
				int spaceIndexB = jsonPointCorrs[i]["spaceIndexB"].asInt();


				if (reverseCorrs) {
					swap(timeIndexA, timeIndexB);
					swap(spaceIndexA, spaceIndexB);
				}

				int indexA = animA->frameToGlobalMap[FrameIndex(timeIndexA,spaceIndexA)];
				int indexB = animB->frameToGlobalMap[FrameIndex(timeIndexB,spaceIndexB)];
				
				if (!jsonPointCorrs[i]["isSpaceCorr"].isNull() && jsonPointCorrs[i]["isSpaceCorr"].asBool()) {

					// Convert single space corr from the file into a set of correspondences over the entire animation (no time limiting on this, but can be implemented if needed)
					VectorN b = animB->vertices[indexB] - animA->vertices[indexA];
					VectorN posA = animA->vertices[indexA];
					//b[DIM] = 0.0;
					LOG(INFO) << "Found space corr:" << timeIndexA << "," << spaceIndexA << " -> " << timeIndexB << "," << spaceIndexB;

					int step = 2*pow(samplingRate, 0.333);
					for (int t = timeIndexA; t < animA->frames.size(); t+=step) {
						Float time = animA->getDepth(t);
						VectorN pos = posA;
						pos[DIM] = time;
						
						VLOG(1) << "Original position:" << posA;
						VLOG(1) << "Searching near:" << pos;


																								
						// Find a point close to the user specified location in each sampled frame
						int pointIndex = animA->findClosestVertex(pos);			
						VLOG(1) << "Found:" << animA->vertices[pointIndex];

						VLOG(1) << "Distance:" << (animA->vertices[pointIndex]-pos).norm();

						if ((animA->vertices[pointIndex]-pos).norm() < 5 * Globals::AvgEdgeLength) {							
							//spaceCorrs[pointIndex] = b;
							pointCorrs[pointIndex] = animB->findClosestVertex(pos + b);

							FrameIndex fIndex = animA->globalToFrameMap[pointIndex];
							LOG(INFO) << "Creating additional space corrs at:" << fIndex.first << "," << fIndex.second << " " << pointIndex;
						}

						
						
					}
										
				}
				else if (!jsonPointCorrs[i]["isLineCorr"].isNull() && jsonPointCorrs[i]["isLineCorr"].asBool()) {
					int lengthA = jsonPointCorrs[i]["lengthA"].asInt();
					int lengthB = jsonPointCorrs[i]["lengthB"].asInt();


					if (reverseCorrs) {
						swap(lengthA, lengthB);						
					}

					float ratio = (1.0*lengthB)/(1.0*lengthA);

					vector<GlobalIndex> lineA, lineB;

					lineA.push_back(indexA);
					lineB.push_back(indexB);

					int curIndex = indexA;					
					bool validLine = true;
					
					LOG(INFO) << "Creating line corrs:" <<  timeIndexA << "," << spaceIndexA << " with length:" << lengthA << " -> " << timeIndexB << "," << spaceIndexB << " with length:" << lengthB ;

					// Trace through the velocity field of A
					for (int k = 0; k < lengthA; k++) {
						//LOG(INFO) << "Tracing A:" << curIndex; 
						int nextIndex = animA->neighbours[curIndex][NextTime].index;
						if (nextIndex != Globals::InvalidNeighbour) {
							lineA.push_back(nextIndex);
							curIndex = nextIndex;
						}
						else
						{
							LOG(ERROR) << "Could not trace line A starting from:" << timeIndexA << "," << spaceIndexA << " with length:" << lengthA  << " ended at " << timeIndexA+k;
							validLine = false;
							break;
						}
					}

					// Trace through velocity field of B
					curIndex = indexB;
					for (int k = 0; k < lengthB; k++) {
						//LOG(INFO) << "Tracing B:" << curIndex; 
						int nextIndex = animB->neighbours[curIndex][NextTime].index;
						if (nextIndex != Globals::InvalidNeighbour) {
							lineB.push_back(nextIndex);
							curIndex = nextIndex;
						}
						else
						{
							LOG(ERROR) << "Could not trace line B starting from:" << timeIndexB << "," << spaceIndexB << " with length:" << lengthB << " ended at " << timeIndexB+k;
							validLine = false;
							break;
						}
					}

					int lineIndexA = 0;
					int lineIndexB = 0;
					FrameIndex fIndexA, fIndexB;	
					pointCorrs[lineA[lineIndexA]] = lineB[lineIndexB];
					fIndexA = animA->globalToFrameMap[lineA[lineIndexA]];
					fIndexB = animB->globalToFrameMap[lineB[lineIndexB]];
					LOG(INFO) << "Created point corr:" << fIndexA.first << "," << fIndexA.second << " -> " << fIndexB.first << "," << fIndexB.second;

					if (validLine) {
						// Create point corrs on first and last points of the line because preserving length is important
											
						
						int s = (int)(2*pow(samplingRate, 0.333)); // Floor the sampling rate

						LOG(INFO) << "Sampling rate:" << s;
						for (int lineIndexA = s; lineIndexA < lengthA - s + 1; lineIndexA += s) {
							int lineIndexB = lineIndexA * ratio; 

							pointCorrs[lineA[lineIndexA]] = lineB[lineIndexB];

							fIndexA = animA->globalToFrameMap[lineA[lineIndexA]];
							fIndexB = animB->globalToFrameMap[lineB[lineIndexB]];
							LOG(INFO) << "Created point corr:" << fIndexA.first << "," << fIndexA.second << " -> " << fIndexB.first << "," << fIndexB.second;
						}						

						pointCorrs[lineA[lineA.size()-1]] = lineB[lineB.size()-1]; // End points

						fIndexA = animA->globalToFrameMap[lineA[lineA.size()-1]];
						fIndexB = animB->globalToFrameMap[lineB[lineB.size()-1]];
						LOG(INFO) << "Created point corr:" << fIndexA.first << "," << fIndexA.second << " -> " << fIndexB.first << "," << fIndexB.second;


					}
					
				}
				else {
					pointCorrs[indexA] = indexB;
				}
				
			}
		}
		LOG(INFO) << "# Point Corrs:" << pointCorrs.size();
		LOG(INFO) << "# Space Corrs:" << spaceCorrs.size();
	}

	/* 
	// Deprecated these kinds of line features in 3D

	Json::Value linesA = root["linesA"];
	if (linesA.isNull()) {
		LOG(WARNING) << "No line features in animation A!";
	}
	else {
		int numLines = linesA.size();
		for (int i = 0; i < numLines; i++) {		
			LineFeature line;
			for (int j = 0; j < linesA[i]["indices"].size(); j++) {				
				line.push_back(linesA[i]["indices"][j].asInt());
			}		
			animA->lineFeatures.push_back(line);
		}
	}

	Json::Value linesB = root["linesB"];
	if (linesA.isNull()) {
		LOG(WARNING) << "No line features in animation B!";
	}
	else {
		int numLines = linesB.size();
		for (int i = 0; i < numLines; i++) {		
			LineFeature line;
			for (int j = 0; j < linesB[i]["indices"].size(); j++) {				
				line.push_back(linesB[i]["indices"][j].asInt());
			}			
			animB->lineFeatures.push_back(line);
		}
	}

	Json::Value jsonLineCorrs = root["lineCorrs"];
	if (jsonLineCorrs.isNull()) {
		LOG(WARNING) << "No line correspondences found!";
	}
	else {
		for (int i = 0; i < jsonLineCorrs.size(); i++) {
			int indexA = jsonLineCorrs[i]["indexA"].asInt();
			int indexB = jsonLineCorrs[i]["indexB"].asInt();
			int offset = jsonLineCorrs[i]["offset"].asInt();

			LineCorr corr(indexB, offset);
			if (lineCorrs.find(indexA) == lineCorrs.end()) {
				lineCorrs[indexA] = vector<LineCorr>();
			}
			lineCorrs[indexA].push_back(corr);			
		}
	}
	*/


	return true;
}