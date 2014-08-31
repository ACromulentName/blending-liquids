
#pragma once

#include "stdafx.h"
#include "Mesh2D.h"
#include "Animation.h"
#include <queue>
#include <set>
#include <fstream>
#include "boost/filesystem.hpp"
#include <functional>
#include "json.h"
using namespace std;

struct RetimingNode {
	vector<int> vertices;
	Float timeDisplacement;
	Float xDisplacement;
	Float yDisplacement;
	Float spaceRadius;
	Float timeRadius;
};

template<int DIM>
struct RetimingConfiguration {
	typedef Eigen::Matrix<Float, DIM+1, 1> VectorN;
	typedef map < int, VectorN > SpaceCorrs;

	RetimingConfiguration() {

	}

	RetimingConfiguration(string configFilename) : RetimingConfiguration() {		
		configPath = configFilename;
		load(configFilename);
	}

	bool load(string configFilename) {
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


		if (root["animA"].isNull()) {
			LOG(ERROR) << "No entry for animation A!";
			return false;
		}
		else
		{
			std::string path = root["animA"]["path"].asString();
			int startFrame = root["animA"]["startFrame"].asInt();
			int endFrame = root["animA"]["endFrame"].asInt();

			bool success = animA->load(path, startFrame, endFrame, false);
			if (!success)
				return false;
		}


		if (root["retimingNodes"].isNull()) {
			LOG(ERROR) << "No retiming nodes found!";
			return false;
		}
		else {
			for (int i = 0; i < root["retimingNodes"].size(); i++) {
				RetimingNode rn;
				
				rn.xDisplacement = root["retimingNodes"][i]["xDisplacement"].asFloat(); 
				rn.timeDisplacement = root["retimingNodes"][i]["timeDisplacement"].asFloat(); // Number of frames
				rn.spaceRadius = root["retimingNodes"][i]["spaceRadius"].asFloat(); // Currently use number of vertices (maybe switch to edge length based version)
				rn.timeRadius = root["retimingNodes"][i]["timeRadius"].asFloat(); // Currently use number of frames

				for (int j = 0; j < root["retimingNodes"][i]["vertices"].size(); j++) {
					int timeIndex = root["retimingNodes"][i]["vertices"][j]["timeIndex"].asInt();
					int spaceIndex = root["retimingNodes"][i]["vertices"][j]["spaceIndex"].asInt();

					VLOG(1) << "Reading re-timing vertex: " << timeIndex << "," << spaceIndex;
					int gIndex = animA->frameToGlobalMap.at(FrameIndex(timeIndex,spaceIndex));
					rn.vertices.push_back(gIndex);					
				}
				retimingNodes.push_back(rn);
			}
		}

		return true;
	}

	vector<RetimingNode> retimingNodes;
	unique_ptr<Animation<DIM>> animA;
	string configPath;
	string configName;

private:
	RetimingConfiguration(const RetimingConfiguration<DIM>&) {}
	RetimingConfiguration& operator=(const RetimingConfiguration<DIM>&) { return *this;}
};

template<>
RetimingConfiguration<2>::RetimingConfiguration() {
	animA = unique_ptr<Animation2D>(new Animation2D());

}


template<>
RetimingConfiguration<3>::RetimingConfiguration() {
	animA = unique_ptr<Animation3D>(new Animation3D());

}