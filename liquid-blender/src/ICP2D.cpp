#include "stdafx.h"
#include "ICP.hpp"

template<>
void ICP<2>::setupRadii() {
	vertexRadius = 1.25 * 2.414 * animA->timeScaling * pow(samplingRate,0.5);
	nodeRadius = 2 * vertexRadius;
}

template<>
void ICP<2>::widenLines() {
	linePointCorrs.clear();
	for (auto it = lineCorrs.begin(); it != lineCorrs.end(); it++) {
		int lineIndexA = it->first;

		vector <int> pointsLineB;
		int width = 2;

		for (int j = 0; j < it->second.size(); j++) {
			int lineIndexB = it->second[j].first;				

			for (auto lineIt = animB->lineFeatures[lineIndexB].begin(); lineIt != animB->lineFeatures[lineIndexB].end(); lineIt++) {
				GlobalIndex pointIndex = *lineIt;
				pointsLineB.push_back(pointIndex);

				int prevIndex = pointIndex;
				int nextIndex = pointIndex;
				for (int w = 0; w < width; w++) {
					pointsLineB.push_back(nextIndex);
					pointsLineB.push_back(prevIndex);
					nextIndex = dynamic_cast<Animation2D*>(animB.get())->getNextSpaceIndex(nextIndex);
					prevIndex = dynamic_cast<Animation2D*>(animB.get())->getPrevSpaceIndex(prevIndex);
				}
			}
		}

		linePointCorrs[lineIndexA] = pointsLineB;

	}
}


template<>
void ICP<2>::sampleUniformTime() {
	
	for (int time = 0; time < animA->frames.size(); time++) {
		// Create a node on this frame?
		if ((enforceEqualLength && ((time == 0) || (time == animA->frames.size() - 1))) ||
			(time % samplingRate == 0)) {

			VisitorMesh2D visitor(*dynamic_pointer_cast<const Mesh2D>(animA->frames[time]));
			int curIndex = 0;
			int count = 0;
			while(true) {
				curIndex = visitor.next();
				if (curIndex == -1) {
					break;
				}

				if (count % samplingRate == 0) {
					FrameIndex frameIndex(time, curIndex);
					auto it =  animA->frameToGlobalMap.find(frameIndex);
					GlobalIndex globalIndex = it->second;
					if (!isNode(globalIndex)) {
						// Create a node here
						createNode(globalIndex);
					}
				}
				count++;				
			}
		}
	}

}






template<>
bool ICP<2>::isNodeUserCorrespondence(int nodeIndex) {
	// TODO: Search point and line correspondences to see if node is part of them	
	GlobalIndex globalIndex = nodes[nodeIndex];

	if (isPointCorrespondence(globalIndex))
		return true;

	if (animA->isLineFeature(globalIndex)) {
		int lineIndex = animA->getLineFeature(globalIndex).first;
		return isLineCorrespondence(lineIndex);
	}

	return false;
}


// Specialization: Line features are unique to 2D
template<>
void ICP<2>::createNodesOnUserCorrs() {

	cout << " User nodes :" << endl;	
	for (auto it = pointCorrs.begin(); it != pointCorrs.end(); it++) {

		if (!isNode(it->first)) {
			createNode(it->first);
		}
	}


	for (auto it = lineCorrs.begin(); it != lineCorrs.end(); it++) {
		int count = 0;
		LineFeature line = animA->lineFeatures[it->first];
		for (auto lineIt = line.begin(); lineIt != line.end(); lineIt++) {
			if (!isNode(*lineIt) && (count % samplingRate == 0) ) {

				LOG(INFO) << "Creating line corr node at" << *lineIt ;
				createNode(*lineIt);				
			}
			count++;
		}
	}
	return;
}


template<> 
void ICP<2>::writeCorrFile(string filename, string configPath, bool writeProjected) {
	fstream f;
	f.open(filename, ios::out);
	f << configPath << endl;
	f << curAnim->vertices.size() << endl;

	for (int i = 0; i < curAnim->vertices.size(); i++) {
		f << curAnim->vertices[i][0] << " " << curAnim->vertices[i][1] << " " << curAnim->vertices[i][2] << endl;
	}
	f.close();

}