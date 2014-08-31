#include "stdafx.h"
#include "ICP.hpp"
#include <sstream>
#include "zlib.h"

template<>
void ICP<3>::setupRadii() {
	vertexRadius = 1.25 * 2.414 * animA->timeScaling * pow(samplingRate,0.333);
	nodeRadius = 1.6 * vertexRadius;
}

template<>
bool ICP<3>::isNodeUserCorrespondence(int nodeIndex) {
	GlobalIndex globalIndex = nodes[nodeIndex];

	if (isPointCorrespondence(globalIndex) || isSpaceCorrespondence(globalIndex))
		return true;

	return false;
}


// Specialization: Line features are unique to 2D
template<>
void ICP<3>::createNodesOnUserCorrs() {
	LOG(INFO) << " User nodes :";

	for (auto it = spaceCorrs.begin(); it != spaceCorrs.end(); it++) {

		if (!isNode(it->first)) {
			createNode(it->first);
			LOG(INFO)<<"Created space corr node at:"<<animA->vertices[it->first].format(CommaFmt);		
		}
	}
	for (auto it = pointCorrs.begin(); it != pointCorrs.end(); it++) {

		if (!isNode(it->first)) {
			createNode(it->first);

			LOG(INFO)<<"Created point corr node at:"<<animA->vertices[it->first].format(CommaFmt);		
		}
	}

}


template<>
void ICP<3>::widenLines() {
	return;
}

template<>
void ICP<3>::writeCorrFile(string filename, string configPath, bool projected) {
	fstream f;
	//f.open(filename, ios::out);

	gzFile gzf;

	gzf = gzopen(filename.c_str(), "wb");

	const char *configStr = configPath.c_str();
	int length = strlen(configStr);

	gzwrite(gzf, &length, sizeof(int));
	//gzwrite(gzf, &configStr, configPath.length());

	gzputs(gzf, configStr);

	//f << configPath << endl;
	int numVertices = curAnim->vertices.size();
	gzwrite(gzf, &numVertices, sizeof(int));
	//f << curAnim->vertices.size() << endl;

	for (int i = 0; i < curAnim->vertices.size(); i++) {
		for (int j = 0; j < 4; j++) {
			if (projected && corrs[i].valid)
				gzwrite(gzf, &corrs[i].p[j], sizeof(Float));
			else
				gzwrite(gzf, &curAnim->vertices[i][j], sizeof(Float));
		}
		//f << curAnim->vertices[i][0] << " " << curAnim->vertices[i][1] << " " << curAnim->vertices[i][2] << " " << curAnim->vertices[i][3] << endl;


	}
	//f.close();
	gzclose(gzf);

}