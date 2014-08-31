#pragma once

#include "stdafx.h"
#include "Mesh2D.h"
#include "Animation.h"
#include <fstream>
#include "boost/filesystem.hpp"
#include <functional>
#include "Configuration.hpp"
#include "zlib.h"
using namespace std;

template<int DIM> class ICP;

template<int DIM> 
void readCorrFile(gzFile& gzf,  vector<Eigen::Matrix<Float, DIM+1, 1>>& corrVertices) {
	typedef Eigen::Matrix<Float, DIM+1, 1> VectorN;
	int numVertices;
	//f >> numVertices;
	gzread(gzf, &numVertices, sizeof(int));
	corrVertices.reserve(numVertices);

	for (int i = 0; i < numVertices; i++) {
		VectorN tmpVertex;
		//f >> tmpVertex[0] >> tmpVertex[1] >> tmpVertex[2] >> tmpVertex[3];
		for (int j = 0; j < DIM+1; j++) 
			gzread(gzf, &tmpVertex[j], sizeof(Float));

		corrVertices.push_back(tmpVertex);
	}
	gzclose(gzf);
}

template<int DIM> void writeCorrFile(string filename, string configPath, Animation<DIM>& anim);


template<int DIM>
void saveError(Configuration<DIM>& cfg, ICP<DIM> &icp, int iter) {
	char filename[1024];
	sprintf(filename, "registrations/%s_%04d.error",cfg.configName.c_str(), iter);
	fstream f(filename, ios::out);

	boost::filesystem::path p = boost::filesystem::canonical(boost::filesystem::path(cfg.configPath));
	const vector<Error>& errors = icp.getErrors();

	f << p.string() << endl;	

	f << iter << endl;

	for (int i = 0; i < iter; i++) {
		f << i << " " << errors[i].avgSmooth << " " << errors[i].avgPoint << "  " << errors[i].avgPlane << " " << 
			errors[i].maxSmooth << " " << errors[i].maxPoint << "  " << errors[i].maxPlane << " " << 
			errors[i].numUnmatchedVertices << endl;
	}
	f.close();
}