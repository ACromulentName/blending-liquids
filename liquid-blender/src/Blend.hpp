#pragma once

#include "stdafx.h"
#include "Mesh2D.h"
#include "Animation.h"
#include <Bench/BenchTimer.h>
#include "zlib.h"
#include "Output.hpp"

// Extract frames given an interpolated animation (writes bobjs to disk)
template<int DIM>
void extractFrames(const unique_ptr<Animation<DIM>>& interpolatedAnim, string path) {
		// Perturb vertices if they lie on frame times to prevent edge cases 
		//interpolatedAnim->smooth(5);

		interpolatedAnim->perturbVertices();
		
		// Compute time bounds for each frame		
		interpolatedAnim->computeTimeBounds();
		interpolatedAnim->setupExtraction();

		int numFrames = interpolatedAnim->frames.size();
		Eigen::BenchTimer timer;
		timer.reset();
		timer.start();
		

#ifdef ENABLE_PARALLEL		
		tbb::parallel_for(size_t(0), size_t(numFrames), size_t(1), [&](int i)						
		//for (int i = 0; i < numFrames; i++)
#else
		for (int i = 0; i < numFrames; i++)
#endif
		{
			cout << "Extracting frame:" << i << endl;

			char framePath[1024];
			sprintf_s(framePath, 1024, "%s/frame_%04d", path.c_str(), i);

			/*
			SDF<DIM> sdf;
			sdf.computeSDF(interpolatedAnim, time, 0.5 * Globals::AvgEdgeLength);			
			MarchingSquares<DIM> ms(sdf);
			ms.extract(framePath);
			*/

			Float time = interpolatedAnim->getDepth(i);
			interpolatedAnim->extractSlice(framePath, time);
		
		}
#ifdef ENABLE_PARALLEL
		);

#endif
		
		timer.stop();

		cout << "Done! Total: " << timer.total() << endl;
}

// Blends using two correspondences
template<int DIM>
bool blend(string corrFiles[2], Float baryU, Float baryV, Float baryW) {
	fstream f;
	typedef Eigen::Matrix<Float, DIM+1, 1> VectorN;
	

	Configuration<DIM> cfg[2];
	vector<VectorN> corrAnims[2];

	gzFile gzf;

	for (int i = 0 ; i < 2; i++) {
		gzf = gzopen(corrFiles[i].c_str(), "rb");		
		int length;

		gzread(gzf, &length, sizeof(int));
		LOG(INFO) << "Reading path with length:" << length;
		char cfgPath[2048];
		gzgets(gzf, cfgPath, length+1);

		LOG(INFO) << "Cfg path:" << cfgPath;		

		bool success = cfg[i].load(cfgPath, true);
		if (!success) 
			return false;
				
		readCorrFile<DIM>(gzf, corrAnims[i]);		
		gzclose(gzf);
	}
	
	unique_ptr<Animation<DIM>> interpolatedAnim(cfg[0].animA->makeCopy());		

	for (int i = 0; i < interpolatedAnim->vertices.size(); i++) {
		interpolatedAnim->vertices[i] = baryU * cfg[0].animA->vertices[i] + baryV * corrAnims[0][i] + baryW * corrAnims[1][i];
	}

	char path[1024];
	sprintf_s(path, 1024, "blend/%s-bary", cfg[0].configName.c_str());
	boost::filesystem::path dir(path);
	boost::filesystem::create_directory(dir);


	extractFrames(interpolatedAnim,path);

	return true;
}

// Blends using current registration
template<int DIM>
bool blend(Configuration<DIM>& cfg, const unique_ptr<Animation<DIM>>& corrAnim, int blendValue) {
	fstream f;
	typedef Eigen::Matrix<Float, DIM+1, 1> VectorN;
	
	try {

		Float blendFactor = blendValue * 0.01;
		 const unique_ptr<Animation<DIM>> interpolatedAnim (cfg.animA->makeCopy());


		for (int i = 0; i < interpolatedAnim->vertices.size(); i++) {
			interpolatedAnim->vertices[i] = (1-blendFactor) * cfg.animA->vertices[i] + (blendFactor) * corrAnim->vertices[i];
		}
		
		int numFrames = cfg.animA->frames.size();
		if (numFrames < cfg.animB->frames.size())
			numFrames = cfg.animB->frames.size();

		char path[1024];
		sprintf_s(path, 1024, "blend/%s-%03d", cfg.configName.c_str(), blendValue);
		boost::filesystem::path dir(path);
		boost::filesystem::create_directory(dir);

		cout << "Extracting frames ... ";
		extractFrames(interpolatedAnim, path);


	}
	catch(exception e) {
		return false;
	}
	return true;
} 


template<int DIM>
bool blend(string corrFile, int blendValue) {
	fstream f;
	typedef Eigen::Matrix<Float, DIM+1, 1> VectorN;
	
	try {
		//f.open(corrFile, ios::in);

		gzFile gzf;
		gzf = gzopen(corrFile.c_str(), "rb");		
		int length;

		gzread(gzf, &length, sizeof(int));
		LOG(INFO) << "Reading path with length:" << length;
		char cfgPath[2048];
		gzgets(gzf, cfgPath, length+1);

		LOG(INFO) << "Cfg path:" << cfgPath;
		Configuration<DIM> cfg;

		bool success = cfg.load(cfgPath, true);
		if (!success) 
			return false;

		vector<VectorN> corrVertices;
		readCorrFile<DIM>(gzf, corrVertices);


		Float blendFactor = blendValue * 0.01;

		// Make a copy of animA		
		unique_ptr<Animation<DIM>> interpolatedAnim(cfg.animA->makeCopy());		
		
		for (int i = 0; i < interpolatedAnim->vertices.size(); i++) {
			interpolatedAnim->vertices[i] = (1-blendFactor) * cfg.animA->vertices[i] + (blendFactor) * corrVertices[i];
		}
		//interpolatedAnim->computeNormals();
		//interpolatedAnim->createKDTree();

		int numFrames = cfg.animA->frames.size();
		if (numFrames < cfg.animB->frames.size())
			numFrames = cfg.animB->frames.size();

		char path[1024];
		sprintf_s(path, 1024, "blend/%s-%03d", cfg.configName.c_str(), blendValue);
		boost::filesystem::path dir(path);
		boost::filesystem::create_directory(dir);

		cout << "Extracting frames ... ";

		// Perturb vertices if they lie on frame times to prevent edge cases 
		interpolatedAnim->perturbVertices();

		
		// Compute time bounds for each frame
		cfg.animA->computeTimeBounds();
		interpolatedAnim->computeTimeBounds();

		interpolatedAnim->setupExtraction();

		Eigen::BenchTimer timer;
		timer.reset();
		timer.start();
		

#ifdef ENABLE_PARALLEL		
		tbb::parallel_for(size_t(0), size_t(numFrames), size_t(1), [&](int i)						
		//for (int i = 0; i < numFrames; i++)
#else
		for (int i = 0; i < numFrames; i++)
#endif
		{
			cout << "Extracting frame:" << i << endl;
			char framePath[1024];
			sprintf_s(framePath, 1024, "%s/frame_%04d", path, i);

			/*
			SDF<DIM> sdf;
			sdf.computeSDF(interpolatedAnim, time, 0.5 * Globals::AvgEdgeLength);			
			MarchingSquares<DIM> ms(sdf);
			ms.extract(framePath);
			*/

			Float time = interpolatedAnim->getDepth(i);
			interpolatedAnim->extractSlice(framePath, time);
		
		}
#ifdef ENABLE_PARALLEL
		);

#endif
		
		timer.stop();

		cout << "Done! Total: " << timer.total() << endl;



	}
	catch(exception e) {
		return false;
	}
	return true;
} 