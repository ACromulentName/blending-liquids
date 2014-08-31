#include "stdafx.h"
#include "Retiming.h"
#include "Animation.h"
#include "Utility.hpp"

template<int DIM>
class Retimer {
	typedef Eigen::Matrix<Float, DIM+1, 1> VectorN;
public:
	Retimer(RetimingConfiguration<DIM>& cfg) {

		curAnim = unique_ptr<Animation<DIM>>(cfg.animA->makeCopy());

		// Compute alphas
		for (int i = 0; i < cfg.retimingNodes.size(); i++) {

			Float rSpace = cfg.retimingNodes[i].spaceRadius ;
			Float rTime = cfg.retimingNodes[i].timeRadius * cfg.animA->timeScaling;

			vector<AnisotropicNeighbour> nbs = doSplitBFS(cfg.animA, cfg.retimingNodes[i].vertices, rTime, rSpace);			

			for (auto nb : nbs) {
				VectorN nbVertex = cfg.animA->vertices[nb.index];
				Float dTime = 1e6;
				for (int k = 0; k < cfg.retimingNodes[i].vertices.size(); k++) {
					VectorN vertex = cfg.animA->vertices[cfg.retimingNodes[i].vertices[k]];
					dTime = min(abs(nbVertex[DIM] - vertex[DIM]),dTime);
				}
				
				Float weightSpace = falloffKernel(nb.dSpace, rSpace);
				Float weightTime = falloffKernel(nb.dTime, rTime);
				
				Float netWeight = sqrt(weightSpace*weightTime);
				netWeight = 0.5 * (weightSpace + weightTime);
				Float alpha = nb.d; // This is the contribution of the current deformation node 

				FrameIndex fIndex = cfg.animA->globalToFrameMap.at(nb.index);
				if (netWeight > 1e-6)
					LOG(INFO) << "Vertex:" << fIndex.first << "," << fIndex.second << " : " << nb.dSpace << " " << weightSpace << " | " << nb.dTime << " " << weightTime; 

				//alpha = 1.0f / (1+alpha);
				alpha = 1.0;
				if (alphas.find(nb.index) == alphas.end()) {
					alphas[nb.index] = alpha;

					for (int k = 0; k < DIM; k++)
						curAnim->vertices[nb.index][k] =  cfg.animA->vertices[nb.index][k];					

					curAnim->vertices[nb.index][0] = alpha * netWeight * cfg.retimingNodes[i].xDisplacement;

					curAnim->vertices[nb.index][DIM] = alpha * netWeight * cfg.animA->getDepth(cfg.retimingNodes[i].timeDisplacement);

					FrameIndex fIndex = cfg.animA->globalToFrameMap.at(nb.index);
					//if (netWeight > 0)
					
				}
				else {
					auto curSumAlphas = alphas[nb.index];
					alphas[nb.index] = curSumAlphas + alpha;

					curAnim->vertices[nb.index][0] += alpha * netWeight * cfg.retimingNodes[i].xDisplacement;
					curAnim->vertices[nb.index][DIM] += alpha * netWeight * cfg.animA->getDepth(cfg.retimingNodes[i].timeDisplacement);
				}
			}			
		}

		for (auto it = alphas.begin(); it != alphas.end(); it++) {			
			//LOG(INFO) << "Net alpha:" << it->first << " : " << it->second;
			{
				Float delta = curAnim->vertices[it->first][0];
				if (it->second > 0)
					delta /= it->second;
				curAnim->vertices[it->first][0] = delta + cfg.animA->vertices[it->first][0];

			}
			{
				Float delta = curAnim->vertices[it->first][DIM];
				if (it->second > 0)
					delta /= it->second;
				curAnim->vertices[it->first][DIM] = delta + cfg.animA->vertices[it->first][DIM];

			}
		}



		// Now extract new frames

		// Perturb vertices if they lie on frame times to prevent edge cases 
		//curAnim->perturbVertices();

		
		// Compute time bounds for each frame
		cfg.animA->computeTimeBounds();
		curAnim->computeTimeBounds();

		curAnim->setupExtraction();

		string path = "retiming";

#ifdef ENABLE_PARALLEL		
		tbb::parallel_for(size_t(0), size_t(cfg.animA->numFrames), size_t(1), [&](int i)						
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

			Float time = curAnim->getDepth(i);
			curAnim->extractSlice(framePath, time);
		
		}
#ifdef ENABLE_PARALLEL
		);

#endif
		
		// Write out positions
		fstream f;
		f.open("retiming/0.corr", ios::out);

		boost::filesystem::path p = boost::filesystem::canonical(boost::filesystem::path(cfg.configPath));
		
		f << p.string() << endl;
		f << curAnim->vertices.size() << endl;

		for (int i = 0; i < curAnim->vertices.size(); i++) {
			f << curAnim->vertices[i][0] << " " << curAnim->vertices[i][1] << " " << curAnim->vertices[i][2] << endl;
		}
		f.close();
	}

	~Retimer() {

	}

	void retime() {

	}

	Float falloffKernel(Float d, Float maxRadius)  {
		if (d > maxRadius) return 0;
		return pow(1.0 - pow(d/maxRadius,2),2);
		//return 1.0f - d/maxRadius;
	}

	unique_ptr<Animation<DIM>> curAnim;
	unordered_map<int, Float > alphas; // Contributions of each deformation node to a vertex
};



