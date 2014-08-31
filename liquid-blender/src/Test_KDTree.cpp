#include "stdafx.h"
#include "kdTree.hpp"
#include "gtest/gtest.h"
#include <Bench/BenchTimer.h>

const int SAMPLES_DIM = 3;
void generateRandomPointCloud(Points3D &samples, const size_t N,const size_t dim, const double max_range = 10.0)
{
	std::cout << "Generating "<< N << " random points...";
	samples.resize(N);
	for (size_t i=0;i<N;i++)
	{
		samples[i].resize(dim);
		for (size_t d=0;d<dim;d++)
			samples[i][d] = max_range * (rand() % 1000) / (1000.0);
	}
	std::cout << "done\n";
}

void kdtree_demo(const size_t nSamples,const size_t dim)
{
	typedef KDTreeVectorOfVectorsAdaptor< std::vector < Vector3 >, Float, 3>  KDTree3D;

	std::vector<Vector3>  samples;

	const double max_range = 20;

	// Generate points:
	generateRandomPointCloud(samples, nSamples,dim, max_range);

	// Query point:
	std::vector<Float> query_pt(dim);
	for (size_t d=0;d<dim;d++)
		query_pt[d] = max_range * (rand() % 1000) / (1000.0);

	// construct a kd-tree index:
	// Dimensionality set at run-time (default: L2)
	// ------------------------------------------------------------
	

	KDTree3D   mat_index(dim , samples, 10 );
	mat_index.index->buildIndex();

	// do a knn search
	const size_t num_results = 3;
	std::vector<size_t>   ret_indexes(num_results);
	std::vector<Float> out_dists_sqr(num_results);

	nanoflann::KNNResultSet<Float> resultSet(num_results);

	resultSet.init(&ret_indexes[0], &out_dists_sqr[0] );
	mat_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

	std::cout << "knnSearch(nn="<<num_results<<"): \n";
	for (size_t i=0;i<num_results;i++)
		std::cout << "ret_index["<<i<<"]=" << ret_indexes[i] << " out_dist_sqr=" << out_dists_sqr[i] << std::endl;
}



TEST(KDTreeTest, RandomlyGeneratedPoints) {
	Eigen::BenchTimer timer;
	timer.start();
	kdtree_demo(1e3 /* samples */, SAMPLES_DIM /* dim */);
	timer.stop();
	LOG(INFO)<<"Time:"<<timer.total();
}