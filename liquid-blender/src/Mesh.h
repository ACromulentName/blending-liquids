#pragma once
#include "stdafx.h"

using namespace std;



template<int DIM>
class Vertex
{
	
	typedef Eigen::Matrix<Float, DIM, 1> VertexVector;

public:
	Vertex() {
		isTinyComponent = false;
	}

	VertexVector pos;
	int id;
	bool isTinyComponent;
};


template<int DIM>
struct Droplet {
	typedef Eigen::Matrix<Float, DIM, 1> VertexVector;
	VertexVector centroid;
	vector<int> vertices; // Local frame indices
};

template<int DIM>
class Mesh 
{
	typedef Eigen::Matrix<Float, DIM+1, 1> VectorN;
	typedef Eigen::Matrix<Float, DIM, 1> VelocityVector;
public:
	Mesh() {}
	//virtual bool load(const string filename);
	virtual bool loadVelocities(string filename) { return false; }
	virtual void cleanup() {}
	virtual ~Mesh() {}

	vector< Vertex<DIM> > vertices;
	vector<VelocityVector > velocities;	
	map < int, int > idToIndexMap; // Maps ids to vertex index 
	vector<int> componentIds; // stores component id of each vertex
	int numVertices;
	vector<Droplet<DIM>> droplets;

};