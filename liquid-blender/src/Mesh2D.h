#pragma once
#include "stdafx.h"

#include "Mesh.h"
using namespace std;

typedef Vertex<2> Vertex2D;

class NeighbouringEdges
{
public:

	int prevEdge;
	int nextEdge;
};


class Edge2D
{
public:
	int a;
	int b; 
};


class Mesh2D : public Mesh<2> 
{
public:
	Mesh2D() : Mesh<2>() { }
	Mesh2D(const string filename);	
	bool load(const string filename);
	virtual bool loadVelocities(string filename);
	virtual void cleanup();
	virtual ~Mesh2D();

	int nextVertexIndex(int index) const { 
		return edges[neighbouringEdges[index].nextEdge].b; 
	}

	int prevVertexIndex(int index) const { 
		return edges[neighbouringEdges[index].prevEdge].a; 
	}

	Vertex<2> const & nextVertex(int index) const {
		return vertices[nextVertexIndex(index)];
	}

	Vertex<2> const & prevVertex(int index) const {
		return vertices[prevVertexIndex(index)];
	}


	vector<NeighbouringEdges> neighbouringEdges;
	vector<Edge2D> edges;
};
