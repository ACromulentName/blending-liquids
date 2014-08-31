#include "stdafx.h"
#include "Mesh2D.h"

#include <fstream>

using namespace std;

Mesh2D::Mesh2D(const string filename) {
	load(filename);
}

bool Mesh2D::load( const string filename ) {
	fstream file;
	try {
		LOG(INFO)<< "Loading mesh:" << filename;
		file.exceptions ( ifstream::eofbit | ifstream::failbit | ifstream::badbit );
		file.open(filename.c_str(), ios::in);	
		if (!file) {
			return false;
		}

		
		file >> numVertices;

		cleanup();

		for (int i = 0; i < numVertices; i++) {
			Vertex2D v;
			NeighbouringEdges nbEdges;
			file >> v.pos[0] >> v.pos[1] >> v.id >> nbEdges.prevEdge >> nbEdges.nextEdge;
			vertices.push_back(v);
			neighbouringEdges.push_back(nbEdges);
			idToIndexMap[v.id] = i;
		}

		int numEdges;
		file >> numEdges;

		for (int i = 0; i < numEdges; i++) {
			Edge2D e;
			int tmpId;
			file >> e.a >> e.b >> tmpId;
			edges.push_back(e);
		}
		file.close();
		
		LOG(INFO)<< "Loaded mesh: " << filename << " successfully.";
		return true;
	}
	catch (exception e) {		
		file.close();
		LOG(ERROR) << "Unable to load mesh:" << filename;
		return false;
	}

	return false;

}

// Load velocities for each vertex (assumes that the mesh has been loaded)
bool Mesh2D::loadVelocities(string filename) {
	ifstream file;
	try {
		LOG(INFO)<< "Loading velocities for: " << filename;
		file.open(filename.c_str(), ios::in);
		if (!file) return false;

		velocities.clear();

		for (int i = 0; i < vertices.size(); i++) {
			int tmpId;
			Vector2 vel;
			file >> tmpId >> vel[0] >> vel[1];
			velocities.push_back(vel);
		}
		file.close();

		LOG(INFO)<< "Loaded velocities for: " << filename << " successfully.";
		return true;
	}
	catch (exception e) {
		file.close();
	}
	return false;
}

// Clean up 
void Mesh2D::cleanup() {
	vertices.clear();
	edges.clear();
	velocities.clear();
	idToIndexMap.clear();	
}

Mesh2D::~Mesh2D() {
	cleanup();
}