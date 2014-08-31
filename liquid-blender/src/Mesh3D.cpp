#include "stdafx.h"
#include "Mesh3D.h"
#include "Utility.hpp"
#include <algorithm>
#include "zlib.h"

using namespace std;

void Mesh3D::cleanup() 
{
	Mesh<3>::cleanup();
	vertices.clear();
	faces.clear();
}


bool Mesh3D::loadBinaryObj( const string filename, Float scale) {
	
	LOG(INFO) << "Mesh3D::loadBinaryObj: Processing file:" << filename;

	gzFile gzf;
	int numVerts, numNorms, numTris;

	// open file
	gzf = gzopen(filename.c_str(), "rb1");
	if (!gzf) {
		LOG(ERROR) << "Cannot open file " << filename;
		return false;
	}
	
	// read vertices
	if (gzread(gzf, &numVerts, 4) != 4) {
		LOG(ERROR) << "Cannot read file " << filename;
		return false;
	}

	LOG(INFO) << "Reading " << numVerts << " vertices";

	for (int i = 0; i<numVerts; i++) {
		Vertex3D v;
		for (int j=0; j<3; j++) {
			float vp;
			if (gzread(gzf, &vp, sizeof(float)) != sizeof(float)) {
				LOG(ERROR) << "Error in file " << filename;
				return false;
			}
			v.pos[j] = vp;
		}
		vertices.push_back(v);
	}
	
	// read normals
	if (gzread(gzf, &numNorms, 4) != 4) {
		LOG(ERROR) << "Error reading normals: " << filename;
		return false;
	}
	
	LOG(INFO) << "Reading " << numNorms << " normals";

	for (int i = 0; i< numNorms; i++) {
		float v[3];
		for (int j=0; j<3; j++) {
			if (gzread(gzf, &v[j], sizeof(float)) != sizeof(float)) {
				LOG(ERROR) << "Error reading normal: " << i;
				return false;			
			}
		}		
	}	

	// read triangles
	if (gzread(gzf, &numTris, 4) != 4) {
		LOG(ERROR) << "Error reading faces: " << filename;
		return false;
	}

	
	for (int i=0; i<numTris; i++) {
		Face3D face;
		

		for (int j=0; j<3; j++) {			
			if (gzread(gzf, &face.v[j], 4) != 4) {
				LOG(ERROR) << "Error reading face: " << i;
				return false;
			}			
			
		}
		
		auto generateFaceHash = [] (int a, int b, int c) -> size_t {
			size_t seed = 0;
			boost::hash_combine(seed, a);
			boost::hash_combine(seed, b);
			boost::hash_combine(seed, c);
			return seed;
		};

		faceHashMap.insert(generateFaceHash(face.v[0],face.v[1],face.v[2]));
		faceHashMap.insert(generateFaceHash(face.v[0],face.v[2],face.v[1]));
		faceHashMap.insert(generateFaceHash(face.v[1],face.v[0],face.v[2]));
		faceHashMap.insert(generateFaceHash(face.v[1],face.v[2],face.v[0]));
		faceHashMap.insert(generateFaceHash(face.v[2],face.v[1],face.v[0]));
		faceHashMap.insert(generateFaceHash(face.v[2],face.v[0],face.v[1]));

		faces.push_back(face);

	}	
	gzclose(gzf);

	LOG(INFO) << "#Faces:" << numTris << " #Verts:" << numVerts;
	
	// Compute components
	buildCornerTable();

	// Now do an O(n) depth-first traversal to find all components
	vector<bool> visited(vertices.size());
	componentIds.reserve(vertices.size());

	int componentId = -1;

	for (int i = 0; i < vertices.size(); i++) {
		if (visited[i]) 
			continue;

		// New component
		componentId++;

		vector<int> curComponentVertices;

		queue<int> Q;
		Q.push(i);

		while (!Q.empty()) {
			int nextIndex = Q.front();
			Q.pop();

			if (!visited[nextIndex]) {
				visited[nextIndex] = true;
				curComponentVertices.push_back(nextIndex);
				componentIds[nextIndex] = componentId;
				vector<int> oneRing = getOneRing(nextIndex);
				for (int j = 0; j < oneRing.size(); j++) {
					Q.push(oneRing[j]);
				}
			}
		}

		const int minVertices = 300;
		if (curComponentVertices.size() < minVertices) {
			Droplet<3> droplet;
			droplet.vertices = curComponentVertices;
			Vector3 centroid(0,0,0);
			for (auto v : curComponentVertices) {
				centroid += vertices[v].pos;
			}
			droplet.centroid = centroid / curComponentVertices.size();
			droplets.push_back(droplet);
		}
	}

	LOG(INFO)<<"# components:" << componentId+1;

	numVertices = numVerts;
	return true;
}

bool Mesh3D::loadObj( const string filename )
{
	fstream file;
	try {
		LOG(INFO)<< "Loading mesh:" << filename;
		file.exceptions (  ifstream::badbit );
		file.open(filename.c_str(), ios::in);	
		if (!file) {
			return false;
		}

		int lineNumber = 0;
		while (!file.eof()) {
			lineNumber++;

			string line;
			getline(file, line);

			vector<string> elems = split(line, ' ', true);
			if (elems.size() <= 0)
				continue;

			if (elems[0] == "#") {
				continue;
			}
			else if (elems[0] == "v") {
				if (elems.size() != 4) {
					LOG(ERROR) << "Error in OBJ format on line:" << lineNumber; 
					file.close();
					return false;
				}
				Vertex3D v;
				for (int i = 1; i < elems.size(); i++) {					
					v.pos[i-1] = boost::lexical_cast<double>(elems[i]);
				}
				vertices.push_back(v);
				
			}
			else if (elems[0] == "vn") {
				continue;
			}
			else if (elems[0] == "vp") {
				continue;
			}
			else if (elems[0] == "f") {
				if (elems.size() != 4) {
					LOG(ERROR) << "Error in OBJ format on line:" << lineNumber; 
					file.close();
					return false;
				}
				Face3D face;
				for (int i = 1; i < elems.size(); i++) {
					vector<string> faceElems = split(elems[i],'/');
					face.v[i-1] = boost::lexical_cast<int>(faceElems[0]) - 1; // Vertex
				}
				faces.push_back(face);
			}

		}
		file.close();

		buildCornerTable();
		LOG(INFO)<< "Loaded mesh: " << filename << " successfully.";
		LOG(INFO)<< "v:" << vertices.size() << " f:" << faces.size();
		numVertices = vertices.size();
		return true;
	}
	catch (exception e) {		
		file.close();
		LOG(ERROR) << "Unable to load mesh: " << filename << endl;
		return false;
	}
	

	return false;
}

// Load velocities for each vertex (assumes that the mesh has been loaded)
bool Mesh3D::loadVelocities(string filename) {
	//ifstream file;
	gzFile gzf = NULL;
	try {
		LOG(INFO)<< "Loading velocities for: " << filename;
		gzf = gzopen(filename.c_str(), "rb");
		//file.open(filename.c_str(), ios::in);
		//if (!file) return false;

		if (gzf == NULL) 
			return false;

		velocities.clear();

		for (int i = 0; i < vertices.size(); i++) {
			int tmpId;
			bool isTinyComponent = false;
			double vel[3];
			//file >> tmpId >> vel[0] >> vel[1] >> vel[2];

			gzread(gzf, &tmpId, sizeof(tmpId));
			gzread(gzf, &vertices[i].isTinyComponent, sizeof(bool));
			gzread(gzf, &vel[0], sizeof(double));
			gzread(gzf, &vel[1], sizeof(double));
			gzread(gzf, &vel[2], sizeof(double));

			Vector3 v;
			for (int j = 0; j < 3; j++) 
				v[j] = (Float)vel[j];

			velocities.push_back(v);
			idToIndexMap[tmpId] = i;
			vertices[i].id = tmpId;

			if (vertices[i].isTinyComponent) {
				//LOG(INFO) << "Vertex :" << i << " is part of a tiny component : " << componentIds[i];
			}
		}
		//file.close();
		gzclose(gzf);



		LOG(INFO)<< "Loaded velocities for: " << filename << " successfully.";
		return true;
	}
	catch (exception e) {
		gzclose(gzf);
		LOG(ERROR)<< "FAILED to load velocities for: " << filename << " successfully.";
		//file.close();
	}
	return false;
}

void Mesh3D::buildCornerTable()
{
	V.clear();
	O.clear();
	C.clear();


	int numCorners = 3 * faces.size(); // 3 corners per face
	V.resize(numCorners); 
	O.resize(numCorners);
	C.resize(vertices.size());

	for (int i = 0; i < numCorners; i++) {
		V[i] = i;
		O[i] = i;
	}
	// Map corners to vertices 
	for (int i = 0; i < faces.size(); i++) {
		for (int j = 0; j < 3; j++) {
			int c = 3*i + j;
			V[c] = faces[i].v[j]; // Because triangle T always has corners (3T, 3T+1, 3T+2)
			C[faces[i].v[j]] = c; 
		}		
	}

	// Create opposite table by inserting each edge and associated corner into a map
	// If the two edges match, the associated corners are opposite to each other.
	// Otherwise, O[c] = c for boundary corners.
	map < pair<int,int>, int > edgeToCornerMap; 
	for (int i = 0; i < numCorners; i++) {
		int v0 = v(n(i));
		int v1 = v(p(i));

		if (v0 > v1)
			swap(v0,v1);

		pair<int,int> edge(v0,v1);
		auto it = edgeToCornerMap.find(edge);
		if (it != edgeToCornerMap.end()) {
			O[i] = it->second;
			O[it->second] = i;
		}
		else {
			edgeToCornerMap[edge] = i;
		}

	}

	edgeToCornerMap.clear();
}

bool Mesh3D::writeBinaryObj(const string filename) {
	try {
		gzFile gzf;
		gzf = gzopen(filename.c_str(), "wb");
		if (!gzf) {
			LOG(ERROR) << "Failed to write binary obj:" << filename;
			return false;
		}

		

		int numVertices = vertices.size();
		gzwrite(gzf, &numVertices, sizeof(numVertices));
		for(size_t i = 0; i < vertices.size(); i++) {
			for(int j = 0; j < 3; j++) {				
				float v = vertices[i].pos[j];
				gzwrite(gzf, &v, sizeof(v)); }
		}

		// Skip normals
		int numNormals = numVertices;
		gzwrite(gzf, &numNormals, sizeof(numNormals));


		// should be the same as Vertices.size		
		for(size_t i=0; i<vertices.size(); i++) {
			Vector3 normal(0,1,0);

			vector<int> oneRing = getOneRing(i);

			int numNbs = oneRing.size();
			
			if (oneRing.size() >= 2 ) {
				int count = 0;
				for (int k = 0; k < oneRing.size() - 1; k++) {
					Vector3 a = vertices[oneRing[k]].pos - vertices[i].pos;
					Vector3 b = vertices[oneRing[(k+1)%oneRing.size()]].pos - vertices[i].pos;				
					Vector3 tmp = b.cross(a);
					tmp.normalize();

					bool validNormal = true;
					if (count > 0) {
						Vector3 lastNormal = normal/count;
						lastNormal.normalize();
						if (tmp.dot(lastNormal) < 0)
							validNormal = false;
					}

					if (validNormal) {
						normal += tmp;
						count++;
					}
				}
				if (count > 0)
					normal /= count;

				normal.normalize();

			}

			for(int j=0; j<3; j++) {
				//float normp = normals[i][j]; // 1.; // mesh->surfaceNodes[i].n[j];
				float normp = normal[j];
				gzwrite(gzf, &normp, sizeof(normp)); }
		}

		int numTris = faces.size();

		gzwrite(gzf, &numTris, sizeof(numTris));
		for(size_t triIndex=0; triIndex < numTris; triIndex++) {
			for(int j = 0; j < 3; j++) { 
				int trip = faces[triIndex].v[j];
				gzwrite(gzf, &trip, sizeof(trip)); 
			} 
		}		
		gzclose( gzf );		
	}
	catch (exception e) {		
		LOG(ERROR) << "Error while writing binary obj";
		return false;
	}
	return true;
}

bool Mesh3D::writeObj(const string filename) {
	try {
		fstream f(filename, ios::out);
		for (int i = 0; i < vertices.size(); i++) {
			f << "v " << vertices[i].pos[0] << " " << vertices[i].pos[1] << " " << vertices[i].pos[2] << endl;
		}

		for (int i = 0; i < faces.size(); i++) {
			f << "f " << faces[i].v[0] + 1 << " " << faces[i].v[1] + 1 << " " << faces[i].v[2] + 1 << endl;
		}
		f.close();
	}
	catch (exception e) {
		return false;
	}
	return true;
}
