#pragma once
#include "stdafx.h"
#include "Mesh.h"
using namespace std;


struct Face3D {
	int v[3]; // index of vertex (ccw)
	Vector3 n; // normal
	
};

typedef Vertex<3> Vertex3D;

class Mesh3D : public Mesh<3> {
public:
	Mesh3D():Mesh<3>() {}
	virtual ~Mesh3D() {}

	virtual bool loadVelocities( string filename);
	virtual void cleanup();


	bool loadObj(const std::string filename);
	bool loadBinaryObj(const string filename, Float scale = 1.0);
	bool writeObj(const string filename);
	bool writeBinaryObj(const string filename);



	vector < Face3D > faces;
	unordered_set<size_t> faceHashMap;


	void buildCornerTable();

	int c(int v) const {
		return C[v];
	}
	
	/** List of vertex ids on boundary that this corner faces */
	vector<int> boundaryVertices(int c) {
		vector<int> result;

		int count = 0;
		int startC = c;
		do {
			result.push_back((int)v(p(c)));
			c = cwToBoundary(c);
			count++;
			if (count > 50) {
				result.clear();
				break;
			}
		} while(c != startC);
		return result;
	}

	int cwToBoundary(int c) {
		int MAX_VALENCE = 50;
		int i = 0;

		c = p(c);
		while(!b(c) && i++ < MAX_VALENCE) c = p(o(c));
		return c;
	}

	int t(int c) const {
		return static_cast<int>(c/3);
	}

	int n(int c) const {
		return 3*t(c) + (c+1)%3; 	
	}
	
	int p(int c) const {
		return n(n(c));
	}

	int v(int c) const {
		return V[c];
	}

	Vertex3D g(int c) const{
		return vertices[v(c)];
	}

	inline bool b(int c) const {
		return (o(c) == c);
	}

	inline int o(int c) const {		
		return O[c];
	}
	
	inline int l(int c) const{
		return o(p(c));
	}

	inline int r(int c) const {
		return o(n(c));
	}

	inline int s(int c) const{
		return p(l(c));
	}
	
	inline int rs(int c) const{
		return n(r(c));
	}

	int visit(bool reset=false) {
		if (reset) {
			visited.clear();

			while (!visitQ.empty())
				visitQ.pop();

			for (int i = 0; i < faces.size(); i++) {
				visited.push_back(false);
			}
		}


		int curCorner = -1;

		while (!visitQ.empty()) {
			curCorner = visitQ.front();		
			visitQ.pop();
			if (!visited[t(curCorner)])
				break;
		} 

		if (visited[t(curCorner)]) 
			curCorner = -1;

		if (curCorner == -1 || visitQ.empty()) {
			int startFace = -1;
			for (int face = 0; face < faces.size(); face++) {
				if (!visited[face]) {
					startFace = face;
					break;
				}
			}

			if (startFace == -1)
				return -1;

			visitQ.push(c(faces[startFace].v[0]));
		}

				
		visited[t(curCorner)] = true;

		if (!visited[t(r(curCorner))])
			visitQ.push(r(curCorner));

		if (!visited[t(l(curCorner))])
			visitQ.push(l(curCorner));

		return t(curCorner);
	}

	vector<int> getOneRing(int vIndex) const {
		vector<int> oneRing;
		int startCorner = c(vIndex);		
		int curCorner = s(startCorner);
		

		// Have to be careful for meshes with boundaries
		int count = 0;
		bool hasBoundary = false;
		while (curCorner != startCorner) {						
			// Stop and swing in reverse direction
			if (v(curCorner) != v(startCorner)) {
				//
				hasBoundary = true;
				break;
			}

			oneRing.push_back(v(n(curCorner)));
			curCorner = s(curCorner);
			count++;
			
			if (count > 20) {
				oneRing.clear();
				break;
			}
			
		}

		if (hasBoundary) {
		curCorner = rs(startCorner);
			while (curCorner != startCorner) {						
				// Stop and swing in reverse direction
				if (v(curCorner) != v(startCorner)) 
					break;

				oneRing.push_back(v(p(curCorner)));
				curCorner = rs(curCorner);
				count++;

				if (count > 20) {
					oneRing.clear();
					break;
				}
			}
		}

		return oneRing;
	}
	
protected:

	// Corner table data structures	 
	vector<int> O, V, C;
	vector<bool> visited;
	queue<int> visitQ;

};






