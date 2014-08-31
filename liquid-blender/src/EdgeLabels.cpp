#include "stdafx.h"
#include "EdgeLabels.h"


bool findEdge(const EdgeLabels &edgeLabels, int a, int b, bool& value) {

	if (edgeLabels.find(Edge(a,b))!=edgeLabels.end()) {
		value = edgeLabels.at(Edge(a,b));
		return true;
	} 
	else if (edgeLabels.find(Edge(b,a))!=edgeLabels.end()){
		value = !edgeLabels.at(Edge(b,a));
		return true;
	}
	return false;

}

bool setEdge(EdgeLabels &edgeLabels, int a, int b, bool value) {

	if (edgeLabels.find(Edge(a,b))!=edgeLabels.end()) {
		edgeLabels[Edge(a,b)] = value;
		return true;
	} 
	else if (edgeLabels.find(Edge(b,a))!=edgeLabels.end()){
		edgeLabels[Edge(b,a)] = !value;
		return true;
	}
	return false;

}

bool isBadFace(const Face3D& face, const EdgeLabels& edgeLabels) {
	bool tmpLabels[3];

	for (int k = 0; k < 3; k++) {
		findEdge(edgeLabels,face.v[k], face.v[(k+1)%3],tmpLabels[k]);				
	}

	return (tmpLabels[0] == tmpLabels[1] && tmpLabels[1] == tmpLabels[2]);

}

// Stochastic
int labelEdgesStochastic(const Mesh3D& mesh, EdgeLabels& edgeLabels, int numIterations) {
	cout << "Labeling edges" << endl;		

	// Create an arbitrary labeling for each edge
	for (int j = 0; j < mesh.faces.size(); j++) {

		for (int k = 0; k < 3; k++) {
			bool value;
			if (!findEdge(edgeLabels, mesh.faces[j].v[k], mesh.faces[j].v[(k+1)%3], value)) {
				edgeLabels[Edge(mesh.faces[j].v[k], mesh.faces[j].v[(k+1)%3])] = rand()%2;
			}

		}		
	}
	vector<int> badFaces;
	for (int i = 0; i < numIterations; i++) {		

		badFaces.clear();
		for (int j = 0; j < mesh.faces.size(); j++) {			
			if (isBadFace(mesh.faces[j], edgeLabels))
				badFaces.push_back(j);					
		}

		if (i == numIterations-1 || badFaces.size() == 0) break;

		for (int k = 0; k < badFaces.size(); k++) {
			if (isBadFace(mesh.faces[badFaces[k]],edgeLabels)) {
				int start = rand()%3;			
				Edge e (mesh.faces[badFaces[k]].v[start], mesh.faces[badFaces[k]].v[(start+1)%3]);
				bool label;
				findEdge(edgeLabels, e.first, e.second, label);
				setEdge(edgeLabels, e.first, e.second, !label);
			}			
		}
	}
	LOG(INFO) << "# bad faces:" << badFaces.size() << " / " << mesh.faces.size();
	return badFaces.size();
}

void labelEdges(const Mesh3D& mesh, EdgeLabels& edgeLabels, int numIterations) {
	for (int j = 0; j < mesh.faces.size(); j++) {
		for (int k = 0; k < 3; k++) {
			bool value;
			if (!findEdge(edgeLabels, mesh.faces[j].v[k], mesh.faces[j].v[(k+1)%3], value)) {

				bool label = false;
				if (mesh.faces[j].v[k] > mesh.faces[j].v[(k+1)%3])
					label = true;
				edgeLabels[Edge(mesh.faces[j].v[k], mesh.faces[j].v[(k+1)%3])] = label;
			}

		}		
	}
}