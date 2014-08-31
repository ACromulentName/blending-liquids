///////////////////////////////////////
// Edge labeling utility functions
///////////////////////////////////////
#pragma once
#include "stdafx.h"
#include "Animation.h"
#include "Mesh3D.h"

// Labels each edge as 0 or 1 (0 => diagonal leaves vertex, 1=> diagonal is connected to future vertex) such that each triangle has 
// at least one 0 and one 1 edge.


bool findEdge(const EdgeLabels &edgeLabels, int a, int b, bool& value);

bool setEdge(EdgeLabels &edgeLabels, int a, int b, bool value);


bool isBadFace(const Face3D& face, const EdgeLabels& edgeLabels);

// Stochastic
int labelEdgesStochastic(const Mesh3D& mesh, EdgeLabels& edgeLabels, int numIterations = 50);

void labelEdges(const Mesh3D& mesh, EdgeLabels& edgeLabels, int numIterations = 50);