#pragma once
#include "stdafx.h"

class NeighbourSmooth {
public:
	virtual int getIndex() = 0;
	virtual Float getWeight() = 0;
};



struct Neighbour : public NeighbourSmooth {
	int index;
	Float d;

	Neighbour() {
		index = Globals::InvalidNeighbour;
		d = 0.0;
	}

	Neighbour(int _index, Float _d) {
		index = _index;
		d = _d;
	}

	virtual bool operator>(const Neighbour& rhs) const {
		return (d > rhs.d);
	}

	virtual bool operator<(const Neighbour& rhs) const {
		return (d < rhs.d);
	}

	virtual Float getWeight() {
		return d;
	}

	virtual int getIndex() {
		return index;
	}
};

struct AnisotropicNeighbour {
	int index;
	Float d;
	Float dSpace;
	Float dTime;

	AnisotropicNeighbour() {
		index = Globals::InvalidNeighbour;
		d = dTime = dSpace = 0.0;
	}



 	AnisotropicNeighbour(int _index,Float _d, Float _dTime, Float _dSpace) {
		index = _index;
		d = _d;
		dTime = _dTime;
		dSpace = _dSpace;
	}

	bool operator>(const AnisotropicNeighbour& rhs) const {
		return (d > rhs.d);
	}

	bool operator<(const AnisotropicNeighbour& rhs) const {
		return (d < rhs.d);
	}

	void getWeight(Float& dT, Float& dS) {
		dT = dTime;
		dSpace = dSpace;
	}

	int getIndex() {
		return index;
	}
};