
#include "stdafx.h"
#include "Utility.hpp"
#include <sstream>
#include "ICP.hpp"

bool readFileIntoString(std::string& filename, std::string& data) {

	std::ifstream f;
	f.exceptions ( std::ifstream::eofbit | std::ifstream::failbit | std::ifstream::badbit );

	try {
		f.open(filename);		
		f.seekg(0, std::ios::end);   
		data.reserve(f.tellg());
		f.seekg(0, std::ios::beg);
		data.assign((std::istreambuf_iterator<char>(f)),
			std::istreambuf_iterator<char>());
		return true;
	}
	catch (std::exception e) {
		return false;
	}

}


std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems, bool removeEmptyStrings) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		if (removeEmptyStrings && item.length() <= 0);
		else
			elems.push_back(item);
	}
	return elems;
}


std::vector<std::string> split(const std::string &s, char delim, bool removeEmptyStrings) {
	std::vector<std::string> elems;
	split(s, delim, elems, removeEmptyStrings);
	return elems;
}

Float determinant(const Matrix3& mat) {
	return (
		mat(0,0) * (mat(1,1)*mat(2,2) - mat(1,2)*mat(2,1)) -
		mat(0,1) * (mat(1,0)*mat(2,2) - mat(1,2)*mat(2,0)) +
		mat(0,2) * (mat(1,0)*mat(2,1) - mat(1,1)*mat(2,0))
		);			
}

