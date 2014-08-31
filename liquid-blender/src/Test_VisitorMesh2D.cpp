#include "stdafx.h"
#include "Utility.hpp"
#include "gtest/gtest.h"

using namespace std;

TEST(VisitorMesh2DTest, Visit) {
	Mesh2D mesh;
	
	mesh.load("tests/dropPool0/dropPool0_0000.mesh");
	VisitorMesh2D visitor(mesh);

	int correctOrder[] = {
		0 ,
		13 ,
		18 ,
		20 ,
		22 ,
		24 ,
		26 ,
		28 ,
		30 ,
		32 ,
		34 ,
		36 ,
		38 ,
		40 ,
		42 ,
		54 ,
		58 ,
		62 ,
		66 ,
		70 ,
		74 ,
		78 ,
		82 ,
		92 ,
		94 ,
		96 ,
		98 ,
		100 ,
		102 ,
		104 ,
		106 ,
		108 ,
		110 ,
		112 ,
		114 ,
		116 ,
		118 ,
		120 ,
		122 ,
		124 ,
		126 ,
		128 ,
		130 ,
		132 ,
		134 ,
		136 ,
		138 ,
		140 ,
		142 ,
		144 ,
		146 ,
		148 ,
		150 ,
		152 ,
		154 ,
		156 ,
		158 ,
		160 ,
		162 ,
		164 ,
		166 ,
		168 ,
		170 ,
		172 ,
		174 ,
		176 ,
		178 ,
		180 ,
		182 ,
		184 ,
		186 ,
		188 ,
		190 ,
		192 ,
		195 ,
		196 ,
		197 ,
		198 ,
		199 ,
		1 ,
		191 ,
		194 ,
		81 ,
		84 ,
		52 ,
		193 ,
		45 ,
		80 ,
		189 ,
		187 ,
		185 ,
		183 ,
		181 ,
		179 ,
		177 ,
		175 ,
		173 ,
		171 ,
		169 ,
		167 ,
		165 ,
		163 ,
		161 ,
		159 ,
		157 ,
		155 ,
		153 ,
		151 ,
		149 ,
		147 ,
		145 ,
		143 ,
		141 ,
		139 ,
		137 ,
		135 ,
		133 ,
		131 ,
		129 ,
		127 ,
		125 ,
		123 ,
		121 ,
		119 ,
		117 ,
		115 ,
		113 ,
		111 ,
		109 ,
		107 ,
		105 ,
		103 ,
		101 ,
		99 ,
		97 ,
		95 ,
		93 ,
		83 ,
		79 ,
		75 ,
		71 ,
		67 ,
		63 ,
		59 ,
		55 ,
		43 ,
		41 ,
		39 ,
		37 ,
		35 ,
		33 ,
		31 ,
		29 ,
		27 ,
		25 ,
		23 ,
		21 ,
		19 ,
		17 ,
		50 ,
		200 ,
		15 ,
		12 ,
		11 ,
		10 ,
		9 ,
		8 ,
		7 ,
		6 ,
		5 ,
		4 ,
		3 ,
		2 ,
		14 ,
		53 ,
		60 ,
		64 ,
		68 ,
		72 ,
		76 ,
		85 ,
		56 ,
		86 ,
		87 ,
		88 ,
		89 ,
		90 ,
		91 ,
		77 ,
		73 ,
		69 ,
		65 ,
		61 ,
		57 ,
		51 ,
		44 ,
		16 ,
		49 ,
		48 ,
		47 ,
		46
	};
	
	int curVertex;
	int count = 0;
	while (true) {		
		curVertex = visitor.next();		
		if (curVertex == -1) break;		
		EXPECT_EQ(correctOrder[count], curVertex);
		count++;
	};

	EXPECT_EQ(mesh.vertices.size(), count);

}