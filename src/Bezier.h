#pragma once

#ifndef BEZIER_H
#define BEZIER_H

#include "primitive.h"
#include <ctime>


#include <vector>
#include <stdio.h>
#include "eigen\Eigen\Dense"
using namespace Eigen;

#define NN 7

struct MyPoint {
	double x;
	double y;
	MyPoint() {}
	MyPoint(double _x, double _y) :x(_x), y(_y) {}
	friend MyPoint operator * (const MyPoint &p, const double t) { return MyPoint(p.x*t, p.y*t); }
	friend MyPoint operator * (const double t, const MyPoint &p) { return MyPoint(p.x*t, p.y*t); }
	friend MyPoint operator * (const MyPoint &p, const int t) { return MyPoint(p.x*t, p.y*t); }
	friend MyPoint operator + (const MyPoint &p1, const MyPoint &p2) { return MyPoint(p1.x + p2.x, p1.y + p2.y); }
	friend MyPoint operator - (const MyPoint &p1, const MyPoint &p2) { return MyPoint(p1.x - p2.x, p1.y - p2.y); }
};


class Bezier : public Primitive {

	Vector3 O, N;//底面中心,底面法向量
	double R;	//曲线与O的偏移量
	double ratio;

	MyPoint P[NN];
	double bxmin, bxmax, bymin, bymax, bzmin, bzmax;
	MyPoint getBezierP(double t);
	MyPoint derivativeP(double t);
	Vector3d computeF(Vector3d x, Vector3 R0, Vector3 Rd);
	Vector3d computeS(double t, double theta);
	Matrix3d dF(Vector3d x, Vector3 R0, Vector3 Rd);
	double intersectBox(Vector3 R0, Vector3 Rd);
	double intersectBottom(Vector3 R0, Vector3 Rd);

public:
	Bezier();
	~Bezier() {}

	int getName() { return 1; }
	void getMesh();
	void Input(std::string, std::stringstream&);
	Collider Collide(Vector3 R0, Vector3 Rd);
	Color GetTexture(Vector3 C);

};

#endif