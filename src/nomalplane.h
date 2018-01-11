#pragma once
#ifndef NOMALPLANE_H
#define NOMALPLANE_H

#include "primitive.h"

class NomalPlane : public Primitive {
	//平面法向量
	Vector3 N, Dx, Dy;
	//原点到平面的距离
	double R;
	Bmp* ntexture;
	double get1(double x);
public:
	NomalPlane();
	~NomalPlane() {}

	Vector3 getN(Vector3 C);
	void Input(std::string, std::stringstream&);
	Collider Collide(Vector3 R0, Vector3 Rd);
	Color GetTexture(Vector3 C);

};

#endif
