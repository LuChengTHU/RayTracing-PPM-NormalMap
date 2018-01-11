#pragma once
#ifndef NOMALPLANE_H
#define NOMALPLANE_H

#include "primitive.h"

class NomalPlane : public Primitive {
	//ƽ�淨����
	Vector3 N, Dx, Dy;
	//ԭ�㵽ƽ��ľ���
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
