#pragma once
#ifndef NOMALBALL_H
#define NOMALBALL_H

#include "primitive.h"

class NomalBall : public Primitive {
	Vector3 O, De, Dc;
	double R;
	Bmp* ntexture;
public:
	NomalBall();
	~NomalBall() {}
	
	Vector3 getN(Vector3 C);
	void Input(std::string, std::stringstream&);
	Collider Collide(Vector3 ray_O, Vector3 ray_V);
	Color GetTexture(Vector3 C);
};

#endif