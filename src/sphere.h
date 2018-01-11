#ifndef SPHERE_H
#define SPHERE_H

#include "primitive.h"

class Sphere : public Primitive {
	Vector3 O , De , Dc;
	double R;

public:
	Sphere();
	~Sphere() {}

	void Input( std::string , std::stringstream& );
	Collider Collide( Vector3 ray_O , Vector3 ray_V );
	Color GetTexture(Vector3 C);
};

#endif