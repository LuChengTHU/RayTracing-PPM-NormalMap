#ifndef PLANE_H
#define PLANE_H

#include "primitive.h"

class Plane : public Primitive {
	Vector3 N , Dx , Dy;
	double R;

public:
	Plane();
	~Plane() {}

	void Input( std::string , std::stringstream& );
	Collider Collide( Vector3 ray_O , Vector3 ray_V );
	Color GetTexture(Vector3 C);
};

#endif
