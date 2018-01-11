#ifndef AREALIGHT_H
#define AREALIGHT_H

#include "light.h"

class AreaLight : public Light {
	Vector3 O , Dx , Dy;
public:
	AreaLight() : Light() {}
	~AreaLight() {}
	
	Vector3 GetO() { return O; }
	void Input( std::string , std::stringstream& );
	LightCollider Collide( Vector3 ray_O , Vector3 ray_V );
	Color GetIrradiance( Collider* collider , Primitive* primitive_head , int shade_quality , int* hash );
	Photon EmitPhoton();
};

#endif
