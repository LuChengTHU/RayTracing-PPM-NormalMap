#ifndef MATERIAL_H
#define MATERIAL_H

#include "color.h"
#include "vector3.h"
#include "bmp.h"

class Material {
public:
	Color color , absor;
	double refl , refr;
	double diff , spec;
	double rindex;
	double drefl;
	Bmp* texture;
	Bmp* bump;

	Material();
	~Material() {}

	void Input( std::string , std::stringstream& );
	double BRDF(Vector3 ray_R, Vector3 N, Vector3 ray_I);
};

#endif
