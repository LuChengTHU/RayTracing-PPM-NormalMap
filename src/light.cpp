#include"light.h"
#include<sstream>
#include<string>
#include<cmath>
#include<cstdlib>
#include<ctime>
#include<iostream>

Light::Light() {
	sample = rand();
	next = NULL;
}

void Light::Input( std::string var , std::stringstream& fin ) {
	if ( var == "color=" ) color.Input( fin );
}

Color Light::CalnIrradiance( Collider* collider , Vector3 V , int* hash ) {
	Primitive* pri = collider->GetPrimitive();
	Color ret = color * pri->GetMaterial()->BRDF(V, collider->N, -collider->I);
	ret /= V.Module2(); //rt adapts to ppm

	if (!ret.IsZeroColor() && hash != NULL)
		*hash = ( *hash + GetSample() ) % HASH_MOD;

	return ret;
}
