#include "plane.h"
#include <cmath>

Plane::Plane() : Primitive() {
	
}

void Plane::Input( std::string var , std::stringstream& fin ) {
	if ( var == "N=" ) N.Input( fin );
	if ( var == "R=" ) fin >> R;
	if ( var == "Dx=" ) Dx.Input( fin );
	if ( var == "Dy=" ) Dy.Input( fin );
	Primitive::Input( var , fin );
}

Collider Plane::Collide( Vector3 ray_O , Vector3 ray_V ) {
	Collider collider;
	ray_V = ray_V.GetUnitVector();
	N = N.GetUnitVector();
	double d = N.Dot( ray_V );
	if ( fabs( d ) < EPS ) return collider;
	double l = ( N * R - ray_O ).Dot( N ) / d;
	if ( l < EPS ) return collider;
	
	collider.crash = true;
	collider.I = ray_V;
	collider.SetPrimitive(this);
	collider.dist = l;
	collider.front = ( d < 0 );
	collider.C = ray_O + ray_V * collider.dist;
	collider.N = ( collider.front ) ? N : -N;
	return collider;
}

double get1(double x) {
	while (x > 1.0)x -= 1.0;
	while (x < 0)x += 1.0;
	return x;
}

Color Plane::GetTexture(Vector3 C) {
	double u = C.Dot( Dx ) / Dx.Module2();
	double v = C.Dot( Dy ) / Dy.Module2();
	u = get1(u);
	v = get1(v);
	return material->texture->GetSmoothColor( u , v );
}
