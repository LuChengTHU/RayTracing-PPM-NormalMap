#include "sphere.h"
#include <cmath>

Sphere::Sphere() : Primitive() {
	De = Vector3( 0 , 0 , 1 );
	Dc = Vector3( 0 , 1 , 0 );
}

void Sphere::Input( std::string var , std::stringstream& fin ) {
	if ( var == "O=" ) O.Input( fin );
	if ( var == "R=" ) fin >> R;
	if ( var == "De=" ) De.Input( fin );
	if ( var == "Dc=" ) Dc.Input( fin );
	Primitive::Input( var , fin );
}

Collider Sphere::Collide( Vector3 ray_O , Vector3 ray_V ) {
	Collider collider;
	ray_V = ray_V.GetUnitVector();
	Vector3 P = ray_O - O;
	double b = -P.Dot( ray_V );
	double det = b * b - P.Module2() + R * R;

	if ( det > EPS ) {
		det = sqrt( det );
		double x1 = b - det  , x2 = b + det;

		if ( x2 < EPS ) return collider;
		if ( x1 > EPS ) {
			collider.dist = x1;
			collider.front = true;
		} else {
			collider.dist = x2;
			collider.front = false;
		} 
	} else 
		return collider;

	collider.crash = true;
	collider.I = ray_V;
	collider.SetPrimitive(this);
	collider.C = ray_O + ray_V * collider.dist;
	collider.N = ( collider.C - O ).GetUnitVector();
	if (!collider.front) collider.N = -collider.N;
	return collider;
}

Color Sphere::GetTexture(Vector3 C) {
	Vector3 I = ( C - O ).GetUnitVector();
	double a = acos( -I.Dot( De ) );
	double b = acos( std::min( std::max( I.Dot( Dc ) / sin( a ) , -1.0 ) , 1.0 ) );
	double u = a / PI , v = b / 2 / PI;
	if ( I.Dot(Dc.Cross(De)) < 0 ) v = 1 - v;
	return material->texture->GetSmoothColor( u , v );
}
