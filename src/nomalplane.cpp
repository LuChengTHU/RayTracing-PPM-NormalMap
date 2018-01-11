#include "nomalplane.h"

NomalPlane::NomalPlane() : Primitive() {}

double NomalPlane::get1(double x) {
	while (x > 1.0)x -= 1.0;
	while (x < 0)x += 1.0;
	return x;
}

Vector3 NomalPlane::getN(Vector3 C)
{
	double u = C.Dot(Dx) / Dx.Module2();
	double v = C.Dot(Dy) / Dy.Module2();
	u = get1(u);
	v = get1(v);
	Color addN = ntexture->GetSmoothColor(u, v) * 2 - Color(1,1,1);
	Vector3 horizontal;
	if (N.x == 0 && N.y == 0)horizontal = Vector3(1, 0, 0);
	else horizontal = Vector3(N.x, N.y, 0);
	Vector3 newN = N * addN.b +   horizontal * addN.r + N.Cross(horizontal).GetUnitVector() * addN.g;
	return newN.GetUnitVector();
}

void NomalPlane::Input(std::string var, std::stringstream &fin)
{
	if (var == "N=") N.Input(fin);
	N = N.GetUnitVector();
	if (var == "R=") fin >> R;
	if (var == "Dx=") Dx.Input(fin);
	if (var == "Dy=") Dy.Input(fin);
	if (var == "NT=")
	{
		std::string file; fin >> file;
		ntexture = new Bmp;
		ntexture->Input(file);
	}
	Primitive::Input(var, fin);
}

Collider NomalPlane::Collide(Vector3 R0, Vector3 Rd)
{
	Collider collider;
	Rd = Rd.GetUnitVector();
	N = N.GetUnitVector();
	double d = N.Dot(Rd);
	if (fabs(d) < EPS) return collider;
	double l = (N * R - R0).Dot(N) / d;
	if (l < EPS) return collider;

	collider.crash = true;
	collider.I = Rd;
	collider.SetPrimitive(this);
	collider.dist = l;
	collider.front = (d < 0);
	collider.C = R0 + Rd * collider.dist;
	collider.N = (collider.front) ? getN(collider.C) : -getN(collider.C);
	return collider;
}

Color NomalPlane::GetTexture(Vector3 C) {
	double u = C.Dot(Dx) / Dx.Module2();
	double v = C.Dot(Dy) / Dy.Module2();
	u = get1(u);
	v = get1(v);
	return material->texture->GetSmoothColor(u, v);
}
