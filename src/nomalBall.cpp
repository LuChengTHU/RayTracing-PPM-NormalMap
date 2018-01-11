#include "nomalBall.h"

NomalBall::NomalBall():Primitive()
{
	De = Vector3(0, 0, 1);
	Dc = Vector3(0, 1, 0);
}

Vector3 NomalBall::getN(Vector3 C)
{
	Vector3 I = (C - O).GetUnitVector();
	double a = acos(-I.Dot(De));
	double b = acos(std::min(std::max(I.Dot(Dc) / sin(a), -1.0), 1.0));
	double u = a / PI, v = b / 2 / PI;
	if (I.Dot(Dc.Cross(De)) < 0) v = 1 - v;
	Color addN = ntexture->GetSmoothColor(u, v) * 2 - Color(1, 1, 1);
	//cout << "N: " << I.x<<" "<<I.y<<" "<<I.z << endl;
	//cout << "add: "<<addN.r<<' '<<addN.g<<' '<<addN.b << endl;
	Vector3 horizontal= C.y*C.y + C.x*C.x < EPS ? Vector3(0, 1, 0) : Vector3(C.y, -C.x, 0).GetUnitVector();
	Vector3 newN = I * addN.b + horizontal*addN.g + I.Cross(horizontal).GetUnitVector() * addN.r;
	//cout << "new: "<<newN.x<<' '<<newN.y<<' '<<newN.z << endl;
	//getchar();
	return newN.GetUnitVector();
}

void NomalBall::Input(std::string var, std::stringstream &fin)
{
	if (var == "O=") O.Input(fin);
	if (var == "R=") fin >> R;
	if (var == "De=") De.Input(fin);
	if (var == "Dc=") Dc.Input(fin);
	if (var == "NT=") {
		std::string file; fin >> file;
		ntexture = new Bmp;
		ntexture->Input(file);
	}
	Primitive::Input(var, fin);
}

Collider NomalBall::Collide(Vector3 ray_O, Vector3 ray_V)
{
	Collider collider;
	ray_V = ray_V.GetUnitVector();
	Vector3 P = ray_O - O;
	double b = -P.Dot(ray_V);
	double det = b * b - P.Module2() + R * R;

	if (det > EPS) {
		det = sqrt(det);
		double x1 = b - det, x2 = b + det;

		if (x2 < EPS) return collider;
		if (x1 > EPS) {
			collider.dist = x1;
			collider.front = true;
		}
		else {
			collider.dist = x2;
			collider.front = false;
		}
	}
	else
		return collider;

	collider.crash = true;
	collider.I = ray_V;
	collider.SetPrimitive(this);
	collider.C = ray_O + ray_V * collider.dist;
	collider.N = getN(collider.C);
	if (!collider.front) collider.N = -collider.N;
	return collider;
}

Color NomalBall::GetTexture(Vector3 C)
{
	Vector3 I = (C - O).GetUnitVector();
	double a = acos(-I.Dot(De));
	double b = acos(std::min(std::max(I.Dot(Dc) / sin(a), -1.0), 1.0));
	double u = a / PI, v = b / 2 / PI;
	if (I.Dot(Dc.Cross(De)) < 0) v = 1 - v;
	return material->texture->GetSmoothColor(u, v);
}

