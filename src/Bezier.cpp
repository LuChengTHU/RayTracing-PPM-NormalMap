#include "Bezier.h"

#define ran() (double(rand() * (RAND_MAX + 1) + rand()) / ((RAND_MAX + 1) * (RAND_MAX + 1) - 1))


MyPoint Bezier::getBezierP(double t)
{
	MyPoint p = (P[6] - 6 * P[5] + 15 * P[4] - 20 * P[3] + 15 * P[2] - 6 * P[1] + P[0])*pow(t, 6) + (6 * P[5] - 30 * P[4] + 60 * P[3] - 60 * P[2] + 30 * P[1] - 6 * P[0])*pow(t, 5) + (15 * P[4] - 60 * P[3] + 90 * P[2] - 60 * P[1] + 15 * P[0])*pow(t, 4) + (20 * P[3] - 60 * P[2] + 60 * P[1] - 20 * P[0])*pow(t, 3) + (15 * P[2] - 30 * P[1] + 15 * P[0])*pow(t, 2) + (6 * P[1] - 6 * P[0])*t + P[0];

	return p;
}

MyPoint Bezier::derivativeP(double t)
{
	MyPoint dp = (P[6] - 6 * P[5] + 15 * P[4] - 20 * P[3] + 15 * P[2] - 6 * P[1] + P[0]) * 6 * pow(t, 5) + (6 * P[5] - 30 * P[4] + 60 * P[3] - 60 * P[2] + 30 * P[1] - 6 * P[0]) * 5 * pow(t, 4) + (15 * P[4] - 60 * P[3] + 90 * P[2] - 60 * P[1] + 15 * P[0]) * 4 * pow(t, 3) + (20 * P[3] - 60 * P[2] + 60 * P[1] - 20 * P[0]) * 3 * pow(t, 2) + (15 * P[2] - 30 * P[1] + 15 * P[0]) * 2 * pow(t, 1) + (6 * P[1] - 6 * P[0]);
	return dp;
}

Vector3d Bezier::computeF(Vector3d x, Vector3 R0, Vector3 Rd)
{
	Vector3 L = R0 + Rd*x(0);
	Vector3d F = Vector3d(L.x, L.y, L.z) - computeS(x(1), x(2));
	return F;
}

Vector3d Bezier::computeS(double t, double theta)
{
	Vector3d S;
	double H = bzmax - bzmin;
	MyPoint p = getBezierP(t);
	S << O.x + (p.x + R)*cos(theta), O.y + (p.x + R)*sin(theta), O.z + H - p.y;
	return S;
}

Matrix3d Bezier::dF(Vector3d x, Vector3 R0, Vector3 Rd)
{
	Matrix3d df;
	MyPoint dp = derivativeP(x(1));
	MyPoint p = getBezierP(x(1));

	//cout << "Line_t:" << x(0) << endl;
	//cout << "Curve_t" << x(1) << endl;
	//cout << "theta" << x(2) << endl;
	//cout << "R0:" << R0.x << " , " << R0.y << " , " << R0.z << endl;
	//cout << "Rd:" << Rd.x << " , " << Rd.y << " , " << Rd.z << endl;
	//cout << "p:" << p.x << ',' << p.y << endl;
	//getchar();


	df << Rd.x, -dp.x*cos(x(2)), (R + p.x)*sin(x(2)),
		Rd.y, -dp.x*sin(x(2)), -(R + p.x)*cos(x(2)),
		Rd.z, dp.y, 0;
	//cout << (R + p.x)*sin(x(2)) << endl;
	//cout << df << endl;
	return df;
}

double Bezier::intersectBox(Vector3 R0, Vector3 Rd)
{
	Rd = Rd.GetUnitVector();
	double txmin, txmax, tymin, tymax, tzmin, tzmax;
	if (Rd.x < EPS) {
		txmin = -INF;
		txmax = INF;
	}
	else {
		txmin = (bxmin - R0.x) / Rd.x;
		txmax = (bxmax - R0.x) / Rd.x;
	}
	if (Rd.y < EPS) {
		tymin = -INF;
		tymax = INF;
	}
	else {
		tymin = (bymin - R0.y) / Rd.y;
		tymax = (bymax - R0.y) / Rd.y;
	}
	if (Rd.z < EPS) {
		tzmin = -INF;
		tzmax = INF;
	}
	else {
		tzmin = (bzmin - R0.z) / Rd.z;
		tzmax = (bzmax - R0.z) / Rd.z;
	}
	if (txmin > txmax) swap(txmin, txmax);
	if (tymin > tymax) swap(tymin, tymax);
	if (tzmin > tzmax) swap(tzmin, tzmax);
	double tmin = max(max(txmin, tymin), tzmin);
	double tmax = min(min(txmax, tymax), tzmax);
	if (bxmin <= R0.x + 0.001 && R0.x <= bxmax + 0.001 && bymin <= R0.y + 0.001 &&R0.y <= bymax + 0.001 && bzmin <= R0.z + 0.001 && R0.z <= bzmax + 0.001)return tmax > EPS ? tmax : 0;
	else if (tmin >= tmax) return 0;
	else return tmin > EPS ? tmin : 0;
}

double Bezier::intersectBottom(Vector3 R0, Vector3 Rd)
{

	Rd = Rd.GetUnitVector();
	double d = N.Dot(Rd);
	if (fabs(d) < EPS) return 0;
	double l = (N * R - R0).Dot(N) / d;
	if (l < EPS) return 0;
	Vector3 OP = R0 + Rd*l - O;
	if (R*R - OP.Module2()< EPS)return 0;
	return l;
}

void Bezier::getMesh()
{
	FILE* fp = fopen("Bezier.obj", "w");
	if (fp == NULL)cout << "NULL";
	std::vector<Vector3d> points(0);          // 存储点
	std::vector<Vector4i> meshes(0);            // 存储面
	float du = 0.01f, dv = 0.01f;             // 自定义密度
	int nu = 1 / du + 1, nv = 1 / dv + 1;   // 密度应该“除”得整

	int pid = 0;                            // 点序号
	for (float u = 0.0f, i = 0; u <= 1.0f; u += du, i++) {
		for (float v = 0.0f, j = 0; v <= 1.0f; v += dv, j++) {
			Vector3d S = computeS(u, 2 * PI*v);
			fprintf(fp, "v %lf %lf %lf\n", S(0), S(1), S(2));
			// 写递归或者DP搞定函数P的计算
			pid++;                              // OBJ格式网格序号从1开始
			if (i != 0 && j != 0) {
				meshes.push_back(Vector4i(pid - nv - 1, pid - nv, pid, pid - 1));

			}
			cout << i << ' ' << j << endl;
		}
	}

	for (int i = 0; i < meshes.size(); i++) {
		fprintf(fp, "f %d %d %d %d\n", meshes[i][0], meshes[i][1], meshes[i][2], meshes[i][3]);
	}
}

Bezier::Bezier()
{
	N = Vector3(0, 0, -1.0);


	//O = Vector3(0, 0, 0);
	//R = 5.0;
	//ratio = 1.0;
	//P[0] = MyPoint(0.4 * ratio, 0.0);
	//P[1] = MyPoint(2.4 * ratio, 13.0 * ratio);
	//P[2] = MyPoint(-4.0 * ratio, 15.0 * ratio);
	//P[3] = MyPoint(-6.0 * ratio, 21.0 * ratio);
	//P[4] = MyPoint(-1.6 * ratio, 22.5 * ratio);
	//P[5] = MyPoint(0.0 * ratio, 38.0 * ratio);

	//bxmin = O.x - R - 2.4 * ratio;
	//bxmax = O.x + R + 2.4 * ratio;
	//bymin = O.y - R - 2.4 * ratio;
	//bymax = O.y + R + 2.4 * ratio;
	//bzmin = O.z;
	//bzmax = O.z + 38.0 * ratio;
}


void Bezier::Input(std::string var, std::stringstream &fin)
{
	if (var == "O=") O.Input(fin);
	if (var == "R=") fin >> R;
	if (var == "ratio=") fin >> ratio;


	//P[0] = MyPoint(2.9 * ratio, 0.0);
	//P[1] = MyPoint(3.9 * ratio, 12.0 * ratio);
	//P[2] = MyPoint(1.3 * ratio, 22.0 * ratio);
	//P[3] = MyPoint(-0.9 * ratio, 40.0 * ratio);
	//P[4] = MyPoint(-3.0 * ratio, 45 * ratio);
	//P[5] = MyPoint(-6.0 * ratio, 48 * ratio);
	//P[6] = MyPoint(0.0 * ratio, 52.0 * ratio);

	P[0] = MyPoint(1.5 * ratio, 0 * ratio);						//204,0 
	P[1] = MyPoint(2.9 * ratio, 13.0 * ratio);				//224,130
	P[2] = MyPoint(-3.5 * ratio, 15.0 * ratio);				//160,150
	P[3] = MyPoint(-5.5 * ratio, 21.0 * ratio);				//140,210
	P[4] = MyPoint(-1.1 * ratio, 22.5 * ratio);				//184,225
	P[5] = MyPoint(-5.5 * ratio, 25 * ratio);				//140,250
	P[6] = MyPoint(0 * ratio, 38.0 * ratio);				//200,380

	double dx = 2.9 + 5.5;
	double dz = 38.0;
	bxmin = O.x - R - dx * ratio;
	bxmax = O.x + R + dx * ratio;
	bymin = O.y - R - dx * ratio;
	bymax = O.y + R + dx * ratio;
	bzmin = O.z;
	bzmax = O.z + dz * ratio;

	Primitive::Input(var, fin);
}

Collider Bezier::Collide(Vector3 R0, Vector3 Rd)
{
	bool interCurve = true;
	Collider collider;
	Rd = Rd.GetUnitVector();
	double t = intersectBox(R0, Rd);
	if (t <= EPS) return collider;

	/*
	Vector3 l = Rd*t + R0 - O;
	double tc;
	double dh = bzmax - bzmin - l.z;
	if (dh < P[1].y)tc = 0.2 * (dh - P[0].y) / (P[1].y - P[0].y);
	else if (dh < P[2].y)tc = 0.2 * (dh - P[1].y) / (P[2].y - P[1].y) + 0.2;
	else if (dh < P[3].y)tc = 0.2 * (dh - P[2].y) / (P[3].y - P[2].y) + 0.4;
	else if (dh < P[4].y)tc = 0.2 * (dh - P[3].y) / (P[4].y - P[3].y) + 0.6;
	else tc = 0.2 * (dh - P[4].y) / (P[5].y - P[4].y) + 0.8;


	double theta = acos((l-Vector3(0,0, (Rd*t + R0).z)).GetUnitVector().Dot(Vector3(1,0,0)));
	if ((l - Vector3(0, 0, (Rd*t + R0).z)).GetUnitVector().Dot(Vector3(0, 1, 0)) <= 0)theta = 2 * PI - theta;
	if(fabs(theta-PI)<PI/8)theta = PI - asin((l - Vector3(0, 0, (Rd*t + R0).z)).GetUnitVector().Dot(Vector3(0, 1, 0)));
	if(fabs(theta)<PI/8)theta = asin((l - Vector3(0, 0, (Rd*t + R0).z)).GetUnitVector().Dot(Vector3(0, 1, 0)));
	*/
	Vector3d x(INF, -1, -1);
	for (int k = 0; k < 10; k++)
	{
		double tc = ran();
		double theta = ran() * 2 * PI;
		double tl = t;
		Vector3d y(tl, ran(), ran() * 2 * PI);
		double dt;
		Vector3d F;
		for (int i = 0;; i++) {
			double pret = y(0);
			Matrix3d df = dF(y, R0, Rd);
			F = computeF(y, R0, Rd);
			Vector3d f = df.inverse()*F;
			//if (fabs(f(1) - x(1)) >= 1) {
			//	f = x + (f - x) / (fabs(f(1) - x(1)) + 1.0);
			//}
			//f = f / 2;
			//if (i == 20)f = f / 2;
			y = y - f;
			dt = y(0) - pret;
			dt = dt > 0 ? dt : -dt;
			if (i > 20)break;
		}
		if (dt > 0.001 || y(0) < EPS || y(1) < 0 || y(1) > 1 || fabs(F(0))+fabs(F(1))+fabs(F(2)) > 0.1)continue;
		else if (x(0) > y(0)) x = y;
	}

	//cout << x(0) << " , " << x(1) << endl;
	//getchar();
	if (x(0) < EPS || x(1) < 0 || x(1) > 1)interCurve = false;
	/*
	double bottom_t = intersectBottom(R0, Rd);
	if (!interCurve&&bottom_t <= EPS) {
	return collider;
	}
	if ((!interCurve && bottom_t > EPS) || (x(0) > bottom_t && bottom_t > EPS)) {
	collider.crash = true;
	collider.I = Rd;
	collider.SetPrimitive(this);
	collider.dist = bottom_t;
	double d = N.Dot(Rd);
	collider.front = (d < 0);
	collider.C = R0 + Rd * collider.dist;
	collider.N = (collider.front) ? N : -N;
	return collider;
	}
	*/
	if (!interCurve)return collider;


	Matrix3d df = dF(x, R0, Rd);

	Vector3d dfdt(df(0, 1), df(1, 1), df(2, 1));
	Vector3d dfdtheta(df(0, 2), df(1, 2), df(2, 2));
	Vector3d Nvec = dfdt.cross(dfdtheta);
	Vector3 N(Nvec(0), Nvec(1), Nvec(2));
	N = N.GetUnitVector();
	double d = N.Dot(Rd);
	collider.crash = true;
	collider.I = Rd;
	collider.SetPrimitive(this);
	collider.dist = x(0);
	collider.front = (d < 0);
	collider.C = R0 + Rd * collider.dist;
	collider.N = (collider.front) ? N : -N;
	collider.u = x(1);
	collider.v = x(2);
	return collider;
}

Color Bezier::GetTexture(Vector3 C)
{
	//if (N.Dot(C - O) < EPS) {
	//	double u = (C - O).Module();
	//	double a = (C - O).Dot(Vector3(1, 0, 0));
	//	double v = acos(a) / (2 * PI);
	//	return material->texture->GetSmoothColor(u, v);
	//}
	double t = C.y;
	while (t > 2 * PI)t = t - 2 * PI;
	while (t < 0)t = t + 2 * PI;
	return material->texture->GetSmoothColor(C.x, t / (2 * PI));
}






/*

class Bezier3 : public Object
{
public:
Matrix3d Px, Py, Pz, J;
double bxmin, bxmax, bymin, bymax, bzmin, bzmax;
double u, v, t;
Bezier3(
Vec p00, Vec p01, Vec p02, Vec p10, Vec p11, Vec p12, Vec p20, Vec p21, Vec p22,
Vec e_, Vec c_, Refl_t refl_, Texture *texture_ = NULL) :
Object(e_, c_, refl_, texture_)
{
Px << p00.x, p01.x, p02.x, p10.x, p11.x, p12.x, p20.x, p21.x, p22.x;
Py << p00.y, p01.y, p02.y, p10.y, p11.y, p12.y, p20.y, p21.y, p22.y;
Pz << p00.z, p01.z, p02.z, p10.z, p11.z, p12.z, p20.z, p21.z, p22.z;
bxmin = bymin = bzmin = 1e10;
bxmax = bymax = bzmax = -1e10;
for (int i = 0; i < 3; i++)
for (int j = 0; j < 3; j++)
{
bxmin = min(bxmin, Px(i, j));
bxmax = max(bxmax, Px(i, j));
bymin = min(bymin, Py(i, j));
bymax = max(bymax, Py(i, j));
bzmin = min(bzmin, Pz(i, j));
bzmax = max(bzmax, Pz(i, j));
}
bxmin--, bymin--, bzmin--;
bxmax++, bymax++, bzmax++;
}
double intersect(const Ray &r)
{
double ox = r.o.x, oy = r.o.y, oz = r.o.z;
double dx = r.d.x, dy = r.d.y, dz = r.d.z;
double tprev = 0, dt = 0;
if ((t = intersectBox(r)) == 0) return 0;
for (int i = 0; ; i++) {
Vector3d X(u, v, t);
Vector3d FF(F(Px, ox, dx, u, v, t), F(Py, oy, dy, u, v, t), F(Pz, oz, dz, u, v, t));
J << dFdu(Px, dx, u, v, t), dFdv(Px, dx, u, v, t), dFdt(Px, dx, u, v, t),
dFdu(Py, dy, u, v, t), dFdv(Py, dy, u, v, t), dFdt(Py, dy, u, v, t),
dFdu(Pz, dz, u, v, t), dFdv(Pz, dz, u, v, t), dFdt(Pz, dz, u, v, t);
X = X - J.inverse() * FF;
u = X(0, 0), v = X(1, 0);
tprev = t;
t = X(2, 0);
dt = tprev - t;
dt = dt < 0 ? -dt : dt;
if (dt < 1e-8)
break;
if (i > 20)
break;
}
double eps = 1e-4;
return (t > eps && u >= 0 && u <= 1 && v >= 0 && v <= 1) ? t : 0;
}

void getSurfaceData(Vec &phit, Vec &nhit, Vec &color) const
{
Vec dfdu(J(0, 0), J(1, 0), J(2, 0));
Vec dfdv(J(0, 1), J(1, 1), J(2, 1));
nhit = (dfdu % dfdv).norm();
if (texture)
{
color = texture->map(u, v);
// color = c;
}
else
color = c;
}

private:
// Vec texture(double u, double v) const
// {
//     return ((int)(8 * u) + (int)(8 * v)) % 2 == 0? Vec() : Vec(1, 1, 1);
// }
double intersectBox(const Ray &r)
{
double txmin = (bxmin - r.o.x) / r.d.x, txmax = (bxmax - r.o.x) / r.d.x;
double tymin = (bymin - r.o.y) / r.d.y, tymax = (bymax - r.o.y) / r.d.y;
double tzmin = (bzmin - r.o.z) / r.d.z, tzmax = (bzmax - r.o.z) / r.d.z;
if (txmin > txmax) swap(txmin, txmax);
if (tymin > tymax) swap(tymin, tymax);
if (tzmin > tzmax) swap(tzmin, tzmax);
double tmin = max(max(txmin, tymin), tzmin);
double tmax = min(min(txmax, tymax), tzmax);
if (tmin >= tmax) return 0;
double eps = 1e-4;
return tmin > 1e-4 ? tmin : 0;
}
double dFdu(Matrix3d &p, double d, double u, double v, double t)
{
Vector3d U(2 * u, 2 - 4 * u, 2 * u - 2);
Vector3d V(v * v, 2 * v * (1 - v), (1 - v) * (1 - v));
return U.transpose() * p * V;
}

double dFdv(Matrix3d &p, double d, double u, double v, double t)
{
Vector3d U(u * u, 2 * u * (1 - u), (1 - u) * (1 - u));
Vector3d V(2 * v, 2 - 4 * v, 2 * v - 2);
return U.transpose() * p * V;
}

double dFdt(Matrix3d &p, double d, double u, double v, double t)
{
return -d;
}

double F(Matrix3d &p, double o, double d, double u, double v, double t)
{
Vector3d U(u * u, 2 * u * (1 - u), (1 - u) * (1 - u));
Vector3d V(v * v, 2 * v * (1 - v), (1 - v) * (1 - v));
return U.transpose() * p * V - o - t * d;
}
};
*/