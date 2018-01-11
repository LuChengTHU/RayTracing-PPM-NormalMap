#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include "color.h"
#include "vector3.h"
#include "material.h"
#include <iostream>
#include <sstream>
#include <string>
#include <algorithm>

extern const double INF;
extern const double EPS;
extern const double PI;
using namespace std;
class Primitive;

class Collider {
	Primitive* pri;
	
public:
	double dist;
	bool crash, front;
	Vector3 N, C, I;
	double u, v;
	
	Collider() {
		pri = NULL;
		crash = false;
	}
	~Collider() {}

	Primitive* GetPrimitive() { return pri; }
	void SetPrimitive(Primitive* _pri) { pri = _pri; }
};

class Primitive {
protected:
	int sample;
	Material* material;
	Primitive* next;

public:
	Primitive();
	Primitive( const Primitive& );
	~Primitive();
	
	void SetSample(int _sample) { sample = _sample; }
	int GetSample() { return sample; }
	void SetMaterial(Material* _material) { material = _material; }
	Material* GetMaterial() { return material; }
	Primitive* GetNext() { return next; }
	void SetNext( Primitive* primitive ) { next = primitive; }
	
	virtual int getName() { return 0; }
	virtual void PreTreatment() {}
	virtual void Input( std::string , std::stringstream& );
	virtual Collider Collide(Vector3 ray_O, Vector3 ray_V) = 0;
	virtual Color GetTexture(Vector3 C);
};

#endif
