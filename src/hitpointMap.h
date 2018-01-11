#ifndef HITPOINTMAP_H
#define HITPOINTMAP_H

#include "vector3.h"
#include "color.h"
#include "primitive.h"
#include "photonmap.h"

extern const double PI;
extern const double INF;

class Collider;

class Hitpoint {
public:
	Vector3 pos, dir, N;
	Primitive* primitive;
	int rc;
	Color weight;
	int plane;
	float R2, maxR2;
	float num, deltaNum;
	Color color;
	
	Hitpoint();
	void CalnIrradiance(Photon*);
};

class HitpointMap {
	int maxHitpoints, storedHitpoints;
	Hitpoint* hitpoints;
	Vector3 boxMin, boxMax;
	double reduction;

	void BalanceSegment(Hitpoint* horg, int index, int st, int en);
	void MedianSplit(Hitpoint* horg, int st, int en, int med, int axis);
	void MaintainHitpointsMaxR2();
	void LocatePhoton(Photon* photon, int index); //called by index = 1

public:
	HitpointMap(int size);
	~HitpointMap();
	
	void SetReduction(double _reduction) { reduction = _reduction; }
	int GetStoredHitpoints() { return storedHitpoints; }
	Hitpoint* GetHitpoints() { return hitpoints; }
	void Store(Hitpoint);
	void MaintainHitpoints();
	void Balance();
	void InsertPhoton(Photon);
};

#endif
