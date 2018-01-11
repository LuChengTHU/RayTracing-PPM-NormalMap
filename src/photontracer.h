#ifndef PHOTONTRACER
#define PHOTONTRACER

#include"scene.h"
#include"photonmap.h"
#include"hitpointMap.h"

extern const int PM_MAX_THREADS;
extern const int MAX_PHOTONTRACING_DEP;

class Photontracer {
	Scene* scene;
	Photonmap* photonmap;
	HitpointMap* hitpointMap;
	
	int iteration;
	bool* completeThread;

	void PhotonTracing( Photon , int dep , bool refracted );
	bool PhotonDiffusion( Collider* , Photon , int dep , bool refracted , double* prob );
	bool PhotonReflection( Collider* , Photon , int dep , bool refracted , double* prob );
	bool PhotonRefraction( Collider* , Photon , int dep , bool refracted , double* prob );
	void Emitting(int threadID, int randID);

public:
	Photontracer();
	~Photontracer() {}

	void SetScene( Scene* _scene ) { scene = _scene; }
	void SetPhotonmap(Photonmap* _photonmap) { photonmap = _photonmap; }
	Photonmap* GetPhotonmap() { return photonmap; }
	void SetHitpointMap(HitpointMap* _hitpointMap) { hitpointMap = _hitpointMap; }
	Photonmap* CalnPhotonmap();
	void Run(int randIDBase = 0);
};

#endif
