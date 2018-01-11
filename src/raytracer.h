#ifndef RAYTRACER_H
#define RAYTRACER_H

#include"bmp.h"
#include"scene.h"
#include"photontracer.h"
#include<string>

extern const int RT_MAX_THREADS;
extern const int MAX_DREFL_DEP;
extern const int MAX_RAYTRACING_DEP;
extern const int HASH_FAC;
extern const int HASH_MOD;

class Raytracer {
	std::string input , output;
	Scene* scene;
	Camera* camera;
	Photonmap* photonmap;
	HitpointMap* hitpointMap;
	
	int H, W;
	int** sample;
	bool* completeThread;

	Color CalnDiffusion( Collider* collider , int* hash, int rc, Color weight);
	Color CalnReflection( Collider* collider , Vector3 ray_V , int dep , bool refracted , int* hash, int rc, Color weight);
	Color CalnRefraction( Collider* collider , Vector3 ray_V , int dep , bool refracted , int* hash, int rc, Color weight);
	Color RayTracing(Vector3 ray_O, Vector3 ray_V, int dep, bool refracted, int* hash, int rc, Color weight);	
	void Sampling(int threadID, int randID);
	void Resampling(int threadID, int randID);
	void MultiThreadSampling(int randIDBase = 0);
	void MultiThreadResampling(int randIDBase = 0);
	void ProgressivePhotonMapping(int SPPMIter = 0);
	void GenerateImage(std::string file);

public:
	Raytracer();
	~Raytracer();
	
	void SetInput( std::string file ) { input = file; }
	void SetOutput( std::string file ) { output = file; }
	void Run();
};

#endif
