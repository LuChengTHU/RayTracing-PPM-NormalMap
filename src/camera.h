#ifndef CAMERA_H
#define CAMERA_H

#include"vector3.h"
#include"color.h"
#include"bmp.h"
#include<string>
#include<sstream>

extern const int KD_MAX_THREADS; //the maxinum threads of Building KDTree
extern const int RT_MAX_THREADS; //the maxinum threads of RayTracing
extern const int PM_MAX_THREADS; //the maxinum threads of PhotonsMapping
extern const std::string STD_ALGORITHM;
extern const int STD_DOF_SAMPLE;
extern const double STD_APERTURE;
extern const double STD_FOCAL_LEN;
extern const double STD_LENS_WIDTH; //the width of lens in the scene
extern const double STD_LENS_HEIGHT;
extern const int STD_IMAGE_WIDTH;
extern const int STD_IMAGE_HEIGTH;
extern const int STD_SHADE_QUALITY; //caln shade :: how many times will be run (*16)
extern const int STD_DREFL_QUALITY; //caln drefl :: how many times will be run (*16)
extern const int STD_MAX_HITPOINTS;
extern const int STD_ITERATIONS;
extern const double STD_REDUCTION;
extern const int STD_MAX_PHOTONS;
extern const int STD_EMIT_PHOTONS;
extern const int STD_SAMPLE_PHOTONS;
extern const double STD_SAMPLE_DIST;

class Camera {
	std::string algorithm;
	Vector3 O , N , Dx , Dy;
	double dofSample, aperture, focalLen;
	double lens_W , lens_H;
	int W , H;
	Color** data;
	double shade_quality;
	double drefl_quality;
	int max_hitpoints;
	int iterations;
	double reduction;
	int max_photons;
	int emit_photons;
	int sample_photons;
	double sample_dist;

public:
	Camera();
	~Camera();
	
	std::string GetAlgorithm() { return algorithm; }
	Vector3 GetO() { return O; }
	Vector3 GetN() { return N; }
	int GetDofSample() { return dofSample; }
	double GetAperture() { return aperture; }
	double GetFocalLen() { return focalLen; }
	int GetW() { return W; }
	int GetH() { return H; }
	Color GetColor(int i, int j) { return data[i][j]; }
	void SetColor( int i , int j , Color color ) { data[i][j] = color; }
	double GetShadeQuality() { return shade_quality; }
	double GetDreflQuality() { return drefl_quality; }
	int GetMaxHitpoints() { return max_hitpoints; }
	int GetIterations() { return iterations; }
	double GetReduction() { return reduction; }
	int GetMaxPhotons() { return max_photons; }
	int GetEmitPhotons() { return emit_photons; }
	int GetSamplePhotons() { return sample_photons; }
	double GetSampleDist() { return sample_dist; }

	Vector3 Emit( double i , double j );
	void DofEmit(double i, double j, Vector3* dof_O, Vector3* dof_V);
	void Initialize();
	void Input( std::string var , std::stringstream& fin );
	void Output( Bmp* );
};

#endif
