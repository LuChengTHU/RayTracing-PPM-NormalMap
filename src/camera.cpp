#include"camera.h"
#include<cstdio>
#include<string>
#include<sstream>
#include<iostream>
#include<cstdlib>
#define ran() (double(rand() * (RAND_MAX + 1) + rand()) / ((RAND_MAX + 1) * (RAND_MAX + 1) - 1))

const int KD_MAX_THREADS = 10;
const int PM_MAX_THREADS = 10;
const int RT_MAX_THREADS = 10;
const std::string STD_ALGORITHM = "SPPM";
const int STD_DOF_SAMPLE = 64;
const double STD_APERTURE = 0;
const double STD_FOCAL_LEN = 1;
const double STD_LENS_WIDTH = 0.88;
const double STD_LENS_HEIGHT = 0.88;
const int STD_IMAGE_WIDTH = 420;
const int STD_IMAGE_HEIGHT = 420;
const int STD_SHADE_QUALITY = 4;
const int STD_DREFL_QUALITY = 20;
const int STD_MAX_HITPOINTS = 4000000;
const int STD_ITERATIONS = 5000;
const double STD_REDUCTION = 0.7;
const int STD_MAX_PHOTONS = 500000;
const int STD_EMIT_PHOTONS = 100000;
const int STD_SAMPLE_PHOTONS = 10;
const double STD_SAMPLE_DIST = 0.1;

Camera::Camera() {
	algorithm = STD_ALGORITHM;
	O = Vector3( 0 , 0 , 0 );
	N = Vector3( 0 , 1 , 0 );
	dofSample = STD_DOF_SAMPLE;
	aperture = STD_APERTURE;
	focalLen = STD_FOCAL_LEN;
	lens_W = STD_LENS_WIDTH;
	lens_H = STD_LENS_HEIGHT;
	W = STD_IMAGE_WIDTH;
	H = STD_IMAGE_HEIGHT;
	shade_quality = STD_SHADE_QUALITY;
	drefl_quality = STD_DREFL_QUALITY;
	max_hitpoints = STD_MAX_HITPOINTS;
	iterations = STD_ITERATIONS;
	reduction = STD_REDUCTION;
	max_photons = STD_MAX_PHOTONS;
	emit_photons = STD_EMIT_PHOTONS;
	sample_photons = STD_SAMPLE_PHOTONS;
	sample_dist = STD_SAMPLE_DIST;
	data = NULL;
}

Camera::~Camera() {
	if ( data == NULL ) {
		for ( int i = 0 ; i < H ; i++ )
			delete[] data[i];
		delete[] data;
	}
}

void Camera::Initialize() {
	N = N.GetUnitVector();
	Dx = N.GetAnVerticalVector();
	Dy = Dx.Cross(N);

	data = new Color*[H];
	for ( int i = 0 ; i < H ; i++ )
		data[i] = new Color[W];
}

Vector3 Camera::Emit( double i , double j ) {
	return N + Dy * lens_H * (i / (H - 1) - 0.5) + Dx * lens_W * (j / (W - 1) - 0.5);
}

void Camera::DofEmit(double i, double j, Vector3* dof_O, Vector3* dof_V) {
	Vector3 focalPoint = O + Emit(i, j) * focalLen;
	double x, y;
	do {
		x = ran() * 2 - 1;
		y = ran() * 2 - 1;
	} while (x * x + y * y > 1);
	*dof_O = O + Dx * aperture * x + Dy * aperture * y;
	*dof_V = (focalPoint - *dof_O).GetUnitVector();
}

void Camera::Input( std::string var , std::stringstream& fin ) {
	if ( var == "algorithm=" ) fin >> algorithm;
	if ( var == "O=" ) O.Input( fin );
	if ( var == "N=" ) N.Input( fin );
	if ( var == "dofSample=" ) fin >> dofSample;
	if ( var == "aperture=" ) fin >> aperture;
	if ( var == "focalLen=" ) fin >> focalLen;
	if ( var == "lens_W=" ) fin >> lens_W;
	if ( var == "lens_H=" ) fin >> lens_H;
	if ( var == "image_W=" ) fin >> W;
	if ( var == "image_H=" ) fin >> H;
	if ( var == "shade_quality=" ) fin >> shade_quality;
	if ( var == "drefl_quality=" ) fin >> drefl_quality;
	if ( var == "max_hitpoints=" ) fin >> max_hitpoints;
	if ( var == "iterations=" ) fin >> iterations;
	if ( var == "reduction=" ) fin >> reduction;
	if ( var == "max_photons=" ) fin >> max_photons;
	if ( var == "emit_photons=" ) fin >> emit_photons;
	if ( var == "sample_photons=" ) fin >> sample_photons;
	if ( var == "sample_dist=" ) fin >> sample_dist;
}

void Camera::Output( Bmp* bmp ) {
	bmp->Initialize( H , W );

	for ( int i = 0 ; i < H ; i++ )
		for ( int j = 0 ; j < W ; j++ )
			bmp->SetColor( i , j , data[i][j] );
}
