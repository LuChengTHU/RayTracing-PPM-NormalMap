#include "material.h"
#include <iostream>
#include <cmath>

Material::Material() {
	color = absor = Color();
	refl = refr = diff = spec = 0;
	rindex = 0;
	drefl = 0;
	texture = bump = NULL;
}

double Material::BRDF(Vector3 ray_R, Vector3 N, Vector3 ray_I) {
	double ret = 0;
	ray_R = ray_R.GetUnitVector();
	ray_I = ray_I.GetUnitVector();
	
	if (diff > EPS && ray_R.Dot(N) > EPS)
		ret += diff * ray_R.Dot(N);
	if (spec > EPS && ray_R.Dot(-ray_I.Reflect(N)) > EPS)
		ret += spec * pow(ray_R.Dot(-ray_I.Reflect(N)), 50);
	
	return ret;
}

void Material::Input( std::string var , std::stringstream& fin ) {
	if ( var == "color=" ) color.Input( fin );
	if ( var == "absor=" ) absor.Input( fin );
	if ( var == "refl=" ) fin >> refl;
	if ( var == "refr=" ) fin >> refr;
	if ( var == "diff=" ) fin >> diff;
	if ( var == "spec=" ) fin >> spec;
	if ( var == "drefl=" ) fin >> drefl;
	if ( var == "rindex=" ) fin >> rindex;
	if ( var == "texture=" ) {
		std::string file; fin >> file;
		texture = new Bmp;
		texture->Input( file );
	}
	if (var == "bump=") {
		std::string file; fin >> file;
		bump = new Bmp;
		bump->Input(file);
	}
}
