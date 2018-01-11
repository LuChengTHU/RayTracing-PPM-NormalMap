#ifndef SCENE_H
#define SCENE_H

#include "color.h"
#include "vector3.h"
#include "plane.h"
#include "sphere.h"
#include "Bezier.h"
#include "nomalBall.h"
#include "nomalplane.h"
#include "areaLight.h"
#include "camera.h"
#include <string>
#include <fstream>
#include <sstream>

class Scene {
	Primitive* primitive_head;
	Light* light_head;
	Camera* camera;
	Color background_color;
	
	void BackgroundInput(std::string var, std::stringstream& fin);

public:
	Scene();
	~Scene();
	
	Primitive* GetPrimitiveHead() { return primitive_head; }
	Light* GetLightHead() { return light_head; }
	Camera* GetCamera() { return camera; }
	Color GetBackgroundColor() { return background_color; }

	void AddPrimitive( Primitive* );
	void AddLight( Light* );
	void CreateScene( std::string file );
	Collider* FindNearestCollide( Vector3 ray_O , Vector3 ray_V );
	LightCollider* FindNearestLight( Vector3 ray_O , Vector3 ray_V );
};

#endif
