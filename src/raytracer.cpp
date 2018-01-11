#include"raytracer.h"
#include<cstdlib>
#include<thread>
#include<cstdio>
#include<mutex>
#include<cmath>
#define ran() (double(rand() * (RAND_MAX + 1) + rand()) / ((RAND_MAX + 1) * (RAND_MAX + 1) - 1))

std::mutex* completeRow;

const int MAX_DREFL_DEP = 2;
const int MAX_RAYTRACING_DEP = 20;
const int HASH_FAC = 7;
const int HASH_MOD = 10000007;

Raytracer::Raytracer() {
	scene = new Scene;
	photonmap = NULL;
	hitpointMap = NULL;
}

Raytracer::~Raytracer() {
	if (scene != NULL)
		delete scene;
}

Color Raytracer::CalnDiffusion( Collider* collider , int* hash, int rc, Color weight) {
	Primitive* pri = collider->GetPrimitive();
	Color color = pri->GetMaterial()->color;
	if (pri->GetMaterial()->texture != NULL)
	{
		if (pri->getName() == 0)
			color = color * pri->GetTexture(collider->C);
		else
			color = color * pri->GetTexture(Vector3(collider->u, collider->v, 0));
	}
		
	Color ret = color * scene->GetBackgroundColor() * pri->GetMaterial()->diff;
	for ( Light* light = scene->GetLightHead() ; light != NULL ; light = light->GetNext() )
		ret += color * light->GetIrradiance( collider , scene->GetPrimitiveHead() , scene->GetCamera()->GetShadeQuality() , hash );
	
	if (camera->GetAlgorithm() == "PM")
		ret += color * photonmap->GetIrradiance( collider , camera->GetSampleDist() , camera->GetSamplePhotons() );
	
	if (camera->GetAlgorithm() == "PPM" || camera->GetAlgorithm() == "SPPM") {
		Hitpoint hitpoint;
		hitpoint.pos = collider->C;
		hitpoint.dir = collider->I;
		hitpoint.N = collider->N;
		hitpoint.primitive = collider->GetPrimitive();
		hitpoint.rc = rc;
		hitpoint.weight = weight * color;
		hitpoint.R2 = photonmap->GetRadius2(collider, camera->GetSampleDist(), camera->GetSamplePhotons());
		hitpointMap->Store(hitpoint);
	}
	
	return ret;
}

Color Raytracer::CalnReflection( Collider* collider , Vector3 ray_V , int dep , bool refracted , int* hash, int rc, Color weight) {
	Primitive* pri = collider->GetPrimitive();
	ray_V = ray_V.Reflect( collider->N );

	if ( pri->GetMaterial()->drefl < EPS || dep > MAX_DREFL_DEP ) {
		Color alpha = pri->GetMaterial()->color * pri->GetMaterial()->refl;
		return RayTracing( collider->C , ray_V , dep + 1 , refracted , hash, rc, weight * alpha) * alpha;
	}

	Vector3 Dx = ray_V.GetAnVerticalVector();
	Vector3 Dy = ray_V.Cross(Dx);
	Dx = Dx.GetUnitVector() * pri->GetMaterial()->drefl;
	Dy = Dy.GetUnitVector() * pri->GetMaterial()->drefl;
	
	int totalSample = camera->GetDreflQuality();
	Color rcol, alpha = pri->GetMaterial()->color * pri->GetMaterial()->refl / totalSample;
	for ( int k = 0 ; k < totalSample ; k++ ) {
		double x , y;
		do {
			x = ran() * 2 - 1;
			y = ran() * 2 - 1;
		} while ( x * x + y * y > 1 );
		x *= pri->GetMaterial()->drefl;
		y *= pri->GetMaterial()->drefl;

		rcol += RayTracing( collider->C , ray_V + Dx * x + Dy * y , dep + MAX_DREFL_DEP , refracted , NULL, rc, weight * alpha);
	}
	return rcol * alpha;
}

Color Raytracer::CalnRefraction( Collider* collider , Vector3 ray_V , int dep , bool refracted , int* hash, int rc, Color weight) {
	Primitive* pri = collider->GetPrimitive();
	double n = pri->GetMaterial()->rindex;
	if ( !refracted ) n = 1 / n;
	
	bool nextRefracted = refracted;
	ray_V = ray_V.Refract( collider->N , n , &nextRefracted );
	
	Color alpha = Color(1, 1, 1) * pri->GetMaterial()->refr;
	if (refracted)
		alpha *= (pri->GetMaterial()->absor * -collider->dist).Exp();
	Color rcol = RayTracing( collider->C , ray_V , dep + 1 , nextRefracted , hash, rc, weight * alpha);
	return rcol * alpha;
}

Color Raytracer::RayTracing( Vector3 ray_O , Vector3 ray_V , int dep , bool refracted , int* hash, int rc, Color weight) {
	if ( dep > MAX_RAYTRACING_DEP ) return Color();
	if ( hash != NULL ) *hash = ( *hash * HASH_FAC ) % HASH_MOD;

	Color ret;
	Collider* collider = scene->FindNearestCollide(ray_O, ray_V);
	LightCollider* lightCollider = scene->FindNearestLight(ray_O, ray_V);
	
	if (lightCollider != NULL) {
		Light* nearest_light = lightCollider->GetLight();
		if (collider == NULL || lightCollider->dist < collider->dist) {
			if ( hash != NULL ) *hash = ( *hash + nearest_light->GetSample() ) % HASH_MOD;
			ret += nearest_light->GetColor() / nearest_light->GetColor().RGBMax();
		}
		delete lightCollider;
	}
	
	if ( collider != NULL ) {
		Primitive* nearest_primitive = collider->GetPrimitive();		
		if ( hash != NULL ) *hash = ( *hash + nearest_primitive->GetSample() ) % HASH_MOD;
		if ( nearest_primitive->GetMaterial()->diff > EPS ) ret += CalnDiffusion( collider , hash, rc, weight);
		if (camera->GetAlgorithm() != "RC") {
			if ( nearest_primitive->GetMaterial()->refl > EPS ) ret += CalnReflection( collider , ray_V , dep , refracted , hash, rc, weight);
			if ( nearest_primitive->GetMaterial()->refr > EPS ) ret += CalnRefraction( collider , ray_V , dep , refracted , hash, rc, weight);
		}
		delete collider;
	}

	if ( dep == 1 ) ret = ret.Confine();
	return ret;
}

void Raytracer::Sampling(int threadID, int randID) {
	srand(randID);
	Vector3 ray_O = camera->GetO();

	for ( int i = 0 ; i < H ; i++ ) {
		if (!completeRow[i].try_lock()) continue;
		for ( int j = 0 ; j < W ; j++ ) {
			if (camera->GetAperture() < EPS) {
				Vector3 ray_V = camera->Emit( i , j );
				sample[i][j] = 0;
				Color color = camera->GetColor(i, j);
				color += RayTracing(ray_O, ray_V, 1, false, &sample[i][j], i * W + j, Color(1, 1, 1));
				camera->SetColor(i, j, color.Confine());
			} else {
				int dofSample = camera->GetDofSample();
				int iteration = (camera->GetAlgorithm() == "SPPM") ? 1 : dofSample;
				Color color = camera->GetColor(i, j);
				for (int k = 0; k < iteration; k++) {
					Vector3 dof_O, dof_V;
					camera->DofEmit(i, j, &dof_O, &dof_V);
					color += RayTracing(dof_O, dof_V, 1, false, NULL, i * W + j, Color(1, 1, 1) / dofSample) / dofSample;
				}
				camera->SetColor(i, j, color.Confine());
			}

			
			if (j == W - 1) {
				printf("Sampling=%d/%d\n", i, H);
				if (((i & 7) == 0) && (camera->GetAlgorithm() == "PPM" || camera->GetAlgorithm() == "SPPM"))
					printf("Stored hitpoints=%d\n", hitpointMap->GetStoredHitpoints());
			}
		}
	}
	completeThread[threadID] = true;
}

void Raytracer::Resampling(int threadID, int randID) {
	srand(randID);
	Vector3 ray_O = camera->GetO();
	
	for ( int i = 0 ; i < H ; i++ ) {
		if (!completeRow[i].try_lock()) continue;
		for ( int j = 0 ; j < W ; j++ ) {
			if (!((i == 0 || sample[i][j] == sample[i - 1][j])&&(i == H - 1 || sample[i][j] == sample[i + 1][j])&&(j == 0 || sample[i][j] == sample[i][j - 1])&&(j == W - 1 || sample[i][j] == sample[i][j + 1]))) {
				Color color = camera->GetColor(i, j) / 5;
				for ( int r = -1 ; r <= 1 ; r++ )
					for ( int c = -1 ; c <= 1 ; c++ ) {
						if (((r + c) & 1) == 0) continue;
						Vector3 ray_V = camera->Emit( i + ( double ) r / 3 , j + ( double ) c / 3 );
						color += RayTracing(ray_O, ray_V, 1, false, NULL, i * W + j, Color(1, 1, 1) / 5) / 5;
					}
				camera->SetColor( i , j , color.Confine() );
			}

			if (j == W - 1) {
				printf("Resampling=%d/%d\n", i, H);
				if (((i & 7) == 0) && (camera->GetAlgorithm() == "PPM" || camera->GetAlgorithm() == "SPPM"))
					printf("Stored hitpoints=%d\n", hitpointMap->GetStoredHitpoints());
			}
		}
	}
	completeThread[threadID] = true;
}

void Raytracer::MultiThreadSampling(int randIDBase) {
	completeThread = new bool[RT_MAX_THREADS];
	for (int i = 0; i < RT_MAX_THREADS; i++)
		completeThread[i] = false;
	
	completeRow = new std::mutex[H];
	for (int i = 0; i < RT_MAX_THREADS; i++) {
		std::thread subThread(&Raytracer::Sampling, this, i, randIDBase + i);
		if (i == RT_MAX_THREADS - 1)
			subThread.join();
		else
			subThread.detach();
	}
	delete[] completeRow;
	
	for (bool end = false; !end; ) {
		end = true;
		for (int i = 0; i < RT_MAX_THREADS; i++)
			if (!completeThread[i]) end = false;
	}
	delete[] completeThread;
}

void Raytracer::MultiThreadResampling(int randIDBase) {
	completeThread = new bool[RT_MAX_THREADS];
	for (int i = 0; i < RT_MAX_THREADS; i++)
		completeThread[i] = false;
	
	completeRow = new std::mutex[H];
	for (int i = 0; i < RT_MAX_THREADS; i++) {
		std::thread subThread(&Raytracer::Resampling, this, i, randIDBase + RT_MAX_THREADS + i);
		if (i == RT_MAX_THREADS - 1)
			subThread.join();
		else
			subThread.detach();
	}
	delete[] completeRow;
	
	for (bool end = false; !end; ) {
		end = true;
		for (int i = 0; i < RT_MAX_THREADS; i++)
			if (!completeThread[i]) end = false;
	}
	delete[] completeThread;
}

void Raytracer::ProgressivePhotonMapping(int SPPMIter) {
	if (camera->GetAlgorithm() != "SPPM")
		GenerateImage("picture_RT.bmp");
	
	int storedHitpoints = hitpointMap->GetStoredHitpoints();
	printf("Stored Hitpoints: %d\n", storedHitpoints);
	printf("HitpointMap Balancing...\n");
	hitpointMap->Balance();
	
	Color** rtColor = new Color*[H];
	for (int i = 0; i < H; i++) {
		rtColor[i] = new Color[W];
		for (int j = 0; j < W; j++)
			rtColor[i][j] = camera->GetColor(i, j);
	}
	
	Photontracer* photontracer = new Photontracer;
	photontracer->SetScene( scene );
	photontracer->SetHitpointMap(hitpointMap);
	for (int iter = 1; iter <= camera->GetIterations(); iter++) {
		photontracer->Run(SPPMIter * camera->GetIterations());
		Hitpoint* hitpoints = hitpointMap->GetHitpoints();
		
		hitpointMap->MaintainHitpoints();
		for (int r = 0; r < H; r++)
			for (int c = 0; c < W; c++)
				camera->SetColor(r, c, rtColor[r][c]);
		
		double minR2 = camera->GetSampleDist(), maxR2 = 0;
		double minNum = 1000000000, maxNum = 0;
		for (int i = 1; i <= storedHitpoints; i++) {
			int r = hitpoints[i].rc / W;
			int c = hitpoints[i].rc % W;
			Color color = camera->GetColor(r, c) + hitpoints[i].color * hitpoints[i].weight * (4.0 / (hitpoints[i].R2 * camera->GetEmitPhotons() * iter));
			camera->SetColor(r, c, color.Confine());
			
			if (hitpoints[i].R2 < minR2) minR2 = hitpoints[i].R2;
			if (hitpoints[i].R2 > maxR2) maxR2 = hitpoints[i].R2;
			if (hitpoints[i].num < minNum) minNum = hitpoints[i].num;
			if (hitpoints[i].num > maxNum) maxNum = hitpoints[i].num;
		}
		printf("Iter=%d, Num=%.2lf~%.2lf, Radius=%.6lf~%.6lf\n", iter, minNum, maxNum, sqrt(minR2), sqrt(maxR2));
		
		if (iter % 10 == 0)
			GenerateImage(output);
	}
	
	if (camera->GetAlgorithm() == "SPPM")
		GenerateImage(output);
	
	delete photontracer;
	for (int i = 0; i < H; i++)
		delete[] rtColor[i];
	delete[] rtColor;
}

void Raytracer::GenerateImage(std::string file) {
	Bmp* bmp = new Bmp;
	camera->Output(bmp);
	bmp->Output(file);
	delete bmp;
}

void Raytracer::Run() {
	scene->CreateScene( input );
	camera = scene->GetCamera();
	
	int SPPMIteration = (camera->GetAlgorithm() == "SPPM") ? camera->GetDofSample() : 1;
	for (int iter = 0; iter < SPPMIteration; iter++) {
		if (camera->GetAlgorithm() == "SPPM")
			printf("SPPM Iteration= %d\n", iter);
		
		if (camera->GetAlgorithm() == "PPM" || camera->GetAlgorithm() == "SPPM") {
			hitpointMap = new HitpointMap(camera->GetMaxHitpoints());
			hitpointMap->SetReduction(camera->GetReduction());
		}
		
		if (camera->GetAlgorithm() == "PM" || camera->GetAlgorithm() == "PPM" || camera->GetAlgorithm() == "SPPM") {
			Photontracer* photontracer = new Photontracer;
			photontracer->SetScene( scene );
			photonmap = photontracer->CalnPhotonmap();
			delete photontracer;
		}
		
		if (camera->GetAlgorithm() == "PM")
			printf("Stored Photons= %d\n", photonmap->GetStoredPhotons());
		
		H = camera->GetH();
		W = camera->GetW();

		sample = new int*[H];
		for ( int i = 0 ; i < H ; i++ )
			sample[i] = new int[W];
		
		MultiThreadSampling(2 * iter * RT_MAX_THREADS);
		if (camera->GetAperture() < EPS) {
			if (camera->GetAlgorithm() == "PPM" || camera->GetAlgorithm() == "SPPM") {
				int storedHitpoints = hitpointMap->GetStoredHitpoints();
				Hitpoint* hitpoints = hitpointMap->GetHitpoints();
				for (int i = 1; i <= storedHitpoints; i++) {
					int r = hitpoints[i].rc / W;
					int c = hitpoints[i].rc % W;
					if (!((r == 0 || sample[r][c] == sample[r - 1][c])&&(r == H - 1 || sample[r][c] == sample[r + 1][c])&&(c == 0 || sample[r][c] == sample[r][c - 1])&&(c == W - 1 || sample[r][c] == sample[r][c + 1])))
						hitpoints[i].weight /= 5;
				}
			}
			MultiThreadResampling((2 * iter + 1) * RT_MAX_THREADS);
		}
		
		for ( int i = 0 ; i < H ; i++ )
			delete[] sample[i];
		delete[] sample;
		
		GenerateImage(output);
		
		if (camera->GetAlgorithm() == "PPM" || camera->GetAlgorithm() == "SPPM")
			ProgressivePhotonMapping(iter);
		
		if (hitpointMap != NULL) {
			delete hitpointMap;
			hitpointMap = NULL;
		}
		if (photonmap != NULL) {
			delete photonmap;
			photonmap = NULL;
		}
	}
}
