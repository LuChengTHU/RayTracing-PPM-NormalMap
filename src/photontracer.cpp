#include"photontracer.h"
#include"scene.h"
#include<cstdlib>
#include<cstdio>
#include<thread>
#define ran() (double(rand() * (RAND_MAX + 1) + rand()) / ((RAND_MAX + 1) * (RAND_MAX + 1) - 1))

const int MAX_PHOTONTRACING_DEP = 20;

Photontracer::Photontracer() {
	iteration = 0;
	scene = NULL;
	photonmap = NULL;
	hitpointMap = NULL;
}

bool Photontracer::PhotonDiffusion( Collider* collider , Photon photon , int dep , bool refracted , double* prob ) {
	Primitive* pri = collider->GetPrimitive();
	Material* material = pri->GetMaterial();
	Color color = material->color;
	if ( material->texture != NULL ) color = color * pri->GetTexture(collider->C);
	double eta = material->diff * color.Power();
	if ( eta <= ran() * ( *prob ) ) {
		*prob -= eta;
		return false;
	}

	photon.dir = collider->N.Diffuse();
	photon.power = photon.power * color / color.Power();
	PhotonTracing( photon , dep + 1 , refracted );
	return true;
}

bool Photontracer::PhotonReflection( Collider* collider , Photon photon , int dep , bool refracted , double* prob ) {
	Primitive* pri = collider->GetPrimitive();
	Material* material = pri->GetMaterial();
	Color color = material->color;
	if ( material->texture != NULL ) color = color * pri->GetTexture(collider->C);
	double eta = material->refl * color.Power();

	if ( eta <= ran() * ( *prob ) ) {
		*prob -= material->refl;
		return false;
	}
	
	photon.dir = photon.dir.Reflect( collider->N );
	if ( material->drefl > EPS ) {
		Vector3 Dx = photon.dir.GetAnVerticalVector();
		Vector3 Dy = photon.dir.Cross(Dx);
		Dx = Dx.GetUnitVector() * material->drefl;
		Dy = Dy.GetUnitVector() * material->drefl;
		double x , y;
		do {
			x = ran() * 2 - 1;
			y = ran() * 2 - 1;
		} while ( x * x + y * y > 1 );
		x *= material->drefl;
		y *= material->drefl;
		photon.dir += Dx * x + Dy * y;
	}
	photon.power = photon.power * color / color.Power();
	PhotonTracing( photon , dep + 1 , refracted );
	return true;
}

bool Photontracer::PhotonRefraction( Collider* collider , Photon photon , int dep , bool refracted , double* prob ) {
	Primitive* pri = collider->GetPrimitive();
	Material* material = pri->GetMaterial();
	double eta = material->refr;
	if ( refracted ) {
		Color trans = (material->absor * -collider->dist).Exp();
		eta *= trans.Power();
		photon.power = photon.power * trans / trans.Power();
	}

	if ( eta <= ran() * ( *prob ) ) {
		*prob -= material->refr;
		return false;
	}
	
	double n = material->rindex;
	if ( !refracted ) n = 1 / n;
	bool nextRefracted = refracted;
	photon.dir = photon.dir.Refract( collider->N , n , &nextRefracted );
	PhotonTracing( photon , dep + 1 , nextRefracted );
	return true;
}

void Photontracer::PhotonTracing( Photon photon , int dep , bool refracted ) {
	if ( dep > MAX_PHOTONTRACING_DEP ) return;
	Collider* collider = scene->FindNearestCollide( photon.pos , photon.dir );
	
	if ( collider != NULL ) {
		Primitive* nearest_primitive = collider->GetPrimitive();
		
		photon.pos = collider->C;
		if ( nearest_primitive->GetMaterial()->diff > EPS && dep > 1 ) {
			if (photonmap != NULL)
				photonmap->Store(photon);
			if (hitpointMap != NULL)
				hitpointMap->InsertPhoton(photon);	
		}
		
		double prob = 1;
		if ( PhotonDiffusion( collider , photon , dep , refracted , &prob ) == false )
		if ( PhotonReflection( collider , photon , dep , refracted , &prob ) == false )
		if ( PhotonRefraction( collider	 , photon , dep , refracted , &prob ) == false );
		delete collider;
	}
}

void Photontracer::Emitting(int threadID, int randID) {
	srand(randID);
	double totalPower = 0;
	for ( Light* light = scene->GetLightHead() ; light != NULL ; light = light->GetNext() )
		totalPower += light->GetColor().Power();
	double photonPower = totalPower / scene->GetCamera()->GetEmitPhotons();

	for ( Light* light = scene->GetLightHead() ; light != NULL ; light = light->GetNext() ) {
		int lightPhotons = (int)(light->GetColor().Power() / photonPower);
		lightPhotons = lightPhotons / PM_MAX_THREADS + ((lightPhotons % PM_MAX_THREADS > threadID) ? 1 : 0);

		for (int i = 0; i < lightPhotons ; i++) {
			if (scene->GetCamera()->GetAlgorithm() == "PM" && threadID == 0 && (i & 65535) == 0)
				printf("Emitted Photons= %d, Stored Photons= %d\n", i * PM_MAX_THREADS, photonmap->GetStoredPhotons());
			Photon photon = light->EmitPhoton();
			photon.power *= totalPower;
			PhotonTracing( photon , 1 , false );
		}
	}
	completeThread[threadID] = true;
}

Photonmap* Photontracer::CalnPhotonmap() {
	photonmap = new Photonmap(scene->GetCamera()->GetMaxPhotons());
	photonmap->SetEmitPhotons(scene->GetCamera()->GetEmitPhotons());
	Run();
	photonmap->Balance();
	return photonmap;
}

void Photontracer::Run(int randIDBase) {
	iteration++;
	completeThread = new bool[PM_MAX_THREADS];
	for (int i = 0; i < PM_MAX_THREADS; i++)
		completeThread[i] = false;

	for (int i = 0; i < PM_MAX_THREADS; i++) {
		std::thread subThread(&Photontracer::Emitting, this, i, (randIDBase + iteration) * PM_MAX_THREADS + i);
		if (i == PM_MAX_THREADS - 1)
			subThread.join();
		else
			subThread.detach();
	}

	for (bool end = false; !end; ) {
		end = true;
		for (int i = 0; i < PM_MAX_THREADS; i++)
			if (!completeThread[i]) end = false;
	}
	delete[] completeThread;
}
