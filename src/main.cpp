#include"raytracer.h"
#include<string>
#include<ctime>
int main() {
	double start_time = clock();
	srand(time(0));
	Raytracer* raytracer = new Raytracer;
	raytracer->SetInput("scene.txt");
	raytracer->SetOutput("picture.bmp");
	raytracer->Run();
	delete raytracer;
	//Bezier b;
	//b.Collide(Vector3(0, -10, 5), Vector3(0, 1, 0));





	double end_time = clock();
	printf("Escaped time: %.5lf\n" , (end_time - start_time)/CLOCKS_PER_SEC);
	return 0;
}
