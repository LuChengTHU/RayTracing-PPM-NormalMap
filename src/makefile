objects = vector3.o color.o bmp.o material.o \
	scene.o \
	camera.o \
	primitive.o plane.o sphere.o Bezier.o NomalBall.o NomalPlane.o\
	light.o areaLight.o \
	photonmap.o hitpointMap.o\
	raytracer.o photontracer.o

all: main.exe clean

main.exe: main.cpp $(objects)
	g++ main.cpp $(objects) -o main.exe -O3

vector3.o: vector3.h vector3.cpp
	g++ -c vector3.cpp -O3

color.o: color.h color.cpp
	g++ -c color.cpp -O3

bmp.o: bmp.h bmp.cpp
	g++ -c bmp.cpp -O3

material.o: material.h material.cpp
	g++ -c material.cpp -O3

scene.o: scene.h scene.cpp
	g++ -c scene.cpp -O3

camera.o: camera.h camera.cpp
	g++ -c camera.cpp -O3

primitive.o: primitive.h primitive.cpp
	g++ -c primitive.cpp -O3

plane.o: plane.h plane.cpp
	g++ -c plane.cpp -O3
	
sphere.o: sphere.h sphere.cpp
	g++ -c sphere.cpp -O3

Bezier.o: Bezier.h Bezier.cpp
	g++ -c Bezier.cpp -O3

NomalBall.o: nomalBall.h nomalBall.cpp
	g++ -c nomalBall.cpp -O3 -std=c++11

NomalPlane.o: nomalplane.h nomalplane.cpp
	g++ -c nomalplane.cpp -O3 -std=c++11

light.o: light.h light.cpp
	g++ -c light.cpp -O3

areaLight.o: areaLight.h areaLight.cpp
	g++ -c areaLight.cpp -O3

photonmap.o: photonmap.h photonmap.cpp
	g++ -c photonmap.cpp -O3

hitpointMap.o: hitpointMap.h hitpointMap.cpp
	g++ -c hitpointMap.cpp -O3 -std=c++11
	
raytracer.o: raytracer.h raytracer.cpp
	g++ -c raytracer.cpp -O3 -std=c++11

photontracer.o: photontracer.h photontracer.cpp
	g++ -c photontracer.cpp -O3 -std=c++11

clean:
	del *.bak *.o
