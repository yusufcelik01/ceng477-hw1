SHELL := /bin/bash
CFLAGS =  -g  -O3  # -funroll-loops
CXXFLAGS = -std=c++17


all:
	g++ $(CFLAGS) $(CXXFLAGS) parser.cpp ppm.cpp notth_raytracer.cpp raytracer_math.cpp tinyxml2.cpp -o raytracer 

thread:
	g++ $(CFLAGS) $(CXXFLAGS) parser.cpp ppm.cpp raytracer.cpp raytracer_math.cpp tinyxml2.cpp -o raytracer  -pthread

test: test.cpp
	g++ -O3 raytracer_math.cpp test.cpp -o test -std=c++17 -g
    
png:
	for foo in *.ppm; do convert "$$foo" "$${foo%.ppm}.png"; done

bash:
	for foo in inputs/*.xml; do ./bash_test.sh ; echo "input/$$foo"; done
	
render_small:
	-./raytracer inputs/simple.xml
	-./raytracer inputs/simple_shading.xml
	-./raytracer inputs/simple_reflectance.xml
	-./raytracer inputs/mirror_spheres.xml



render__all:
	for foo in inputs/*xml; do bar="time_$${foo#inputs/}"; (time ./raytracer "$$foo") 2> "$${bar%.xml}.txt" ;  done
